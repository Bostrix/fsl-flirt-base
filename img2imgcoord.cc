/*  img2imgcoord.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2000 University of Oxford  */

/*  CCOPYRIGHT  */

#include <string>
#include <iostream>
#include <fstream>
#include <unistd.h>
#define WANT_STREAM
#define WANT_MATH

#include "newmatap.h"
#include "newmatio.h"
#include "mjimage.h"
#include "miscmaths.h"
#include "miscimfns.h"
#include "generalio.h"

#ifndef NO_NAMESPACE
 using namespace MISCMATHS;
 using namespace MISCIMFNS;
 using namespace MJIMAGE;
 using namespace NEWMAT;
 using namespace GENERALIO;
#endif

////////////////////////////////////////////////////////////////////////////
// the real defaults are provided in the function parse_command_line

class globaloptions {
public:
  string destfname;
  string srcfname;
  string xfmfname;
  string coordfname;
  bool mm;
  int verbose;
public:
  globaloptions();
  ~globaloptions() {};
};

globaloptions globalopts;


globaloptions::globaloptions()
{
  // set up defaults
  destfname = "";
  srcfname = "";
  coordfname = "";
  xfmfname = "";
  mm = false;
}


////////////////////////////////////////////////////////////////////////////

// Parsing functions for command line parameters

void print_usage(int argc, char *argv[])
{
  cout << "Usage: " << argv[0] << " [options] <filename containing coordinates>\n\n"
       << "  Options are:\n"
       << "        -src <Source Image filename>         (NB: volume, not timeseries)\n"
       << "        -dest <Destination Image filename>   (NB: volume, not timeseries)\n"
       << "        -xfm <Source to Destination transform filename>\n"
       << "        -vox                                 (all coordinates in voxels - default)\n"
       << "        -mm                                  (all coordinates in mm)\n"
       << "        -help\n\n"
       << " Note that the first three options are compulsory\n";
}


void parse_command_line(int argc, char* argv[])
{
  if(argc<2){
    print_usage(argc,argv);
    exit(1);
  }


  int n=1;
  string arg;
  char first;

  while (n<argc) {
    arg=argv[n];
    if (arg.size()<1) { n++; continue; }
    first = arg[0];
    if (first!='-') {
      globalopts.coordfname = arg;
      n++;
      continue;
    }
    
    // put options without arguments here
    if ( arg == "-help" ) {
      print_usage(argc,argv);
      exit(0);
    } else if ( arg == "-vox" ) {
      globalopts.mm = false;
      n++;
      continue;
    } else if ( arg == "-mm" ) {
      globalopts.mm = true;
      n++;
      continue;
    } else if ( arg == "-v" ) {
      globalopts.verbose = 5;
      n++;
      continue;
    }

    if (n+1>=argc) 
      { 
	cerr << "Lacking argument to option " << arg << endl;
	break; 
      }

    // put options with 1 argument here
    if ( arg == "-dest") {
      globalopts.destfname = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-src") {
      globalopts.srcfname = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-xfm") {
      globalopts.xfmfname = argv[n+1];
      n+=2;
      continue;
    } else { 
      cerr << "Unrecognised option " << arg << endl;
      exit(-1);
    } 

  }  // while (n<argc)

//    if (globalopts.coordfname.size()<1) {
//      cerr << "Input coordinate file not found\n\n";
//      print_usage(argc,argv);
//      exit(2);
//    }
  if ((globalopts.srcfname.size()<1)) {
    cerr << "ERROR:: Source volume filename not found\n\n";
  }
  if ((globalopts.destfname.size()<1)) {
    cerr << "ERROR:: Destination volume filename not found\n\n";
  }
}

////////////////////////////////////////////////////////////////////////////

void print_info(const volume& vol, const string& name) {
  cout << name << ":: SIZE = " << vol.xsize() << " x " << vol.ysize() 
       << " x " << vol.zsize() << endl;
  cout << name << ":: DIMS = " << vol.getx() << " x " << vol.gety() 
       << " x " << vol.getz() << " mm" << endl << endl;
}  

////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[])
{
  parse_command_line(argc,argv);


  volume srcvol, destvol;
    // read volumes
  if (read_volume_hdr_only(srcvol,globalopts.srcfname)<0) {
    cerr << "Cannot read Source volume" << endl;
    return -1;
  }
  if (read_volume_hdr_only(destvol,globalopts.destfname)<0) {
    cerr << "Cannot read Destination volume" << endl;
    return -1;
  }
    
  if (globalopts.verbose>3) {
    print_info(destvol,"Destination Volume");
    cout << " origin = " << destvol.avw_origin.t() << endl << endl;
    print_info(srcvol,"Source Volume");
    cout << " origin = " << srcvol.avw_origin.t() << endl;
  }

  // read matrices
  Matrix affmat(4,4);
  int returnval;
  returnval = read_matrix(affmat,globalopts.xfmfname,srcvol,destvol);
  if (returnval<0) {
    cerr << "Cannot read transform file" << endl;
    return -2;
  }
    
  if (globalopts.verbose>3) {
    cout << " affmat =" << endl << affmat << endl << endl;
  }

  // Let Volume 2 be Source and Volume 1 be Destination
  //  notate variables as (v=vox, w=world, f=flirt, m=medx, t=dest)
  Matrix vf2w2(4,4), vf1w1(4,4);
  Identity(vf2w2);
  vf2w2(1,1) = srcvol.getx();
  vf2w2(2,2) = srcvol.gety();
  vf2w2(3,3) = srcvol.getz();
  Identity(vf1w1);
  vf1w1(1,1) = destvol.getx();
  vf1w1(2,2) = destvol.gety();
  vf1w1(3,3) = destvol.getz();

  // the swap matrices convert flirt voxels to medx voxels
  Matrix swapy1(4,4), swapy2(4,4);
  Identity(swapy1);  Identity(swapy2);
  swapy1(2,2) = -1.0;
  swapy2(2,2) = -1.0;
  swapy1(2,4) = destvol.ysize()-1.0;
  swapy2(2,4) = srcvol.ysize()-1.0;

  Matrix destvox2world, srcvox2world;
  destvox2world = vf1w1 * swapy1;
  srcvox2world = vf2w2 * swapy2;
  if (globalopts.verbose>3) {
    cout << " destvox2world =" << endl << destvox2world << endl << endl;
    cout << " srcvox2world =" << endl << srcvox2world << endl;
  }

  ColumnVector srccoord(4), destcoord(4), oldsrc(4);
  srccoord = 0;
  destcoord = 0;
  srccoord(4)=1;
  destcoord(4)=1;
  oldsrc = 0;  // 4th component set to 0, so that initially oldsrc -ne srccoord

  cout << "Coordinates in Destination volume (in mm):" << endl;

  if (globalopts.coordfname.size()>1) {
    ifstream matfile(globalopts.coordfname.c_str());
    if (!matfile) { 
      cerr << "Could not open matrix file " << globalopts.coordfname << endl;
      return -1;
    }
    
    while (!matfile.eof()) {
      for (int j=1; j<=3; j++) {
	matfile >> destcoord(j);
      }
      if (globalopts.mm) {  // in mm
	destcoord = vf1w1 * destvox2world.i() * affmat * srcvox2world * vf2w2.i() * srccoord;
      } else { // in voxels
	destcoord = destvox2world.i() * affmat * srcvox2world * srccoord; 
      }
      cout << destcoord(1) << "  " << destcoord(2) << "  " << destcoord(3) << endl;
    }
    
    matfile.close();
  } else {
    cout << "Please type in Source coordinates";
    if (globalopts.mm) { 
      cout << " (in mm) :" << endl;
    } else { 
      cout << " (in voxels) :" << endl; 
    }
    while (!cin.eof()) {
      for (int j=1; j<=3; j++) {
	cin >> srccoord(j);
      }
      if (oldsrc == srccoord)  return 0;
      oldsrc = srccoord;
      if (globalopts.mm) {  // in mm
	destcoord = vf1w1 * destvox2world.i() * affmat * srcvox2world * vf2w2.i() * srccoord;
      } else { // in voxels
	destcoord = destvox2world.i() * affmat * srcvox2world * srccoord; 
      }
      cout << destcoord(1) << "  " << destcoord(2) << "  " << destcoord(3) << endl;
    }
  }


  return 0;
}








