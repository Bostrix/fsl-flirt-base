/*  img2talcoord.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2000 University of Oxford  */

/*  CCOPYRIGHT  */

#include <string>
#include <iostream>
#include <fstream>
#include <unistd.h>
#define WANT_STREAM
#define WANT_MATH

#include "newmat/newmatap.h"
#include "newmat/newmatio.h"
#include "miscmaths/miscmaths.h"
#include "newimage/newimageall.h"

#ifndef NO_NAMESPACE
 using namespace MISCMATHS;
 using namespace NEWMAT;
 using namespace NEWIMAGE;
#endif

////////////////////////////////////////////////////////////////////////////
// the real defaults are provided in the function parse_command_line

class globaloptions {
public:
  string talfname;
  string imgfname;
  string xfmfname;
  string coordfname;
  bool mm;
  int verbose;
  bool medx;
public:
  globaloptions();
  ~globaloptions() {};
};

globaloptions globalopts;


globaloptions::globaloptions()
{
  // set up defaults
  talfname = "";
  imgfname = "";
  coordfname = "";
  xfmfname = "";
  mm = false;
  medx=true;
}


////////////////////////////////////////////////////////////////////////////

// Parsing functions for command line parameters

void print_usage(int argc, char *argv[])
{
  cout << "Usage: " << argv[0] << " [options] <filename containing coordinates>\n\n"
       << "  Options are:\n"
       << "        -tal <Talairach volume filename>\n"
       << "        -img <example IMG volume filename>   (NB: volume, not timeseries)\n"
       << "        -xfm <IMG to Talairach transform filename>\n"
       << "        -vox                                 (input coordinates in voxels - default)\n"
       << "        -mm                                  (input coordinates in mm)\n"
       << "        -flirt                               (use flirt, not medx, coordinate conventions)\n"    
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
    } else if ( arg == "-flirt" ) {
      globalopts.medx = false;
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
    if ( arg == "-tal") {
      globalopts.talfname = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-img") {
      globalopts.imgfname = argv[n+1];
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
  if ((globalopts.imgfname.size()<1)) {
    cerr << "ERROR:: IMG volume filename not found\n\n";
  }
  if ((globalopts.talfname.size()<1)) {
    cerr << "ERROR:: Talairach volume filename not found\n\n";
  }
}

////////////////////////////////////////////////////////////////////////////

void print_info(const volume<float>& vol, const string& name) {
  cout << name << ":: SIZE = " << vol.xsize() << " x " << vol.ysize() 
       << " x " << vol.zsize() << endl;
  cout << name << ":: DIMS = " << vol.xdim() << " x " << vol.ydim() 
       << " x " << vol.zdim() << " mm" << endl << endl;
}  

////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[])
{
  parse_command_line(argc,argv);


  volume<float> imgvol, talvol;
    // read volumes
  if (read_volume_hdr_only(imgvol,globalopts.imgfname)<0) {
    cerr << "Cannot read IMG volume" << endl;
    return -1;
  }
  if (read_volume_hdr_only(talvol,globalopts.talfname)<0) {
    cerr << "Cannot read Talairach volume" << endl;
    return -1;
  }
    
  if (globalopts.verbose>3) {
    print_info(talvol,"Talairach Volume");
    cout << " origin = " << talvol.getorigin().t() << endl << endl;
    print_info(imgvol,"IMG Volume");
    cout << " origin = " << imgvol.getorigin().t() << endl;
  }

  // read matrices
  Matrix affmat(4,4);
  int returnval;
  returnval = read_matrix(affmat,globalopts.xfmfname,imgvol,talvol);
  if (returnval<0) {
    cerr << "Cannot read transform file" << endl;
    return -2;
  }
    
  if (globalopts.verbose>3) {
    cout << " affmat =" << endl << affmat << endl << endl;
  }

  // Let Volume 2 be IMG and Volume 1 be Talairach
  //  notate variables as (v=vox, w=world, f=flirt, m=medx, t=tal)
  Matrix vf2w2(4,4), vf1w1(4,4), vt1vm1(4,4);
  Identity(vf2w2);
  vf2w2(1,1) = imgvol.xdim();
  vf2w2(2,2) = imgvol.ydim();
  vf2w2(3,3) = imgvol.zdim();
  Identity(vf1w1);
  vf1w1(1,1) = talvol.xdim();
  vf1w1(2,2) = talvol.ydim();
  vf1w1(3,3) = talvol.zdim();
  get_outputusermat(globalopts.talfname,vt1vm1);
  vt1vm1 = vt1vm1.i();

  // the swap matrices convert flirt voxels to medx voxels
  Matrix swapy1(4,4), swapy2(4,4);
  Identity(swapy1);  Identity(swapy2);
  if (globalopts.medx) {
    swapy1(2,2) = -1.0;
    swapy2(2,2) = -1.0;
    swapy1(2,4) = talvol.ysize()-1.0;
    swapy2(2,4) = imgvol.ysize()-1.0;
  }

  Matrix talvox2world, imgvox2world;
  talvox2world = vf1w1 * swapy1 * vt1vm1;
  imgvox2world = vf2w2 * swapy2;
  if (globalopts.verbose>3) {
    cout << " talvox2world =" << endl << talvox2world << endl << endl;
    cout << " imgvox2world =" << endl << imgvox2world << endl;
  }

  ColumnVector imgcoord(4), talcoord(4), oldimg(4);
  imgcoord = 0;
  talcoord = 0;
  imgcoord(4)=1;
  talcoord(4)=1;
  oldimg = 0;  // 4th component set to 0, so that initially oldimg -ne imgcoord

  cout << "Coordinates in Talairach volume (in mm):" << endl;

  if (globalopts.coordfname.size()>1) {
    ifstream matfile(globalopts.coordfname.c_str());
    if (!matfile) { 
      cerr << "Could not open matrix file " << globalopts.coordfname << endl;
      return -1;
    }
    
    while (!matfile.eof()) {
      for (int j=1; j<=3; j++) {
	matfile >> imgcoord(j);
      }
      if (globalopts.mm) {  // in mm
	talcoord = talvox2world.i() * affmat * imgvox2world * vf2w2.i() * imgcoord;
      } else { // in voxels
	talcoord = talvox2world.i() * affmat * imgvox2world * imgcoord; 
      }
      cout << talcoord(1) << "  " << talcoord(2) << "  " << talcoord(3) << endl;
    }
    
    matfile.close();
  } else {
    cout << "Please type in IMG coordinates";
    if (globalopts.mm) { 
      cout << " (in mm) :" << endl;
    } else { 
      cout << " (in voxels) :" << endl; 
    }
    while (!cin.eof()) {
      for (int j=1; j<=3; j++) {
	cin >> imgcoord(j);
      }
      if (oldimg == imgcoord)  return 0;
      oldimg = imgcoord;
      if (globalopts.mm) {  // in mm
	talcoord = talvox2world.i() * affmat * imgvox2world * vf2w2.i() * imgcoord;
      } else { // in voxels
	talcoord = talvox2world.i() * affmat * imgvox2world * imgcoord; 
      }
      cout << talcoord(1) << "  " << talcoord(2) << "  " << talcoord(3) << endl;
    }
  }


  return 0;
}








