/*  tal2epicoord.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

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
  string talfname;
  string epifname;
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
  talfname = "";
  epifname = "";
  coordfname = "";
  xfmfname = "";
  mm = true;
}


////////////////////////////////////////////////////////////////////////////

// Parsing functions for command line parameters

void print_usage(int argc, char *argv[])
{
  cout << "Usage: " << argv[0] << " [options] <filename containing coordinates>\n\n"
       << "  Options are:\n"
       << "        -tal <Talairach volume filename>\n"
       << "        -epi <example EPI volume filename>   (NB: volume, not timeseries)\n"
       << "        -xfm <EPI to Talairach transform filename>\n"
       << "        -mm                                  (outputs coordinates in mm - default)\n"
       << "        -vox                                 (outputs coordinates in voxels)\n"
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
    if ( arg == "-tal") {
      globalopts.talfname = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-epi") {
      globalopts.epifname = argv[n+1];
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

  if (globalopts.coordfname.size()<1) {
    cerr << "Input coordinate file not found\n\n";
    print_usage(argc,argv);
    exit(2);
  }
  if ((globalopts.epifname.size()<1)) {
    cerr << "ERROR:: EPI volume filename not found\n\n";
  }
  if ((globalopts.talfname.size()<1)) {
    cerr << "ERROR:: Talairach volume filename not found\n\n";
  }
}

////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[])
{
  parse_command_line(argc,argv);


  volume epivol, talvol;
    // read volumes
  if (read_volume_hdr_only(epivol,globalopts.epifname)<0) {
    cerr << "Cannot read EPI volume" << endl;
    return -1;
  }
  if (read_volume_hdr_only(talvol,globalopts.talfname)<0) {
    cerr << "Cannot read Talairach volume" << endl;
    return -1;
  }
    
  if (globalopts.verbose>3) {
    print_volume_info(talvol,"Talairach Volume");
    cout << " origin = " << talvol.avw_origin.t() << endl << endl;
    print_volume_info(epivol,"EPI Volume");
    cout << " origin = " << epivol.avw_origin.t() << endl;
  }

  // read matrices
  Matrix affmat(4,4);
  int returnval;
  returnval = read_matrix(affmat,globalopts.xfmfname,epivol,talvol);
  if (returnval<0) {
    cerr << "Cannot read transform file" << endl;
    return -2;
  }
    
  if (globalopts.verbose>3) {
    cout << " affmat =" << endl << affmat << endl << endl;
  }

  // apply inverse 
  affmat = affmat.i();

  if (globalopts.verbose>3) {
    cout << " Inverse affmat =" << endl << affmat << endl << endl;
  }

  // Let Volume 2 be EPI and Volume 1 be Talairach
  //  notate variables as (v=vox, w=world, f=flirt, m=medx, t=tal)
  Matrix vf2w2(4,4), vf1w1(4,4), vt1vm1(4,4);
  Identity(vf2w2);
  vf2w2(1,1) = epivol.getx();
  vf2w2(2,2) = epivol.gety();
  vf2w2(3,3) = epivol.getz();
  Identity(vf1w1);
  vf1w1(1,1) = talvol.getx();
  vf1w1(2,2) = talvol.gety();
  vf1w1(3,3) = talvol.getz();
  get_outputusermat(globalopts.talfname,vt1vm1);
  vt1vm1 = vt1vm1.i();

  // the swap matrices convert flirt voxels to medx voxels
  Matrix swapy1(4,4), swapy2(4,4);
  Identity(swapy1);  Identity(swapy2);
  swapy1(2,2) = -1.0;
  swapy2(2,2) = -1.0;
  swapy1(2,4) = talvol.ysize()-1.0;
  swapy2(2,4) = epivol.ysize()-1.0;

  Matrix talvox2world, epivox2world;
  talvox2world = vf1w1 * swapy1 * vt1vm1;
  epivox2world = vf2w2 * swapy2;
  if (globalopts.verbose>3) {
    cout << " talvox2world =" << endl << talvox2world << endl << endl;
    cout << " epivox2world =" << endl << epivox2world << endl;
  }

  ColumnVector epicoord(4), talcoord(4);
  epicoord = 0;
  talcoord = 0;
  epicoord(4)=1;
  talcoord(4)=1;

  if (globalopts.mm) {
    cout << "Coordinates in EPI volume (in mm):" << endl;
  } else {
    cout << "Coordinates in EPI volume (in voxels):" << endl;
  }

  ifstream matfile(globalopts.coordfname.c_str());
  if (!matfile) { 
    cerr << "Could not open matrix file " << globalopts.coordfname << endl;
    return -1;
  }

  while (!matfile.eof()) {
    for (int j=1; j<=3; j++) {
      matfile >> talcoord(j);
    }
    if (globalopts.mm) {  // in mm
      epicoord = vf2w2 * epivox2world.i() * affmat * talvox2world * talcoord;
    } else { // in voxels
      epicoord = epivox2world.i() * affmat * talvox2world * talcoord; 
    }
    cout << epicoord(1) << "  " << epicoord(2) << "  " << epicoord(3) << endl;
  }

  matfile.close();

  return 0;
}








