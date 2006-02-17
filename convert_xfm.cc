/*  convert_xfm.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 1999-2006 University of Oxford  */

/*  CCOPYRIGHT  */

#include <string>
#include <iostream>
#include <fstream>
#include <unistd.h>
#define WANT_STREAM
#define WANT_MATH

string version="2.0"

#include "newmatap.h"
#include "newmatio.h"
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"

 using namespace MISCMATHS;
 using namespace NEWMAT;
 using namespace NEWIMAGE;

////////////////////////////////////////////////////////////////////////////
// the real defaults are provided in the function parse_command_line

class globaloptions {
public:
  string testfname;
  string reffname;
  string outputmatascii;
  string initmatfname;
  string concatfname;
  string fixfname;
  string intervolfname;
  string xfm_type;
  int verbose;
  bool inverse;
  bool matonly;
public:
  globaloptions();
  ~globaloptions() {};
};

globaloptions globalopts;


globaloptions::globaloptions()
{
  // set up defaults
  reffname = "";

  testfname = "";
  outputmatascii = "";
  initmatfname = "";
  concatfname = "";
  fixfname = "";
  intervolfname = "";
  xfm_type = "a";
  inverse = false;
  matonly = true;
}


////////////////////////////////////////////////////////////////////////////

// Parsing functions for command line parameters

void print_usage(int argc, char *argv[])
{
  cout << "convert_xfm (Version " << version << ")" << endl
       << "Tool for manipulating FSL transformation matrices" << endl
       << "Copyright(c) 1999, University of Oxford (Mark Jenkinson)" << endl
       << endl 
       << "Usage: " << argv[0] << " [options] <input-matrix-filename>" << endl
       << "  e.g. " << argv[0] << " -omat <outmat> -inverse <inmat>" << endl
       << "       " << argv[0] << " -omat <outmat_AtoC> -concat <mat_BtoC> <mat_AtoB>" << endl << endl
       << "  Available options are:" << endl
       << "        -omat <matrix-filename>            (4x4 ascii format)" << endl
    //       << "        -ominc <matrix-filename>           (MINC format)" << endl
       << "        -concat <second-matrix-filename>" << endl
       << "        -fixscaleskew <second-matrix-filename>" << endl
       << "        -inverse                           (Reference image must be the one originally used)" << endl
       << "        -help" << endl;
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
      globalopts.initmatfname = arg;
      n++;
      continue;
    }
    
    // put options without arguments here
    if ( arg == "-help" ) {
      print_usage(argc,argv);
      exit(0);
    } else if ( arg == "-inverse" ) {
      globalopts.inverse = true;
      n++;
      continue;
    } else if ( arg == "-matonly" ) {
      globalopts.matonly = true;
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
    if ( arg == "-ref") {
      globalopts.reffname = argv[n+1];
      globalopts.matonly = false;
      n+=2;
      continue;
    } else if ( arg == "-in") {
      globalopts.testfname = argv[n+1];
      globalopts.matonly = false;
      n+=2;
      continue;
    } else if ( arg == "-concat") {
      globalopts.concatfname = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-fixscaleskew") {
      globalopts.fixfname = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-middlevol") {
      globalopts.intervolfname = argv[n+1];
      globalopts.matonly = false;
      n+=2;
      continue;
    } else if ( arg == "-omat") {
      globalopts.outputmatascii = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-verbose") {
      globalopts.verbose = atoi(argv[n+1]);
      n+=2;
      continue;
    } else { 
      cerr << "Unrecognised option " << arg << endl;
      exit(-1);
    } 

  }  // while (n<argc)

  if (globalopts.initmatfname.size()<1) {
    cerr << "Input matrix filename not found" << endl << endl;
    print_usage(argc,argv);
    exit(2);
  }
  if (!globalopts.matonly) {
    cerr << "WARNING:: Using old style options - please update usage" << endl;
  }
  if ((!globalopts.matonly) && (globalopts.testfname.size()<1)) {
    cerr << "ERROR:: Inputvol filename not found" << endl << endl;
  }
  if ((!globalopts.matonly) && (globalopts.reffname.size()<1)) {
    cerr << "ERROR:: Reference volume filename not found" << endl << endl;
  }
}

////////////////////////////////////////////////////////////////////////////

int vector2affine(const ColumnVector& params, Matrix& aff)
{
  // order of parameters is 3 rotation + 3 translation + 3 scales + 3 skews
  // angles are in radians

  ColumnVector centre(3);
  centre = 0;
  compose_aff(params,12,centre,aff,construct_rotmat_euler);
  return 0;
}  


int affmat2vector(const Matrix& aff, ColumnVector& params)
{
  ColumnVector centre(3);
  centre = 0;
  decompose_aff(params,aff,centre,rotmat2euler);
  return 0;
}


////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[])
{
  parse_command_line(argc,argv);


  volume<float> testvol, refvol, intervol;
  if (! globalopts.matonly) {
    // read volumes
    if (read_volume_hdr_only(testvol,globalopts.testfname)<0) {
      cerr << "Cannot read input volume" << endl;
      return -1;
    }
    if (read_volume_hdr_only(refvol,globalopts.reffname)<0) {
      cerr << "Cannot read reference volume" << endl;
      return -1;
    }
    if (globalopts.intervolfname.size()>=1) {
      if (read_volume_hdr_only(intervol,globalopts.intervolfname)<0) {
	cerr << "Cannot read intermediary volume" << endl;
	return -1;
      }
    } else {
      intervol = refvol;
    }
    
    if (globalopts.verbose>3) {
      print_volume_info(refvol,"Reference Volume");
      print_volume_info(testvol,"Input Volume");
    }
  }

  // read matrices
  Matrix affmat(4,4);
  int returnval;
  if (globalopts.matonly)
    returnval = read_ascii_matrix(affmat,globalopts.initmatfname);
  else 
    returnval = read_matrix(affmat,globalopts.initmatfname,testvol,intervol);
  if (returnval<0) {
    cerr << "Cannot read input-matrix" << endl;
    return -2;
  }
    

  if (globalopts.fixfname.size() >= 1) {
    Matrix fixmat(4,4);
    if (globalopts.matonly)
      returnval = read_ascii_matrix(fixmat,globalopts.fixfname);
    else 
      returnval = read_matrix(fixmat,globalopts.fixfname,intervol,refvol);
    
    if (returnval<0) {
      cerr << "Cannot read fixscaleskew-matrix" << endl;
      return -3;
    } else {
      if (globalopts.verbose>2) {
	cout << "Initial matrix:" << endl << affmat << endl;
	cout << "Fix Scale-Skew matrix:" << endl << fixmat << endl;
      }
      // do the work of combining scale/skew from fix and rest from init
      ColumnVector initp(12), fixp(12), combp(12);
      affmat2vector(affmat,initp);
      affmat2vector(fixmat,fixp);
      combp.SubMatrix(1,6,1,1) = initp.SubMatrix(1,6,1,1);
      combp.SubMatrix(7,12,1,1) = fixp.SubMatrix(7,12,1,1);
      vector2affine(combp,affmat);
    }
  }

  
  if (globalopts.concatfname.size() >= 1) {
    Matrix concatmat(4,4);
    if (globalopts.matonly)
      returnval = read_ascii_matrix(concatmat,globalopts.concatfname);
    else 
      returnval = read_matrix(concatmat,globalopts.concatfname,intervol,refvol);
    
    if (returnval<0) {
      cerr << "Cannot read concat-matrix" << endl;
      return -3;
    } else {
      if (globalopts.verbose>2) {
	cout << "Initial matrix:" << endl << affmat << endl;
	cout << "Second matrix:" << endl << concatmat << endl;
      }
      affmat = concatmat * affmat;
    }
  }
  
  // apply inverse (if requested)
  if (globalopts.inverse) {
    affmat = affmat.i();
  }
  
  
  // Write outputs
  if (globalopts.outputmatascii.size() >= 1) {
    write_ascii_matrix(affmat,globalopts.outputmatascii);
  }
  
  if (globalopts.verbose>0) {
    cout << affmat << endl;
  }

  return 0;
}


