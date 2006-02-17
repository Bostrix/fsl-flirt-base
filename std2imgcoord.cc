/*  std2imgcoord.cc

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
  string stdfname;
  string imgfname;
  string xfmfname;
  string coordfname;
  bool usestd;
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
  stdfname = "";
  imgfname = "";
  coordfname = "";
  xfmfname = "";
  verbose = 0;
  usestd = false;
  mm = true;
}


////////////////////////////////////////////////////////////////////////////

// Parsing functions for command line parameters

void print_usage(int argc, char *argv[])
{
  cout << "Usage: " << argv[0] << " [options] <filename containing coordinates>\n\n"
       << "e.g.   " << argv[0] << " -img <invol> -std <standard image> -xfm <img2standard mat file> <coordinate file>\n"
       << "       " << argv[0] << " -img <invol> <coordinate file>\n"
       << "       " << argv[0] << " -img <invol> - \n\n"
       << "  Options are:\n"
       << "        -std <standard volume filename>\n"
       << "        -img <example image filename>   (NB: 3D image, not timeseries)\n"
       << "        -xfm <image to standard transform filename>\n"
       << "        -mm                                  (outputs coordinates in mm - default)\n"
       << "        -vox                                 (outputs coordinates in voxels)\n"
       << "        -v                                   (verbose output)\n"
       << "        -verbose                             (more verbose output)\n"
       << "        -help\n\n"
       << " Notes:\n"
       << "  (1) if '-' is used as coordinate filename then coordinates are read from standard input\n"
       << "  (2) the -img option is compulsory\n";
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
      cerr << "WARNING::Using outdated options, please update to new usage" << endl;
      // do nothing anymore
      n++;
      continue;
    } else if ( arg == "-v" ) {
      globalopts.verbose = 1;
      n++;
      continue;
    } else if ( arg == "-verbose" ) {
      globalopts.verbose = 5;
      n++;
      continue;
    } else if ( arg == "-") {
      globalopts.coordfname = "-";
      n++;
      continue;
    }

    if (n+1>=argc) 
      { 
	cerr << "Lacking argument to option " << arg << endl;
	break; 
      }

    // put options with 1 argument here
    if ( arg == "-std") {
      globalopts.stdfname = argv[n+1];
      globalopts.usestd = true;
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

  if ((globalopts.imgfname.size()<1)) {
    cerr << "ERROR:: image filename not found\n\n";
  }
  if ((globalopts.usestd) && (globalopts.stdfname.size()<1)) {
    cerr << "ERROR:: standard image filename not found\n\n";
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


  volume<float> imgvol, stdvol;
    // read volumes
  if (read_volume_hdr_only(imgvol,globalopts.imgfname)<0) {
    cerr << "Cannot read input image" << endl;
    return -1;
  }
  if (globalopts.usestd && read_volume_hdr_only(stdvol,globalopts.stdfname)<0) {
    cerr << "Cannot read standard image" << endl;
    return -1;
  }
    
  if (globalopts.verbose>3) {
    if (globalopts.usestd) {
      print_info(stdvol,"standard image");
    }
    print_info(imgvol,"input image");
  }

  // read matrices
  Matrix affmat(4,4);
  int returnval;
  bool use_sform=false;
  if (globalopts.xfmfname.length()>0) {
    returnval = read_matrix(affmat,globalopts.xfmfname,imgvol,stdvol);
    use_sform = false;
    if (returnval<0) {
      cerr << "Cannot read transform file" << endl;
      return -2;
    }
  } else {
    use_sform = true;
    affmat = Identity(4);
  }

    
  if (globalopts.verbose>3) {
    cout << " affmat =" << endl << affmat << endl << endl;
  }

  if (globalopts.verbose>3) {
    cout << " Inverse affmat =" << endl << affmat.i() << endl << endl;
  }


  /////////////// SET UP MATRICES ////////////////

  Matrix vox2std(4,4);

  if (use_sform) {

    // set the main matrix
    vox2std = imgvol.sform_mat();

    if (imgvol.sform_code()==NIFTI_XFORM_UNKNOWN) { 
      if (globalopts.verbose>0) {
	cerr << "WARNING:: standard coordinates not set in img" << endl; 
      }
      // using sform = sampling_mat
      vox2std = imgvol.sampling_mat();
    }
  } else {

    // set the main matrix
    vox2std = stdvol.sform_mat() * stdvol.sampling_mat().i() * affmat * imgvol.sampling_mat();

    if (stdvol.sform_code()==NIFTI_XFORM_UNKNOWN) { 
      if (globalopts.verbose>0) {
	cerr << "WARNING:: standard coordinates not set in standard image" << endl; 
      }
      // using sform = sampling_mat
      vox2std = affmat * imgvol.sampling_mat();
    }
    
    if (globalopts.verbose>3) {
      cout << " stdvox2world =" << endl << stdvol.sform_mat() << endl << endl;
    }


  }


  // initialise coordinate vectors

  ColumnVector imgcoord(4), stdcoord(4), oldstd(4);
  imgcoord = 0;
  stdcoord = 0;
  imgcoord(4)=1;
  stdcoord(4)=1;
  oldstd = 0;  // 4th component set to 0, so that initially oldstd -ne stdcoord


  bool use_stdin = false;
  if ( (globalopts.coordfname=="-") || (globalopts.coordfname.size()<1)) {
    use_stdin = true;
  }

  if (globalopts.verbose>0) {
    if (globalopts.mm) {
      cout << "Coordinates in input image (in mm):" << endl;
    } else {
      cout << "Coordinates in input image (in voxels):" << endl;
    }
  }


  // set up coordinate reading (from file or stdin) //
  ifstream matfile(globalopts.coordfname.c_str());

  if (use_stdin) {
    if (globalopts.verbose>0) {
      cout << "Please type in standard coordinates :" << endl;
    } 
  } else {
    if (!matfile) { 
      cerr << "Could not open matrix file " << globalopts.coordfname << endl;
      return -1;
    }
  }
  
  // loop around reading coordinates and displaying output


 while ( (use_stdin && (!cin.eof())) || ((!use_stdin) && (!matfile.eof())) ) {
    for (int j=1; j<=3; j++) {
      if (use_stdin) { cin >> stdcoord(j); }
      else { matfile >> stdcoord(j); }
    } 
    if  (use_stdin) {
      // this is in case the pipe continues to input a stream of zeros
      if (oldstd == stdcoord)  return 0;
      oldstd = stdcoord;
    }

    if (globalopts.mm) {  // in mm
      imgcoord = imgvol.sampling_mat() * vox2std.i() * stdcoord;
    } else { // in voxels
      imgcoord = vox2std.i() * stdcoord; 
    }

    cout << imgcoord(1) << "  " << imgcoord(2) << "  " << imgcoord(3) << endl;
  }
  
  if (!use_stdin) { matfile.close(); }

  return 0;
}








