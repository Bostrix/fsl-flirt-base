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
  string talfname;
  string imgfname;
  string xfmfname;
  string coordfname;
  bool usetal;
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
  verbose = 0;
  usetal = false;
  mm = false;
  medx=false;
}


////////////////////////////////////////////////////////////////////////////

// Parsing functions for command line parameters

void print_usage(int argc, char *argv[])
{
  cout << "Usage: " << argv[0] << " [options] <filename containing coordinates>\n"
       << "e.g.   " << argv[0] << " -img <invol> -tal <standard image> -xfm <img2standard mat file> <coordinate file>\n"
       << "       " << argv[0] << " -img <invol> <coordinate file>\n"
       << "       " << argv[0] << " -img <invol> - \n\n"
       << "  Options are:\n"
       << "        -img <example IMG volume filename>   (NB: volume, not timeseries)\n"
       << "        -tal <Talairach volume filename>\n"
       << "        -xfm <IMG to Talairach transform filename>\n"
       << "        -vox                                 (input coordinates in voxels - default)\n"
       << "        -mm                                  (input coordinates in mm)\n"
       << "        -flirt                               (use flirt mm coordinate conventions - default)\n"    
       << "        -medx                                (use medx  mm coordinate conventions)\n" 
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
      globalopts.medx = false;
      n++;
      continue;
    } else if ( arg == "-medx" ) {
      globalopts.medx = true;
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
    if ( arg == "-tal") {
      globalopts.talfname = argv[n+1];
      globalopts.usetal = true;
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
    cerr << "ERROR:: IMG volume filename not found\n\n";
  }
  if ((globalopts.usetal) && (globalopts.talfname.size()<1)) {
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
  if (globalopts.usetal && (read_volume_hdr_only(talvol,globalopts.talfname)<0)) {
    cerr << "Cannot read Talairach volume" << endl;
    return -1;
  }
    
  if (globalopts.verbose>3) {
    if (globalopts.usetal) {
      print_info(talvol,"Talairach Volume");
    }
    print_info(imgvol,"IMG Volume");
  }

  // read matrices
  Matrix affmat(4,4);
  int returnval;
  bool use_sform=false;
  if (globalopts.xfmfname.length()>0) {
    returnval = read_matrix(affmat,globalopts.xfmfname,imgvol,talvol);
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

  /////////////// SET UP MATRICES ////////////////

  Matrix vox2tal(4,4);

  if (use_sform) {

    // set the main matrix
    vox2tal = imgvol.sform_mat();
 
    if (imgvol.sform_code()==NIFTI_XFORM_UNKNOWN) { 
      if (globalopts.verbose>0) {
	cerr << "WARNING:: Standard coordinates not set in img" << endl; 
 	// using sampling_mat instead of sform
        vox2tal = imgvol.sampling_mat();
      }
    }
  } else {
    
    // set the main matrix
    vox2tal = talvol.sform_mat() * talvol.sampling_mat().i() * affmat * imgvol.sampling_mat();
    
    if (talvol.sform_code()==NIFTI_XFORM_UNKNOWN) { 
      if (globalopts.verbose>0) {
	cerr << "WARNING:: Standard coordinates not set in tal image" << endl; 
 	// using sampling_mat instead of sform
        vox2tal = affmat * imgvol.sampling_mat();
      }
    }
    
    if (globalopts.verbose>3) {
      cout << " talvox2world =" << endl << talvol.sform_mat() << endl << endl;
    }
    
    // bloody medx conventions
    if (globalopts.medx) {
      Matrix swapy2(4,4);
      Identity(swapy2);
      swapy2(2,2) = -1.0;
      swapy2(2,4) = imgvol.ysize()-1.0;
      vox2tal = vox2tal * swapy2;
    }
    
  }

  // initialise coordinate vectors

  ColumnVector imgcoord(4), talcoord(4), oldimg(4);
  imgcoord = 0;
  talcoord = 0;
  imgcoord(4)=1;
  talcoord(4)=1;
  oldimg = 0;  // 4th component set to 0, so that initially oldimg -ne imgcoord


  bool use_stdin = false;
  if ( (globalopts.coordfname=="-") || (globalopts.coordfname.size()<1)) {
    use_stdin = true;
  }

  if (globalopts.verbose>0) {
    cout << "Coordinates in Talairach volume (in mm):" << endl;
  }

  // set up coordinate reading (from file or stdin) //
  ifstream matfile(globalopts.coordfname.c_str());

  if (use_stdin) {
    if (globalopts.verbose>0) {
      cout << "Please type in IMG coordinates";
      if (globalopts.mm) { 
	cout << " (in mm) :" << endl;
      } else { 
	cout << " (in voxels) :" << endl; 
      }
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
      if (use_stdin) { cin >> imgcoord(j); }
      else { matfile >> imgcoord(j); }
    } 
    if  (use_stdin) {
      // this is in case the pipe continues to input a stream of zeros
      if (oldimg == imgcoord)  return 0;
      oldimg = imgcoord;
    }
    
    if (globalopts.mm) {  // in mm
      talcoord = vox2tal * imgvol.sampling_mat().i() * imgcoord;
    } else { // in voxels
      talcoord = vox2tal * imgcoord; 
    }
    cout << talcoord(1) << "  " << talcoord(2) << "  " << talcoord(3) << endl;
  }
  
  if (!use_stdin) { matfile.close(); }
  
  return 0;
}








