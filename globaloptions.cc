/*  globaloptions.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */

#include "globaloptions.h"

globaloptions* globaloptions::gopt = NULL;

void globaloptions::parse_command_line(int argc,char** argv,
				       const string &p_version)
{
  version = p_version;

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
      inputfname = arg;
      n++;
      cerr << "WARNING: change in option usage" << endl << endl;
      cerr << "To specify the input volume the option -in should be used" 
	   << endl << "Accepting the filename for now, but please update "
	   << "to new syntax in future." << endl << endl;
      continue;
    }
    
    // put options without arguments here
    if ( arg == "-help" ) {
      print_usage(argc,argv);
      exit(0);
    } else if ( arg == "-version") {
      print_version();
      exit(0);
    } else if ( arg == "-applyxfm" || arg == "-applynonisoxfm" ) {
      do_optimise = false;
      iso = false;
      nosave = false;
      n++;
      continue;
    } else if ( arg == "-i") {
      interactive = true;
      n++;
      continue;
    } else if ( arg == "-nosearch") {
      searchrx(1) = 0;
      searchrx(2) = 0;
      searchry(1) = 0;
      searchry(2) = 0;
      searchrz(1) = 0;
      searchrz(2) = 0;
      n++;
      continue;
    } else if ( arg == "-nosave") {
      nosave = true;
      n++;
      continue;
    } else if ( arg == "-noresample") {
      resample = false;
      n++;
      continue;
    } else if ( arg == "-forcescaling") {
      force_scaling = true;
      n++;
      continue;
    } else if ( arg == "-2D") {
      mode2D = true;
      n++;
      continue;
    } else if ( arg == "-debugsave") {
      nosave = false;
      n++;
      continue;
   } else if ( arg == "-v" ) {
      verbose = 1;
      n++;
      continue;
    }

    if (n+1>=argc) 
      { 
	cerr << "Lacking argument to option " << arg << endl;
	exit(-1);
      }

    // put options with 1 argument here
    if ( (arg == "-o") || (arg == "-out") ) {
      outputfname = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-ref") {
      reffname = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-in") {
      inputfname = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-init") {
      initmatfname = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-schedule") {
      schedulefname = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-refweight") {
      refweightfname = argv[n+1];
      useweights = true;
      n+=2;
      continue;
    } else if ( arg == "-inweight") {
      testweightfname = argv[n+1];
      useweights = true;
      n+=2;
      continue;
    } else if ( arg == "-omat") {
      outputmatascii = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-omedx") {
      outputmatmedx = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-bins") {
      no_bins = atoi(argv[n+1]);
      n+=2;
      continue;
    } else if ( arg == "-dof") {
      no_params = atoi(argv[n+1]);
      dof = atoi(argv[n+1]);
      n+=2;
      continue;
    } else if ( arg == "-sincwidth") {
      sincwidth = atoi(argv[n+1]);
      n+=2;
      continue;
    } else if ( arg == "-paddingsize") {
      paddingsize = atof(argv[n+1]);
      n+=2;
      continue;
    } else if ( arg == "-applyisoxfm" ) {
      isoscale = atof(argv[n+1]);
      do_optimise = false;
      iso = true;
      nosave = false;
      n+=2;
      continue;
    } else if ( arg == "-basescale") {
      force_basescale = true;
      basescale = atof(argv[n+1]);
      n+=2;
      continue;
    } else if ( arg == "-minsampling") {
      min_sampling = atof(argv[n+1]);
      n+=2;
      continue;
    } else if ( arg == "-coarsesearch") {
      coarsedelta = atof(argv[n+1])*M_PI/180.0;
      n+=2;
      continue;
    } else if ( arg == "-finesearch") {
      finedelta = atof(argv[n+1])*M_PI/180.0;
      n+=2;
      continue;
    } else if ( arg == "-verbose") {
      verbose = atoi(argv[n+1]);
      n+=2;
      continue;
    } else if ( arg == "-datatype") {
      {
	forcedatatype = true;
	string dataarg = argv[n+1];
	if (dataarg == "double") {
	  datatype = DT_DOUBLE;
	} else if (dataarg == "float") {
	  datatype = DT_FLOAT;
	} else if (dataarg == "int") {
	  datatype = DT_SIGNED_INT;
	} else if (dataarg == "short") {
	  datatype = DT_SIGNED_SHORT;
	} else if (dataarg == "char") {
	  datatype = DT_UNSIGNED_CHAR;
	} else {
	  cerr << "Unrecognised data type: " << dataarg << endl;
	  exit(-1);
	}
      }
      n+=2;
      continue;
    } else if ( arg == "-cost") {
      {
	string costarg = argv[n+1];
	if (costarg == "mutualinfo") {
	  maincostfn = MutualInfo;
	} else if (costarg == "corratio") {
	  maincostfn = CorrRatio;
	} else if (costarg == "woods") {
	  maincostfn = Woods;
	} else if (costarg == "normcorr") {
	  maincostfn = NormCorr;
	} else if (costarg == "normmi") {
	  maincostfn = NormMI;
	} else if (costarg == "leastsq") {
	  maincostfn = LeastSq;
	} else {
	  cerr << "Unrecognised cost function type: " << costarg << endl;
	  exit(-1);
	}
      }
      n+=2;
      continue;
    } else if ( arg == "-searchcost") {
      {
	string costarg = argv[n+1];
	if (costarg == "mutualinfo") {
	  searchcostfn = MutualInfo;
	} else if (costarg == "corratio") {
	  searchcostfn = CorrRatio;
	} else if (costarg == "woods") {
	  searchcostfn = Woods;
	} else if (costarg == "normcorr") {
	  searchcostfn = NormCorr;
	} else if (costarg == "normmi") {
	  searchcostfn = NormMI;
	} else if (costarg == "leastsq") {
	  searchcostfn = LeastSq;
	} else {
	  cerr << "Unrecognised cost function type: " << costarg << endl;
	  exit(-1);
	}
      }
      n+=2;
      continue;
    } else if ( arg == "-interp") {
      {
	string interparg = argv[n+1];
	if (interparg == "trilinear") {
	  interpmethod = TriLinear;
	} else if (interparg == "nearestneighbour") {
	  interpmethod = NearestNeighbour;
	} else if (interparg == "sinc") {
	  interpmethod = NEWIMAGE::Sinc;
	} else {
	  cerr << "Unrecognised interpolation method: " << interparg << endl;
	  exit(-1);
	}
      }
      n+=2;
      continue;
    } else if ( arg == "-sincwindow") {
      {
	string winarg = argv[n+1];
	if (winarg == "rectangular") {
	  sincwindow = Rect;
	} else if (winarg == "hanning") {
	  sincwindow = Hanning;
	} else if (winarg == "blackman") {
	  sincwindow = Blackman;
	} else {
	  cerr << "Unrecognised sinc window: " << winarg << endl;
	  exit(-1);
	}
      }
      n+=2;
      continue;
    } else if ( arg == "-anglerep" ) {
      {
	string anglearg = argv[n+1];
	if (anglearg == "quaternion") {
	  anglerep = Quaternion;
	} else if (anglearg == "euler") {
	  anglerep = Euler;
	} else {
	  cerr << "Unrecognised angle representation: " << anglearg << endl;
	  exit(-1);
	}
      }
      n+=2;
      continue;
    }

    if (n+2>=argc) 
      { 
	cerr << "Lacking argument to option " << arg << endl;
	exit(-1);
      }
    
    
    // put options with 2 arguments here
    if ( arg == "-searchrx" ) {
      searchrx(1) = Min(atof(argv[n+1]),atof(argv[n+2]))*M_PI/180.0;
      searchrx(2) = Max(atof(argv[n+1]),atof(argv[n+2]))*M_PI/180.0;
      n+=3;
      continue;
    } else if ( arg == "-searchry" ) {
      searchry(1) = Min(atof(argv[n+1]),atof(argv[n+2]))*M_PI/180.0;
      searchry(2) = Max(atof(argv[n+1]),atof(argv[n+2]))*M_PI/180.0;
      n+=3;
      continue;
    } else if ( arg == "-searchrz" ) {
      searchrz(1) = Min(atof(argv[n+1]),atof(argv[n+2]))*M_PI/180.0;
      searchrz(2) = Max(atof(argv[n+1]),atof(argv[n+2]))*M_PI/180.0;
      n+=3;
      continue;
    } else { 
      cerr << "Unrecognised option " << arg << endl;
      exit(-1);
    } 

    

  }  // while (n<argc)

  if (inputfname.size()<1) {
    cerr << "ERROR:: Input volume filename not found\n\n";
    print_usage(argc,argv);
    exit(2);
  }

  if (reffname.size()<1) {
    cerr << "ERROR:: Reference volume filename not found\n\n";
    print_usage(argc,argv);
    exit(2);
  }
}

void globaloptions::print_usage(int argc, char *argv[])
{
  print_version();
  cout << endl;
  cout << "Usage: " << argv[0] << " [options] -in <inputvol> -ref <refvol> -out <outputvol>\n"
       << "       " << argv[0] << " [options] -in <inputvol> -ref <refvol> -omat <outputmatrix>\n"
       << "       " << argv[0] << " [options] -in <inputvol> -ref <refvol> -applyxfm -init <matrix> -out <outputvol>\n\n"
       << "  Available options are:\n"
       << "        -in  <inputvol>                    (no default)\n"
       << "        -ref <refvol>                      (no default)\n"
       << "        -init <matrix-filname>             (input 4x4 affine matrix - "
                                        << "autodetects format ascii/medx)\n"
       << "        -omat <matrix-filename>            (output in 4x4 ascii format)\n"
       << "        -omedx <matrix-filename>           (output in MEDx format)\n"
       << "        -out, -o <outputvol>               (default is none)\n"
       << "        -datatype {char,short,int,float,double}                    (force output data type)\n"
       << "        -cost {mutualinfo,corratio,normcorr,normmi,leastsq}        (default is corratio)\n"
       << "        -searchcost {mutualinfo,corratio,normcorr,normmi,leastsq}  (default is corratio)\n"
       << "        -anglerep {quaternion,euler}       (default is euler)\n"
       << "        -interp {trilinear,nearestneighbour,sinc}  (final interpolation: def - trilinear)\n"
       << "        -sincwidth <full-width in voxels>  (default is 7)\n"
       << "        -sincwindow {rectangular,hanning,blackman}\n"
       << "        -bins <number of histogram bins>   (default is "
                                        << no_bins << ")\n"
       << "        -dof  <number of transform dofs>   (default is "
                                        << dof << ")\n"
       << "        -noresample                        (do not change input sampling)\n"
       << "        -forcescaling                      (force rescaling even for low-res images)\n"
       << "        -minsampling <vox_dim>             (set minimum voxel dimension for sampling (in mm))\n"
       << "        -basescale <scale>                 (sets the scaling for the final scaling - default is 1.0)\n"
       << "        -applyxfm                          (applies transform (no optimisation) - requires -init)\n"
       << "        -applyisoxfm <scale>               (as applyxfm but forces isotropic resampling)\n"
       << "        -paddingsize <number of voxels>    (for applyxfm: interpolates outside image by size)\n"
       << "        -searchrx <min_angle> <max_angle>  (angles in degrees: default is -90 90)\n" 
       << "        -searchry <min_angle> <max_angle>  (angles in degrees: default is -90 90)\n" 
       << "        -searchrz <min_angle> <max_angle>  (angles in degrees: default is -90 90)\n" 
       << "        -nosearch                          (sets all angular search ranges to 0 0)\n" 
       << "        -coarsesearch <delta_angle>        (angle in degrees: default is 60)\n" 
       << "        -finesearch <delta_angle>          (angle in degrees: default is 18)\n" 
       << "        -schedule <schedule-file>          (replaces default schedule)\n"
       << "        -refweight <volume>                (use weights for reference volume)\n"
       << "        -inweight <volume>                 (use weights for input volume)\n"
       << "        -2D                                (use 2D rigid body mode - ignores dof)\n"
       << "        -verbose <num>                     (0 is least and default)\n"
       << "        -v                                 (same as -verbose 1)\n"
       << "        -i                                 (pauses at each stage: default is off)\n"
       << "        -version                           (prints version number)\n"
       << "        -help\n";
}


void globaloptions::print_version()
{
  cout << "FLIRT version " << version << endl;
}


