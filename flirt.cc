#include <string>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <time.h>
#include <vector>
#include <algorithm>
#define WANT_STREAM
#define WANT_MATH

#include "newmatap.h"
#include "newmatio.h"
#include "mjavwio.h"
#include "AvwRead.h"
#include "miscimfns.h"
#include "miscmaths.h"
#include "interpolation.h"
#include "nrcode.h"
#include "mjimage.h"
#include "costfns.h"
#include "generalio.h"
#include "defaultschedule.h"

#ifndef NO_NAMESPACE
 using namespace MISCMATHS;
 using namespace MISCIMFNS;
 using namespace COSTFNS;
 using namespace INTERPOLATION;
 using namespace NUMREC;
 using namespace MJAVWIO;
 using namespace NEWMAT;
 using namespace MJIMAGE;
 using namespace GENERALIO;
#endif

////////////////////////////////////////////////////////////////////////////

// GLOBAL REQUIREMENTS


  enum costfns { Woods, CorrRatio, MutualInfo, NormCorr, NormMI };
  enum anglereps { Euler, Quaternion };
  typedef std::vector<RowVector> MatVec;
  typedef MatVec* MatVecPtr;

  std::vector<MatVec> globalusrmat(26);
  MatVec globalsearchoptmat, globalpreoptsearchmat;


// the real defaults are provided in the function parse_command_line

class globaloptions {
public:
  string inputfname;
  string outputfname;
  string reffname;
  string outputmatascii;
  string outputmatmedx;
  string initmatfname;
  Matrix initmat;

  string schedulefname;

  imagepair *impair;
  ColumnVector refparams;
  Matrix parammask;
  int no_params;
  int dof;
  int searchdof;
  int no_bins;
  costfns maincostfn;
  costfns searchcostfn;
  costfns currentcostfn;
  anglereps anglerep;
  float min_sampling;

  ColumnVector searchrx;
  ColumnVector searchry;
  ColumnVector searchrz;
  float coarsedelta;
  float finedelta;

  string planereffname;
  string planetestfname;

  int verbose;
  bool interactive;
  bool do_optimise;
  bool measure_cost;
  int defaultplanes;
  bool nosave;
  bool iso;

  int single_param;
public:
  globaloptions();
  ~globaloptions() {};
};

globaloptions globalopts;


globaloptions::globaloptions()
{
  // set up defaults
  outputfname = "";
  reffname = "/usr/local/medx/fmrib/data/structural/mni/average_305-16SI";

  inputfname = "";
  outputmatascii = "";
  outputmatmedx = "";
  initmatfname = "";
  initmat = identity(4);

  schedulefname = "";
  
  impair = 0;
  refparams.ReSize(12);
  refparams << 0.0 << 0.0 << 0.0 << 0.0 << 0.0 << 0.0 << 1.0 << 1.0 << 1.0
	    << 0.0 << 0.0 << 0.0;
  parammask.ReSize(12,12);
  identity(parammask);
  no_params = 12;
  dof = 12;
  searchdof = 12;
  no_bins = 256;
  maincostfn = CorrRatio;
  searchcostfn = CorrRatio;
  currentcostfn = CorrRatio;
  anglerep = Euler;
  min_sampling = 1.0;

  searchrx.ReSize(2);
  searchry.ReSize(2);
  searchrz.ReSize(2);
  searchrx << -M_PI/2.0 << M_PI/2.0;
  searchry << -M_PI/2.0 << M_PI/2.0;
  searchrz << -M_PI << M_PI;
  coarsedelta = 60.0*M_PI/180.0;
  finedelta = 18.0*M_PI/180.0;

  planereffname = "";
  planetestfname = "";

  interactive = false;
  do_optimise = true;
  measure_cost = false;
  defaultplanes = 2;
  nosave = true;
  iso = true;

  single_param = -1;
}



// In addition, each imagepair class contains...
//  All stored inside the class imagepair
//
//  refvol    - The reference volume (const)
//  testvol  - The test image volume (const)
//  outputvol - A volume used to store the transformed test image
//  refmin, refmax - the min/max float value of the intensities in refvol
//  testmin, testmax - the min/max float value of the intensities in testvol
//  plnp - an array of -p*log(p) values used in mutual information
//  hist - an array used to store the joint histogram of the images
//           (updated by mutual information - hence its non-const-ness)


////////////////////////////////////////////////////////////////////////////

// Parsing functions for command line parameters

void print_usage(int argc, char *argv[])
{
  cout << "Usage: " << argv[0] << " [options] <testvol>\n\n"
       << "  Available options are:\n"
       << "        -ref <refvol>                      (default is "
                                        << globalopts.reffname << ")\n"
       << "        -init <matrix-filname>             (4x4 affine matrix - "
                                        << "autodetects format ascii/medx)\n"
       << "        -omat <matrix-filename>            (4x4 ascii format)\n"
       << "        -omedx <matrix-filename>           (MEDx format)\n"
       << "        -out, -o <outputvol>               (default is none)\n"
       << "        -cost {mutualinfo,woods,corratio,normcorr,normmi}        (default is corratio)\n"
       << "        -searchcost {mutualinfo,woods,corratio,normcorr,normmi}  (default is corratio)\n"
       << "        -anglerep {quaternion,euler}       (default is euler)\n"
       << "        -bins <number of histogram bins>   (default is "
                                        << globalopts.no_bins << ")\n"
       << "        -dof  <number of transform dofs>   (default is "
                                        << globalopts.dof << ")\n"
       << "        -applyxfm                          (applies init - "
                                        << "no optimisation)\n"
       << "        -applynonisoxfm                    (as applyxfm but no isotropic resampling)\n"
       << "        -measurecost                       (calculates cost function"
                                        << " - no optimisation)\n"
       << "        -searchrx <min_angle> <max_angle>  (angles in degrees: default is -90 90)\n" 
       << "        -searchry <min_angle> <max_angle>  (angles in degrees: default is -90 90)\n" 
       << "        -searchrz <min_angle> <max_angle>  (angles in degrees: default is -180 180)\n" 
       << "        -coarsesearch <delta_angle>        (angle in degrees: default is 60)\n" 
       << "        -finesearch <delta_angle>          (angle in degrees: default is 18)\n" 
       << "        -schedule <schedule-file>          (replaces default schedule)\n"
       << "        -findplanes                        (default is xy bounding planes)\n"
       << "        -planeref <filename>               (default is none)\n"
       << "        -planetest <filename>              (default is none)\n"
       << "        -plotcost\n"
       << "        -verbose <num>                     (0 is least and default)\n"
       << "        -v                                 (same as -verbose 5)\n"
       << "        -i                                 (pauses at each stage: default is off)\n"
       << "        -nosave                            (do not save intermediate volumes - default)\n"
       << "        -debugsave                         (save any intermediate volumes)\n"
       << "        -help\n";
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
      globalopts.inputfname = arg;
      n++;
      continue;
    }
    
    // put options without arguments here
    if ( arg == "-help" ) {
      print_usage(argc,argv);
      exit(0);
    } else if ( arg == "-applyxfm" ) {
      globalopts.do_optimise = false;
      globalopts.nosave = false;
      n++;
      continue;
    } else if ( arg == "-applynonisoxfm" ) {
      globalopts.do_optimise = false;
      globalopts.iso = false;
      globalopts.nosave = false;
      n++;
      continue;
    } else if ( arg == "-measurecost" ) {
      globalopts.measure_cost = true;
      n++;
      continue;
    } else if ( arg == "-plotcost") {
      globalopts.single_param = 1;
      n++;
      continue;
    } else if ( arg == "-findplanes") {
      globalopts.defaultplanes = 0;
      n++;
      continue;
    } else if ( arg == "-i") {
      globalopts.interactive = true;
      n++;
      continue;
    } else if ( arg == "-nosave") {
      globalopts.nosave = true;
      n++;
      continue;
    } else if ( arg == "-debugsave") {
      globalopts.nosave = false;
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
    if ( (arg == "-o") || (arg == "-out") ) {
      globalopts.outputfname = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-ref") {
      globalopts.reffname = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-init") {
      globalopts.initmatfname = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-schedule") {
      globalopts.schedulefname = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-omat") {
      globalopts.outputmatascii = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-omedx") {
      globalopts.outputmatmedx = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-bins") {
      globalopts.no_bins = atoi(argv[n+1]);
      n+=2;
      continue;
    } else if ( arg == "-dof") {
      globalopts.no_params = atoi(argv[n+1]);
      globalopts.dof = atoi(argv[n+1]);
      n+=2;
      continue;
    } else if ( arg == "-coarsesearch") {
      globalopts.coarsedelta = atof(argv[n+1])*M_PI/180.0;
      n+=2;
      continue;
    } else if ( arg == "-finesearch") {
      globalopts.finedelta = atof(argv[n+1])*M_PI/180.0;
      n+=2;
      continue;
    } else if ( arg == "-verbose") {
      globalopts.verbose = atoi(argv[n+1]);
      n+=2;
      continue;
    } else if ( arg == "-planeref") {
      globalopts.planereffname = argv[n+1];
      globalopts.defaultplanes--;
      n+=2;
      continue;
    } else if ( arg == "-planetest") {
      globalopts.planetestfname = argv[n+1];
      globalopts.defaultplanes--;
      n+=2;
      continue;
    } else if ( arg == "-cost") {
      {
	string costarg = argv[n+1];
	if (costarg == "mutualinfo") {
	  globalopts.maincostfn = MutualInfo;
	} else if (costarg == "corratio") {
	  globalopts.maincostfn = CorrRatio;
	} else if (costarg == "woods") {
	  globalopts.maincostfn = Woods;
	} else if (costarg == "normcorr") {
	  globalopts.maincostfn = NormCorr;
	} else if (costarg == "normmi") {
	  globalopts.maincostfn = NormMI;
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
	  globalopts.searchcostfn = MutualInfo;
	} else if (costarg == "corratio") {
	  globalopts.searchcostfn = CorrRatio;
	} else if (costarg == "woods") {
	  globalopts.searchcostfn = Woods;
	} else if (costarg == "normcorr") {
	  globalopts.searchcostfn = NormCorr;
	} else if (costarg == "normmi") {
	  globalopts.searchcostfn = NormMI;
	} else {
	  cerr << "Unrecognised cost function type: " << costarg << endl;
	  exit(-1);
	}
      }
      n+=2;
      continue;
    } else if ( arg == "-anglerep" ) {
      {
	string anglearg = argv[n+1];
	if (anglearg == "quaternion") {
	  globalopts.anglerep = Quaternion;
	} else if (anglearg == "euler") {
	  globalopts.anglerep = Euler;
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
      globalopts.searchrx(1) = Min(atof(argv[n+1]),atof(argv[n+2]))*M_PI/180.0;
      globalopts.searchrx(2) = Max(atof(argv[n+1]),atof(argv[n+2]))*M_PI/180.0;
      
      n+=3;
      continue;
    } else if ( arg == "-searchry" ) {
      globalopts.searchry(1) = Min(atof(argv[n+1]),atof(argv[n+2]))*M_PI/180.0;
      globalopts.searchry(2) = Max(atof(argv[n+1]),atof(argv[n+2]))*M_PI/180.0;
      n+=3;
      continue;
    } else if ( arg == "-searchrz" ) {
      globalopts.searchrz(1) = Min(atof(argv[n+1]),atof(argv[n+2]))*M_PI/180.0;
      globalopts.searchrz(2) = Max(atof(argv[n+1]),atof(argv[n+2]))*M_PI/180.0;
      n+=3;
      continue;
    } else { 
      cerr << "Unrecognised option " << arg << endl;
      exit(-1);
    } 

    

  }  // while (n<argc)

  if (globalopts.inputfname.size()<1) {
    cerr << "Input filename not found\n\n";
    print_usage(argc,argv);
    exit(2);
  }
}


//------------------------------------------------------------------------//
// Some interfaces to generalio

int safe_save_volume(const volume& source, const string& filename)
{
  if (!globalopts.nosave) {
    save_volume(source,filename);
  }
  return 0;
}


void save_matrix_data(const Matrix& matresult, const volume& initvol, 
		      const volume& finalvol)
{
    write_ascii_matrix(matresult,globalopts.outputmatascii);
    write_medx_matrix(matresult,globalopts.outputmatmedx,
		      initvol,finalvol,"a",globalopts.reffname);
}

void save_global_data(const Matrix& matresult, const volume& initvol,
		      const volume& finalvol)
{
    // save the data now
    safe_save_volume(globalopts.impair->testvol,"inputvol");
    safe_save_volume(globalopts.impair->refvol,"refvol");
    volume outputvol;
    outputvol = globalopts.impair->refvol;
    filled_affine_transform(outputvol,globalopts.impair->testvol,
		     matresult * globalopts.initmat);
    safe_save_volume(outputvol,globalopts.outputfname.c_str());
    save_matrix_data(matresult * globalopts.initmat,initvol,finalvol);
}

void save_global_data(const Matrix& matresult) {
  save_global_data(matresult,globalopts.impair->testvol,
		   globalopts.impair->testvol);
}


float costfn(const Matrix& matresult);

void save_workspace_and_pause(const Matrix& matresult) {
  Tracer tr("save_workspace_and_pause");
  costfn(matresult);      
  save_global_data(matresult);
  if (globalopts.interactive) {
    cerr << "Enter a number to continue..."; // forces the transform to be done
    int dummy; 
    cin >> dummy;
  }
}

//////////////////////////////////////////////////////////////////////////
// Plane-fit stuff

void set_default_planes(volume& vol)
{
  // the planept points should be in VOXEL coordinates
  //  but all saved values should be in MEASUREMENT coordinates
  ColumnVector axis1(3), planept1(3), axis2(3), planept2(3);
  axis1 << 0.0 << 0.0 << 1.0;
  axis2 << 0.0 << 0.0 << -1.0;
  planept1 << 0.0 << 0.0 << -0.5;
  planept2 << 0.0 << 0.0 << (float) vol.slices() - 0.5;
  
  // convert to measurement space
  axis1 = vol.sampling_matrix().SubMatrix(1,3,1,3) * axis1;
  axis2 = vol.sampling_matrix().SubMatrix(1,3,1,3) * axis2;
  planept1 = vol.sampling_matrix().SubMatrix(1,3,1,3) * planept1
    + vol.sampling_matrix().SubMatrix(1,3,4,4);
  planept2 = vol.sampling_matrix().SubMatrix(1,3,1,3) * planept2
    + vol.sampling_matrix().SubMatrix(1,3,4,4);

  // reinforce normalisation of axes
  axis1 = axis1/norm2(axis1);
  axis2 = axis2/norm2(axis2);
  vol.setplanes(axis1,planept1,axis2,planept2);
}


void fit_planes(volume& vol, const string& fname)
{
  // the planept points should be in VOXEL coordinates
  //  but all saved values should be in MEASUREMENT coordinates
  ColumnVector axis1(3), planept1(3), axis2(3), planept2(3);
  int rval=0;
  if (fname.size()>1) {
    Matrix planeparams(4,3);
    rval = read_ascii_matrix(planeparams,fname);
    if (rval>=0) {
      axis1 = planeparams.SubMatrix(1,1,1,3).t();
      planept1 = planeparams.SubMatrix(2,2,1,3).t();
      axis2 = planeparams.SubMatrix(3,3,1,3).t();
      planept2 = planeparams.SubMatrix(4,4,1,3).t();
    }
  } else {
    if ((globalopts.defaultplanes>0) || (rval<0)) {
      // do nothing as the previously set default planes account for padding
      vol.getplanes(axis1,planept1,axis2,planept2);
    } else {
      // make a gradient volume
      if (globalopts.verbose >= 1) 
	cout << "Calculating Gradient" << endl;
      float p90=0.0;
      volume testgrad;
      
      gradient(vol,testgrad,true);      
      if (globalopts.verbose >= 1) 
	cout << "Calculating plane fit for vol" << endl;
      p90=calculate_planefit(axis1,planept1,axis2,planept2,testgrad,
			     globalopts.verbose);
      // convert to measurement space
      axis1 = vol.sampling_matrix().SubMatrix(1,3,1,3) * axis1;
      axis2 = vol.sampling_matrix().SubMatrix(1,3,1,3) * axis2;
      planept1 = vol.sampling_matrix().SubMatrix(1,3,1,3) * planept1
	+ vol.sampling_matrix().SubMatrix(1,3,4,4);
      planept2 = vol.sampling_matrix().SubMatrix(1,3,1,3) * planept2
	+ vol.sampling_matrix().SubMatrix(1,3,4,4);
    }
  }

  // reinforce normalisation of axes
  axis1 = axis1/norm2(axis1);
  axis2 = axis2/norm2(axis2);
  vol.setplanes(axis1,planept1,axis2,planept2);
  if (globalopts.verbose>=3) {
    cout << "Plane parameters:\n axis1    = " << axis1.t()
	 << " planept1 = " << planept1.t() 
	 << " axis2    = " << axis2.t() 
	 << " planept2 = " << planept2.t();
  }
}


//////////////////////////////////////////////////////////////////////////

// OPTIMISATION SUPPORT


int vector2affine(const ColumnVector& params, int n, const ColumnVector& centre,
		  Matrix& aff)
{
  if (n<=0) return 0;
  // order of parameters is 3 rotation + 3 translation + 3 scales + 3 skews
  // angles are in radians

  switch (globalopts.anglerep) 
    {
    case Euler:
      compose_aff(params,n,centre,aff,construct_rotmat_euler);
      break;
    case Quaternion:
      compose_aff(params,n,centre,aff,construct_rotmat_quat);
      break;
    default:
      cerr << "Invalid Rotation Representation" << endl;
      return -1;
    }
  return 0;
}  


int vector2affine(const ColumnVector& params, int n, Matrix& aff)
{
  return vector2affine(params,n,globalopts.impair->testvol.cog(),aff);
}


int vector2affine(const float params[], int n, Matrix& aff)
{
  ColumnVector vparams(12);
  for (int i=1; i<=n; i++) {
    vparams(i) = params[i];
  }
  return vector2affine(vparams,n,aff);
}


int affmat2vector(const Matrix& aff, int n, const ColumnVector& centre,
		  ColumnVector& params)
{
  switch (globalopts.anglerep) 
    {
    case Euler:
      decompose_aff(params,aff,centre,rotmat2euler);
      break;
    case Quaternion:
      decompose_aff(params,aff,centre,rotmat2quat);
      break;
    default:
      cerr << "Invalid Rotation Representation" << endl;
    }
  return 0;
}


int affmat2vector(const Matrix& aff, int n, ColumnVector& params)
{
  return affmat2vector(aff,n,globalopts.impair->testvol.cog(),params);
}


void set_param_basis(Matrix &parambasis, int no_params)
{
  parambasis = 0.0;
  for (int i=1; i<=no_params; i++) {
    parambasis(i,i)=1.0;
  }
}

void set_param_tols(ColumnVector &param_tol, int no_params)
{
      // Tolerances are: 0.57 degrees (0.005 radians), 0.2mm translation
      //    0.02 scale and 0.001 skew
  float diagonal[12]={0.005, 0.005, 0.005, 0.2, 0.2, 0.2, 0.002, 0.002, 0.002,
  		      0.001, 0.001, 0.001};
  if (param_tol.Nrows()<no_params) {
    param_tol.ReSize(no_params);
  }
  for (int i=1; i<=no_params; i++) {
    param_tol(i)=diagonal[i-1];
  }
}


void set_simplex_start(Matrix &pts, int no_params, const ColumnVector& params)
{
  Matrix parambasis(no_params,no_params);
  set_param_basis(parambasis,no_params);
  for (int i=1; i<=no_params; i++) {
    pts(1,i) = params(i);
  }
  for (int j=1; j<=no_params; j++) {
    for (int i=1; i<=no_params; i++) {
      pts(j+1,i) = params(i) + parambasis(j,i);
    }
  }
}


void initialise_params(ColumnVector& params)
{
  Real paramsf[12] = {0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0};
  params << paramsf;
}


void simplex(ColumnVector& params, int no_params, ColumnVector& param_tol,
	      int *no_its, float *fans, float (*costfunc)(float []))
{
  // sets up the initial parameters and calls the amoeba simplex routine
  initialise_params(params);
  Matrix simplex_pts(no_params+1,no_params);
  set_simplex_start(simplex_pts,no_params,params);
  float y[13];
  for (int i=1; i<=no_params; i++) { y[i] = params(i); }
  float ptol[13];
  for (int i=1; i<=no_params; i++) { ptol[i] = param_tol(i); }
  //amoeba(simplex_pts,y,no_params,ptol,costfunc,no_its);
  for (int i=1; i<=no_params; i++) { params(i) = simplex_pts(1,i); }
  *fans = y[1];
}


void powell_opt(ColumnVector& params, int no_params, ColumnVector& param_tol, 
		int *no_its, float *fans, float (*costfunc)(float []), 
		int itmax)
{
  // sets up the initial parameters and calls the powell optimisation routine
  if (params.MaximumAbsoluteValue() < 0.001)  initialise_params(params);
  { 
    Matrix affmattst(4,4);
    vector2affine(params,no_params,affmattst);
    if (globalopts.verbose>=5) {
      cout << "Starting with : " << endl << affmattst;
    }
  }
  Matrix parambasis(no_params,no_params);
  set_param_basis(parambasis,no_params);
  float ptol[13];
  for (int i=1; i<=no_params; i++) { ptol[i] = param_tol(i); }
  
  // the optimisation call
  powell(params,parambasis,no_params,ptol,no_its, fans, costfunc, itmax);
}


void optimise(ColumnVector& params, int no_params, ColumnVector& param_tol, 
	      int *no_its, float *fans, float (*costfunc)(float []), 
	      int itmax=4)
{
  powell_opt(params,no_params,param_tol,no_its,fans,costfunc,itmax);
}

////////////////////////////////////////////////////////////////////////////

// OPTIMISATION SUPPORT (cost function interfaces)


float costfn(const Matrix& uninitaffmat)
{
  Tracer tr("costfn");
  Matrix affmat = uninitaffmat * globalopts.initmat;  // apply initial matrix
  float retval = 0.0;
  switch (globalopts.currentcostfn) 
    {
    case NormCorr:
      retval = 1.0 - fabs(normcorr(globalopts.impair,affmat));  // MAXimise corr
      break;
    case CorrRatio:
      retval = 1.0 - corr_ratio(globalopts.impair,affmat);  // MAXimise corr
      break;
    case Woods:
      retval = woods_fn(globalopts.impair,affmat);  // minimise variance/mean
      break;
    case MutualInfo:
      retval = -mutual_info(globalopts.impair,affmat);  // MAXimise info
      break;
    case NormMI:
      retval = -normalised_mutual_info(globalopts.impair,affmat);  // MAXimise
      break;
    default:
      cerr << "Invalid cost function type" << endl;
      break;
    }
  return retval;
}


float costfn(const ColumnVector& params)
{
  Tracer tr("costfn");
  Matrix affmat(4,4);
  vector2affine(params,globalopts.no_params,affmat);
  float retval = costfn(affmat);
  return retval;
}
  

float costfn(float params[])
{
  Tracer tr("costfn");
  Matrix affmat(4,4);
  vector2affine(params,globalopts.no_params,affmat);
  float retval = costfn(affmat);
  if (globalopts.verbose>=5) {
    cout << globalopts.impair->count++ << " : ";
    cout << retval << " :: ";
    for (int i=1; i<=globalopts.no_params; i++) cout << params[i] << " ";
    cout << endl;
  }
  return retval;
}

//----------------------------------------------------------------------//


void params12toN(ColumnVector& params)
{
  // Convert the full 12 dof param vector to a small param vector
  Tracer tr("params12toN");
  ColumnVector nparams;
  nparams = pinv(globalopts.parammask)*(params - globalopts.refparams);
  params = nparams;
}


void paramsNto12(ColumnVector& params)
{
  // Convert small param vector to full 12 dof param vector
  Tracer tr("paramsNto12");
  ColumnVector param12;
  param12 = globalopts.parammask*params + globalopts.refparams;
  params = param12;
}



float subset_costfn(const ColumnVector& params)
{
  Tracer tr("subset_costfn");
  ColumnVector param12;
  param12 = params;
  paramsNto12(param12);
  float retval = costfn(param12);
  if (globalopts.verbose>=7) {
    cout << globalopts.impair->count++ << " : ";
    cout << retval << " :: " << param12.t() << endl;
  }
  return retval;
}


float subset_costfn(float params[])
{
  Tracer tr("subset_costfn");
  ColumnVector newparams(globalopts.parammask.Ncols());
  for (int n=1; n<=globalopts.parammask.Ncols(); n++) {
    newparams(n) = params[n];
  }
  return subset_costfn(newparams);
}


////////////////////////////////////////////////////////////////////////////

// DIAGNOSTIC - plots the cost function (using gnuplot)

void costgnuplot(const ColumnVector& param0)
{
  Tracer tr("costgnuplot");
  volume validmask;
  int no_points = 20;
  ColumnVector dirn(12);
  float dim=0.0, min=-1.0, max=1.0;
  while (1) {
    dim = 0.0;
    while ((dim<1.0) || (dim>globalopts.no_params)) {
      cout << "Enter the dimension to traverse (1 to "
	   << globalopts.no_params << ", -1 to exit): ";
      cin >> dim;
      if (dim==-1.0)  return;
    }
    dirn = 0.0;
    dirn((int) dim) = 1.0;

    cout << "Enter minimum and maximum limits for deviations from zero: ";
    cin >> min;
    cin >> max;

    output_costplot(param0,dirn,globalopts.no_params,min,max,no_points,costfn,
		    "costplot.txt");
    { 
      ofstream fstr("plotcost");
      fstr << "plot 'costplot.txt'" << endl << "pause -1" << endl;
      fstr.close();
    }
    cout << endl << "Hit return (in this window) to continue" << endl;
    system("gnuplot plotcost");
  }
}


void costplot(void) {
  Tracer tr("costplot");
  int sp = globalopts.single_param;
  if (sp<=0) return;
  float scalefactor = 1.0, fans=0.0;
  ColumnVector params_8(12), param_tol(12);
  Matrix paramres(2*globalopts.no_params+2,100);

  // loop while the user wants to continue
  while (sp>0) {
    set_param_tols(param_tol,12);
    float numoffset = (float) (paramres.Ncols()/2);

    cout << "Enter a central " <<globalopts.no_params<< "-paramater set\n";
    for (int n=1; n<=globalopts.no_params; n++) {
      cin >> globalopts.refparams(n);
    }
    cout << "Enter a scaling factor for the tolerance ";
    cin >> scalefactor;

    // calculate the cost function along each parameter axis
    for (int single_param=1; single_param <= globalopts.no_params;
	 single_param++) {
      for (int num=1; num<=paramres.Ncols(); num++) {
	params_8 = globalopts.refparams;
	params_8(single_param) += 
	  ((float) num-numoffset)*scalefactor*param_tol(single_param);
	fans = costfn(params_8);
	paramres(2*single_param-1,num) = params_8(single_param);
	paramres(2*single_param,num) = fans;
      }
      cerr << "*";
    }

    // now do the joint scaling (all scales equal) case
    for (int num=1; num<=paramres.Ncols(); num++) {
      params_8 = globalopts.refparams;
      params_8(7) += ((float) num-numoffset)*scalefactor*param_tol(7);
      params_8(8) = params_8(7);
      params_8(9) = params_8(7);
      fans = costfn(params_8);
      paramres(2*globalopts.no_params+1,num) = params_8(7);
      paramres(2*globalopts.no_params+2,num) = fans;
    }
    cerr << endl;

    // miscellaneous trace writes
    /*
      params_8 = globalopts.refparams;
      fans = costfn(params_8);
      globalopts.impair->outputvol.create_valid_mask(validmask);
      safe_save_volume(validmask,"validmask");
      safe_save_volume(globalopts.impair->outputvol,"t2tst");
    */
    cerr << "Writing results to costplot.txt" << endl;
    write_ascii_matrix(paramres,"costplot.txt");


    cout << "\nAnother central paramater? (enter -1 to exit) ";
    cin >> sp;
  }
}


//------------------------------------------------------------------------//


float estimate_scaling(const volume& vol) {
  Tracer tr("estimate_scaling");
  return Min(Min(vol.getx(),vol.gety()),vol.getz());
}

float estimate_scaling() {
  Tracer tr("estimate_scaling");
  return estimate_scaling(globalopts.impair->refvol);
}




void find_cost_minima(Matrix& bestpts, const volume& cost) {
  Tracer tr("find_cost_minima");
  volume minv(cost);
  ColumnVector bestpt(3);
  minv = 0.0;
  bestpt = 0.0;
  int xb=cost.columns(), yb=cost.rows(), zb=cost.slices(), minima=0;
  // calculate the number of neighbouring voxels *less* than the centre one
  for (int z=0; z<zb; z++) {
    for (int y=0; y<yb; y++) {
      for (int x=0; x<xb; x++) {
	for (int zo=-1; zo<=1; zo++) {
	  for (int yo=-1; yo<=1; yo++) {
	    for (int xo=-1; xo<=1; xo++) {
	      if (cost.in_bounds(x+xo,y+yo,z+zo)) {
		if (cost(x+xo,y+yo,z+zo) < cost(x,y,z)) {
		  minv(x,y,z) += 1.0;
		}
	      }
	    }
	  }
	}
	if (cost.in_bounds(x,y,z)) {
	  if (cost(x,y,z) < cost((int)bestpt(1),(int)bestpt(2),(int)bestpt(3)))
	    {
	      bestpt(1) = (float) x;
	      bestpt(2) = (float) y;
	      bestpt(3) = (float) z;
	    }
	  if ( cost.in_bounds(x+1,y+1,z+1) && (minv(x,y,z) < 0.5) )
	    minima++;
	}
      }
    }
  }
  //px=0; py=0; pz=0; 
  int idx=1;
  if (minima<=0) {
    bestpts.ReSize(1,16);
    bestpts(1,1) = bestpt(1);
    bestpts(1,2) = bestpt(2);
    bestpts(1,3) = bestpt(3);
    return;
  }
  bestpts.ReSize(minima,16);
  for (int z=0; z<zb; z++) {
    for (int y=0; y<yb; y++) {
      for (int x=0; x<xb; x++) {
	//if (fabs(minv(x,y,z) - minv(px,py,pz)) < 0.5) {
	//if (cost(x,y,z) < cost(px,py,pz)) {
	    // use new point if it is the same minv but lower cost
	        //px=x; py=y; pz=z;
	    //}
	//}
	//if (minv(x,y,z) < (minv(px,py,pz)-0.5)) {
	  // always prefer a lower minv (regardless of the cost)
	  // NB: -0.5 used to ensure that it is a lower integer
	      //px=x; py=y; pz=z;
	  //}
	// trace write
	if ( cost.in_bounds(x,y,z) && cost.in_bounds(x+1,y+1,z+1) ) {
	  if (minv(x,y,z) < 0.5) {
	    bestpts(idx,1) = (float) x; 
	    bestpts(idx,2) = (float) y; 
	    bestpts(idx,3) = (float) z;
	    idx++;
	    if (globalopts.verbose>=3)
	      cerr << "COST minima at : " << x << "," << y << "," << z << endl;
	  }
	}
      }
    }  
  }
}
  

int round(const float val) {
  Tracer tr("round");
  return (int) (val + 0.5);
}


void set_rot_sampling(ColumnVector& rots, float lowerbound, float upperbound) {
  Tracer tr("set_rot_sampling");
  // this function RELIES on the number of rows in the vector rots being set
  int nmax = rots.Nrows();
  if (nmax==1) { rots(1) = (upperbound + lowerbound) / 2.0;  return; }
  for (int n=1; n<=nmax; n++) {
    rots(n) = lowerbound + ((float) (n-1)) * (upperbound - lowerbound) / 
                                                        ((float) (nmax-1));
  }
}
    

void set_rot_samplings(ColumnVector& rxcoarse, ColumnVector& rycoarse,
		       ColumnVector& rzcoarse, ColumnVector& rxfine,
		       ColumnVector& ryfine, ColumnVector& rzfine) {
  Tracer tr("set_rot_samplings");
  //int coarsesize = 4, finesize = 11;
  // sets the number of rows (angle samples) for each axis
  rxcoarse.ReSize(round((globalopts.searchrx(2) - globalopts.searchrx(1))
			/globalopts.coarsedelta)+1);
  rycoarse.ReSize(round((globalopts.searchry(2) - globalopts.searchry(1))
			/globalopts.coarsedelta)+1);
  rzcoarse.ReSize(round((globalopts.searchrz(2) - globalopts.searchrz(1))
			/globalopts.coarsedelta)+1);
  rxfine.ReSize(round((globalopts.searchrx(2) - globalopts.searchrx(1))
			/globalopts.finedelta)+1);
  ryfine.ReSize(round((globalopts.searchry(2) - globalopts.searchry(1))
			/globalopts.finedelta)+1);
  rzfine.ReSize(round((globalopts.searchrz(2) - globalopts.searchrz(1))
			/globalopts.finedelta)+1);
  // now get the appropriate angle sample values
  set_rot_sampling(rxcoarse,globalopts.searchrx(1),globalopts.searchrx(2));
  set_rot_sampling(rycoarse,globalopts.searchry(1),globalopts.searchry(2));
  set_rot_sampling(rzcoarse,globalopts.searchrz(1),globalopts.searchrz(2));
  set_rot_sampling(rxfine,globalopts.searchrx(1),globalopts.searchrx(2));
  set_rot_sampling(ryfine,globalopts.searchry(1),globalopts.searchry(2));
  set_rot_sampling(rzfine,globalopts.searchrz(1),globalopts.searchrz(2));
  if (globalopts.verbose>=4) {
    cout << "Coarse rotation samplings are:\n" << rxcoarse.t() << rycoarse.t() 
	 << rzcoarse.t() << " and fine rotation samplings are:\n"
	 << rxfine.t() << ryfine.t() << rzfine.t() << endl;
  }
}



void search_cost(Matrix& paramlist, volume& costs, volume& tx, 
		 volume& ty, volume& tz, volume& scale) {
  Tracer tr("search_cost");
  int storedverbose = globalopts.verbose;
  int storeddof = globalopts.dof;
  anglereps useranglerep = globalopts.anglerep;
  globalopts.verbose -= 2;
  globalopts.anglerep = Euler;  // a workaround hack
  globalopts.currentcostfn = globalopts.searchcostfn;
  globalopts.dof = globalopts.searchdof;

  ColumnVector coarserx, coarsery, coarserz, finerx, finery, finerz;
  set_rot_samplings(coarserx,coarsery,coarserz,finerx,finery,finerz);

  // set up the type of parameter subset (via the global mask)
  // here 3 translations and 1 (common) scaling are used
  if (globalopts.dof>6) {
    globalopts.parammask.ReSize(12,4);  // was 3
    globalopts.parammask = 0.0;
    globalopts.parammask(7,1) = 1.0;  // didn't used to exist
    globalopts.parammask(8,1) = 1.0;  // didn't used to exist
    globalopts.parammask(9,1) = 1.0;  // didn't used to exist
    globalopts.parammask(4,2) = 1.0;
    globalopts.parammask(5,3) = 1.0;
    globalopts.parammask(6,4) = 1.0;
  } else {
    globalopts.parammask.ReSize(12,3);
    globalopts.parammask = 0.0;
    globalopts.parammask(4,1) = 1.0;
    globalopts.parammask(5,2) = 1.0;
    globalopts.parammask(6,3) = 1.0;
  }
	
  ColumnVector param_tol, param_tol0(12), param_tol1(12), params_8(12);
  globalopts.no_params = 12; // necessary for any subset_costfn call
  param_tol0 = globalopts.refparams;
  params12toN(param_tol0);
  set_param_tols(param_tol1,12);
  param_tol1 = estimate_scaling()*param_tol1;  
  param_tol1 = param_tol1 + globalopts.refparams;
  params12toN(param_tol1);
  param_tol = param_tol1 - param_tol0;

  // search coarsely, optimising each point and storing the results
  int no_its=0;
  float fans=0.0, rx,ry,rz;
  Matrix affmat(4,4);
  ColumnVector trans(3), testv(4), testv2(4);
  tx.reinitialize(coarserx.Nrows(),coarsery.Nrows(),coarserz.Nrows());
  ty.reinitialize(coarserx.Nrows(),coarsery.Nrows(),coarserz.Nrows());
  tz.reinitialize(coarserx.Nrows(),coarsery.Nrows(),coarserz.Nrows());
  scale.reinitialize(coarserx.Nrows(),coarsery.Nrows(),coarserz.Nrows());
  // fix the reference parameter (starting estimates)
  globalopts.refparams = 0.0;
  globalopts.refparams(7) = 1.0;
  globalopts.refparams(8) = 1.0;
  globalopts.refparams(9) = 1.0;
  // set the initial translation (to align cog's)
  trans = globalopts.impair->refvol.cog() - globalopts.impair->testvol.cog();
  globalopts.refparams(4) = trans(1);
  globalopts.refparams(5) = trans(2);
  globalopts.refparams(6) = trans(3);
  for (int ix=0; ix<coarserx.Nrows(); ix++) {
    for (int iy=0; iy<coarsery.Nrows(); iy++) {
      for (int iz=0; iz<coarserz.Nrows(); iz++) {
	rx = coarserx(ix+1);
	ry = coarsery(iy+1);
	rz = coarserz(iz+1);
	globalopts.refparams(1) = rx;
	globalopts.refparams(2) = ry;
	globalopts.refparams(3) = rz;
	params_8 = globalopts.refparams;
	if (globalopts.verbose>=4) {
	  cout << "Starting with " << params_8.t();
	  cout << "  and tolerance " << param_tol.t();
	}
	params12toN(params_8);
	optimise(params_8,globalopts.parammask.Ncols(),
		 param_tol,&no_its,&fans,subset_costfn);
	paramsNto12(params_8);
	tx(ix,iy,iz) = params_8(4);
	ty(ix,iy,iz) = params_8(5);
	tz(ix,iy,iz) = params_8(6);
	scale(ix,iy,iz) = params_8(7);
	
	if (globalopts.verbose>=4) {
	  cout << " dearranged: " << params_8.t();
	}
      }
      if (globalopts.verbose>=2) cerr << "*";
    }
  }
  if (globalopts.verbose>=2) cerr << endl;

  // scale = 1.0;  // for now disallow non-unity scalings

  if (globalopts.verbose>=4) {
    safe_save_volume(tx,"tx");
    safe_save_volume(ty,"ty");
    safe_save_volume(tz,"tz");
    safe_save_volume(scale,"sc");
  }

  // use the optimised translation (and scale) to interpolate parameters
  //  and find the costs at many more sites (without excessive time)
  float xf, yf, zf, txv, tyv, tzv, scv;
  float factorx = ((float) coarserx.Nrows()-1) 
                      / Max((float) 1.0,((float) finerx.Nrows()-1));
  float factory = ((float) coarsery.Nrows()-1)
                      / Max((float) 1.0,((float) finery.Nrows()-1));
  float factorz = ((float) coarserz.Nrows()-1)
                      / Max((float) 1.0,((float) finerz.Nrows()-1));
  costs.reinitialize(finerx.Nrows(),finery.Nrows(),finerz.Nrows());
  for (int ix=0; ix<finerx.Nrows(); ix++) {
    for (int iy=0; iy<finery.Nrows(); iy++) {
      for (int iz=0; iz<finerz.Nrows(); iz++) {
	rx = finerx(ix+1);
	ry = finery(iy+1);
	rz = finerz(iz+1);
	globalopts.refparams(1) = rx;
	globalopts.refparams(2) = ry;
	globalopts.refparams(3) = rz;
	xf = ((float) ix)*factorx;
	yf = ((float) iy)*factory;
	zf = ((float) iz)*factorz;
	txv = tri_interpolation(tx,xf,yf,zf);
	tyv = tri_interpolation(ty,xf,yf,zf);
	tzv = tri_interpolation(tz,xf,yf,zf);
	scv = tri_interpolation(scale,xf,yf,zf);
	if ((scv<0.5) || (scv>2.0))  scv = 1.0;
	if (globalopts.dof<=6) scv = 1.0;
	globalopts.refparams(4) = txv;
	globalopts.refparams(5) = tyv;
	globalopts.refparams(6) = tzv;
	globalopts.refparams(7) = scv;
	globalopts.refparams(8) = scv;
	globalopts.refparams(9) = scv;
	params_8 = globalopts.refparams;
	costs(ix,iy,iz) = costfn(params_8);
      }
      if (globalopts.verbose>=2) cerr << "*";
    }
  }
  if (globalopts.verbose>=2) cerr << endl;

  if (globalopts.verbose>=4) {
    safe_save_volume(costs,"costs");
  }

  // now search to find the points that are local minima
  //globalopts.no_params = 7;
  //int ix,iy,iz;
  //Matrix bestpts;
  //find_cost_minima(bestpts,costs);


  // for each local minima found optimise the translation and scale
  //set_param_tols(param_tol,12);
  //param_tol = 8*param_tol;  
  //globalopts.no_params = 12;  // temp test
  //params12toN(param_tol);
  /*
  for (int n=1; n<=bestpts.Nrows(); n++) {
    ix = (int) bestpts(n,1);
    iy = (int) bestpts(n,2);
    iz = (int) bestpts(n,3);
    cerr << "Cost minima at : " << ix << "," << iy << "," << iz << endl;
    rx = finerx(ix+1);
    ry = finery(iy+1);
    rz = finerz(iz+1);
    xf = ((float) ix)*factorx;
    yf = ((float) iy)*factory;
    zf = ((float) iz)*factorz;
    txv = tri_interpolation(tx,xf,yf,zf);
    tyv = tri_interpolation(ty,xf,yf,zf);
    tzv = tri_interpolation(tz,xf,yf,zf);
    scv = tri_interpolation(scale,xf,yf,zf);
    if ((scv<0.5) || (scv>2.0))  scv = 1.0;
    if (globalopts.dof<=6) scv = 1.0;
    params_8 = 0.0;
    params_8(1) = rx;  params_8(2) = ry;  params_8(3) = rz;
    params_8(4) = txv; params_8(5) = tyv; params_8(6) = tzv;
    params_8(7) = scv; params_8(8) = scv; params_8(9) = scv;
    globalopts.refparams = params_8;
    params12toN(params_8);
    optimise(params_8,globalopts.parammask.Ncols(),
	     param_tol,&no_its,&fans,subset_costfn);
    paramsNto12(params_8);
    bestpts(n,4) = fans;
    bestpts.SubMatrix(n,n,5,16) = params_8.t();
  }
  if (globalopts.verbose>=3) {
    cout << "Costs (4th column) are:\n" << bestpts << endl;
  }

  // find the single best cost
  int bestidx = 1;
  for (int n=1; n<=bestpts.Nrows(); n++) {
    if (bestpts(n,4) < bestpts(bestidx,4)) {
      bestidx = n;
    }
  }
  */

  //int bestidx = 1;
  float costmin, costmax;
  get_min_max(costs,costmin,costmax);
  // the following assumes that costmin > 0
  float factor = 0.2;
  float costthresh = Min(costmin + factor*(costmax-costmin),
			 get_percentile(costs,20));
  // avoid the percentile giving the costmin (or less)
  if (costthresh <= costmin)  costthresh = Max(costmin*1.0001,costmin*0.9999);
  if (globalopts.verbose>=4) {
    cout << "Cost threshold = " << costthresh << " and minimum is "
	 << costmin << endl;
  }

  // for all costs less than 150% of the best cost so far, optimise...
  //float costthresh = bestpts(bestidx,4) * 1.5;
  int numsubcost = 0;
  for (int ix=0; ix<costs.columns(); ix++) {
    for (int iy=0; iy<costs.rows(); iy++) {
      for (int iz=0; iz<costs.slices(); iz++) {
	if (costs(ix,iy,iz) <= costthresh) {
	  numsubcost++;
	}
      }
    }
  }
  if (numsubcost<=0) {
    cerr << "WARNING: Found 0 or less sub-threshold costs" << endl;
    numsubcost = 1;
  }
  Matrix bestparams(numsubcost,13);
  int n=1;
  for (int ix=0; ix<finerx.Nrows(); ix++) {
    for (int iy=0; iy<finery.Nrows(); iy++) {
      for (int iz=0; iz<finerz.Nrows(); iz++) {
	if (costs(ix,iy,iz) < costthresh) {
	  rx = finerx(ix+1);
	  ry = finery(iy+1);
	  rz = finerz(iz+1); 
	  xf = ((float) ix)*factorx;
	  yf = ((float) iy)*factory;
	  zf = ((float) iz)*factorz;
	  txv = tri_interpolation(tx,xf,yf,zf);
	  tyv = tri_interpolation(ty,xf,yf,zf);
	  tzv = tri_interpolation(tz,xf,yf,zf);
	  scv = tri_interpolation(scale,xf,yf,zf);
	  if ((scv<0.5) || (scv>2.0))  scv = 1.0;
	  if (globalopts.dof<=6) scv = 1.0;
	  params_8 = 0.0;
	  params_8(1) = rx;  params_8(2) = ry;  params_8(3) = rz;
	  params_8(4) = txv; params_8(5) = tyv; params_8(6) = tzv;
	  params_8(7) = scv; params_8(8) = scv; params_8(9) = scv;
	  globalopts.refparams = params_8;
	  params12toN(params_8);
	  optimise(params_8,globalopts.parammask.Ncols(),
		   param_tol,&no_its,&fans,subset_costfn);
	  paramsNto12(params_8);
	  costs(ix,iy,iz) = fans;
	  bestparams(n,1) = fans;
	  bestparams.SubMatrix(n,n,2,13) = params_8.t();
	  n++;
	  if (globalopts.verbose>=3) {
	    cout << "(" << ix << "," << iy << "," << iz << ") => " << fans
		 << " with " << params_8.t();
	  }
	}
      }
    }
  }

  if (globalopts.verbose>=3) {
    safe_save_volume(costs,"costs");
    cout << "Costs (1st column) are:\n" << bestparams << endl;
  }

  // find the single best cost
  /*
    bestidx = 1;
    for (int n=1; n<=bestparams.Nrows(); n++) {
    if (bestparams(n,1) < bestparams(bestidx,1)) {
    bestidx = n;
    }
    }

    params = bestparams.SubMatrix(bestidx,bestidx,2,13).t();
    if (globalopts.verbose>=3) {
    cout << "Chosen parameter: " << params.t() << endl;
    }
  */

  // find the cost minima and return these
  Matrix bestpts;
  find_cost_minima(bestpts,costs);
  paramlist.ReSize(bestpts.Nrows(),12);
  int ix,iy,iz;
  for (int n=1; n<=paramlist.Nrows(); n++) {
    ix = (int) bestpts(n,1);
    iy = (int) bestpts(n,2);
    iz = (int) bestpts(n,3);
    if (globalopts.verbose>=3) 
      cerr << "Cost minima at : " << ix << "," << iy << "," << iz << endl;
    rx = finerx(ix+1);
    ry = finery(iy+1);
    rz = finerz(iz+1);
    xf = ((float) ix)*factorx;
    yf = ((float) iy)*factory;
    zf = ((float) iz)*factorz;
    txv = tri_interpolation(tx,xf,yf,zf);
    tyv = tri_interpolation(ty,xf,yf,zf);
    tzv = tri_interpolation(tz,xf,yf,zf);
    scv = tri_interpolation(scale,xf,yf,zf);
    if ((scv<0.5) || (scv>2.0))  scv = 1.0;
    if (globalopts.dof<=6) scv = 1.0;
    params_8 = 0.0;
    params_8(1) = rx;  params_8(2) = ry;  params_8(3) = rz;
    params_8(4) = txv; params_8(5) = tyv; params_8(6) = tzv;
    params_8(7) = scv; params_8(8) = scv; params_8(9) = scv;
    paramlist.SubMatrix(n,n,1,12) = params_8.t();
  }

  if (globalopts.verbose>=3) {
    cout << "Chosen parameters:\n" << paramlist << endl;
  }

  globalopts.anglerep = useranglerep;
  globalopts.verbose = storedverbose;
  globalopts.dof = storeddof;
}


////////////////////////////////////////////////////////////////////////////

float measure_cost(Matrix& affmat, int input_dof)
{
  Tracer tr("measure_cost");
  // the most basic strategy - just do a single optimisation run at the
  //  specified dof
  int dof=input_dof;
  if (dof<6) { 
    cerr << "Erroneous dof " << dof << " : using 6 instead\n"; 
    dof=6; 
  }
  if (dof>12) {
    cerr << "Erroneous dof " << dof << " : using 12 instead\n"; 
    dof=12;
  }

  return costfn(affmat);
}  

////////////////////////////////////////////////////////////////////////////

int optimise_strategy1(Matrix& matresult, float& fans, int input_dof, 
		       int max_iterations=4)
{
  Tracer tr("optimise_strategy1");
  // the most basic strategy - just do a single optimisation run at the
  //  specified dof
  int dof=input_dof;
  if (dof<6) { 
    cerr << "Erroneous dof " << dof << " : using 6 instead\n"; 
    dof=6; 
  }
  if (dof>12) {
    cerr << "Erroneous dof " << dof << " : using 12 instead\n"; 
    dof=12;
  }

  ColumnVector params(12), param_tol(12);
  int no_its=0;
  globalopts.no_params = dof;
  set_param_tols(param_tol,12);  // 12 used to be dof
  param_tol = param_tol * estimate_scaling();
  affmat2vector(matresult,dof,params);
  optimise(params,dof,param_tol,&no_its,&fans,costfn,max_iterations);
  vector2affine(params,dof,matresult);
  return no_its;
}  

//-------------------------------------------------------------------------//

void optimise_strategy2(Matrix& matresult)
{
  Tracer tr("optimise_strategy2");
  // initially do a 7 dof optimisation, then increase to 9 and 12 
  //  (if called for)  - this is best used on the 2x subsampled image
  int dof, no_its=0;
  float costval=0.0;

  dof = Min(globalopts.dof,7);
  no_its = optimise_strategy1(matresult,costval,dof);

  if ((no_its<3) || (globalopts.dof>dof)) {
    dof = Min(globalopts.dof,9);
    optimise_strategy1(matresult,costval,dof,1);
  }

  if (globalopts.dof>=10) {
    dof = Min(globalopts.dof,12);
    optimise_strategy1(matresult,costval,dof,2);
  }
}

//-------------------------------------------------------------------------//

float rms_deviation(const Matrix& affmat1, const Matrix& affmat2) 
{
  float rmax = 80.0;
  Tracer trcr("rms_deviation");
  Matrix adiff(3,3);
  adiff = affmat1.SubMatrix(1,3,1,3) - affmat2.SubMatrix(1,3,1,3);
  ColumnVector tr(3);
  tr = affmat1.SubMatrix(1,3,4,4) - affmat2.SubMatrix(1,3,4,4);
  float rms = (tr.t() * tr).AsScalar() + rmax*rmax/5.0*Trace(adiff.t()*adiff);
  return rms;
}


int sorted_posn(const float costval, const Matrix& opt_matrixlist)
{
  Tracer tr("sorted_posn");
  int n=1;
  if (opt_matrixlist.Nrows()<1) { return 1; }
  for (n=1; n<=opt_matrixlist.Nrows(); n++) {
    if (costval < opt_matrixlist(n,1)) {
      return n;
    }
  }
  return opt_matrixlist.Nrows()+1;
}
 

void delete_row(Matrix& mat, const int row) 
{
  Tracer tr("delete_row");
  if ((row<1) || (row>mat.Nrows()))  return;
  int cols = mat.Ncols();
  if (cols<1) return;
  Matrix temp(mat.Nrows()-1,cols);
  if (row>=2) {
    temp.SubMatrix(1,row-1,1,cols) = mat.SubMatrix(1,row-1,1,cols);
  }
  if (row<mat.Nrows()) {
    temp.SubMatrix(row,mat.Nrows()-1,1,cols) 
      = mat.SubMatrix(row+1,mat.Nrows(),1,cols);
  }
  mat = temp;
}
 

void insert_row(Matrix& mat, const int row, const Matrix& rowmat)
{
  Tracer tr("insert_row");
  Matrix temp;
  if ((row<1) || (row>(mat.Nrows()+1)))  return;
  int cols = mat.Ncols(), nmax = mat.Nrows();
  if (cols<1) return;
  if ((rowmat.Ncols()!=cols) || (rowmat.Nrows()!=1)) { return; }
  temp.ReSize(nmax + 1,cols);
  if ((nmax>0) && (row>=2)) {
    temp.SubMatrix(1,row-1,1,cols) = mat.SubMatrix(1,row-1,1,cols);
  }
  if ((nmax>0) && (row<=nmax)) {
    temp.SubMatrix(row+1,nmax+1,1,cols) = mat.SubMatrix(row,nmax,1,cols);
  }
  temp.SubMatrix(row,row,1,cols) = rowmat;
  mat = temp;
}


int add2list(const Matrix& rowmat, Matrix& opt_matrixlist, float rms_min=1.0) 
{
  Tracer tr("add2list");
  if (rowmat.Nrows()!=1) {
    cerr << "WARNING (add2list): cannot add matrix with " << rowmat.Nrows() 
	 << " to matrix\n";
    return -1;
  }
  Matrix mat1(4,4), mat2(4,4);
  reshape(mat1,rowmat.SubMatrix(1,1,2,17),4,4);
  float costval = rowmat(1,1);
  int nmax = opt_matrixlist.Nrows();
  for (int n=1; n<=nmax; n++) {
    reshape(mat2,opt_matrixlist.SubMatrix(n,n,2,17),4,4);
    if (rms_deviation(mat1,mat2)<rms_min) {
      if (costval < opt_matrixlist(n,1)) {
	// if it is close but better, then delete, and add in appropriate posn
	delete_row(opt_matrixlist,n);
	break;
      } else {
	// it is close but worse, so don't add it
	return -1;
      }
    }
  }
  // otherwise add this one to maintain sorted order (ascending) of cost
  int pos = sorted_posn(costval,opt_matrixlist);
  insert_row(opt_matrixlist,pos,rowmat);
  return pos;
}


void optimise_strategy3(Matrix& opt_matrixlist)
{
  Tracer tr("optimise_strategy3");
  // the basic low-level search strategy - it searches the cost function
  //  space (across rotations)  then selects the best candidates and optimises
  //  these too   : best used for the 8x subsampled (far too slow otherwise)
  // the returned matrixlist contains in each row the cost and reshaped matrix 
  //  for both optimised, and pre-optimised positions, ranked in ascending
  //  order of the OPTIMISED cost
  
  
  Matrix paramlist;
  volume costs,tx,ty,tz,scale;
  globalopts.currentcostfn = globalopts.searchcostfn;
  search_cost(paramlist,costs,tx,ty,tz,scale);
  
  int dof = Min(globalopts.dof,7);
  ColumnVector params_8(12);
  float costval=0.0, rms_min = estimate_scaling();
  Matrix reshapedmat(1,16), matresult(4,4), premat(4,4), matrow(1,34);
  opt_matrixlist.ReSize(0,34);
  // freely optimise each member of the parameter list (allow rotns to vary)
  int verbose = globalopts.verbose;
  globalopts.verbose-=2;
  if (verbose>=3) { 
    cout << "After free optimisation, parameters are:" << endl;
  }
  for (int n=1; n<=paramlist.Nrows(); n++) {
    params_8 = paramlist.SubMatrix(n,n,1,12).t();
    vector2affine(params_8,dof,matresult);
    premat = matresult;
    optimise_strategy1(matresult,costval,dof);
    // form the optimised and pre-optimised costs and matrices into a row
    reshape(reshapedmat,matresult,1,16);
    matrow(1,1) = costval;
    matrow.SubMatrix(1,1,2,17) = reshapedmat;
    matrow(1,18) = costfn(premat);
    reshape(reshapedmat,premat,1,16);
    matrow.SubMatrix(1,1,19,34) = reshapedmat;
    // add the row to the optimised list (ordered by optimised cost)
    add2list(matrow,opt_matrixlist,rms_min);
    if (verbose>=3) { 
      affmat2vector(matresult,12,params_8);
      cout << costval << " :: " << params_8.t(); 
    }
  }
  globalopts.verbose = verbose;
  globalopts.currentcostfn = globalopts.maincostfn;

  if (globalopts.verbose>=3) {
    cout << "Parameters (1st column costs) are:\n" << opt_matrixlist << endl;
  }
}

//-------------------------------------------------------------------------//

void set_perturbations(Matrix& delta, Matrix& perturbmask)
{
  Tracer tr("set_perturbations");
  // set the perturbations to be applied to the pre-optimised
  //  solutions (used in optimise_strategy4)
  // the perturbations are given by the elementwise product of
  //  delta (a single column) and the individual columns of perturbmask
  // the number of columns of perturbmask set the number of perturbations
  //  to be tried - it should always include a zero column = unperturbed case
  ColumnVector coarserx, coarsery, coarserz, finerx, finery, finerz;
  set_rot_samplings(coarserx,coarsery,coarserz,finerx,finery,finerz);
  ColumnVector param_tol(12);
  set_param_tols(param_tol,12);
  param_tol = param_tol * estimate_scaling();
  // set the magnitude of the variations
  float delscale, delrx, delry, delrz;
  delrx = 3.0*param_tol(1); 
  delry = 3.0*param_tol(2);
  delrz = 3.0*param_tol(3);
  if (finerx.Nrows()>1) { delrx = 0.5*(finerx(2) - finerx(1)); }
  if (finery.Nrows()>1) { delry = 0.5*(finery(2) - finery(1)); }
  if (finerz.Nrows()>1) { delrz = 0.5*(finerz(2) - finerz(1)); }
  delscale=0.0;  
  if (globalopts.dof>=7) delscale = 0.1;  // scale only set if allowed
  delta.ReSize(12,1);
  delta = 0.05;  // translation default
  delta(1,1) = delrx; delta(2,1) = delry;  delta(3,1) = delrz;
  delta(7,1) = delscale;  delta(8,1) = delscale;  delta(9,1) = delscale;
  // set the directions of the variations 
  // (+/- for each rotation and  4 scale settings: +/- 1 and 2 times delscale)
  perturbmask.ReSize(12,11);
  perturbmask = 0.0;
  perturbmask(1,1) = 1.0;
  perturbmask(1,2) = -1.0;
  perturbmask(2,3) = 1.0;
  perturbmask(2,4) = -1.0;
  perturbmask(3,5) = 1.0;
  perturbmask(3,6) = -1.0;
  perturbmask(7,7) = 1.0; 
  perturbmask(8,7) = 1.0; 
  perturbmask(9,7) = 1.0;
  perturbmask(7,8) = 2.0; 
  perturbmask(8,8) = 2.0; 
  perturbmask(9,8) = 2.0;
  perturbmask(7,9) = -1.0; 
  perturbmask(8,9) = -1.0; 
  perturbmask(9,9) = -1.0;
  perturbmask(7,10)=-2.0; 
  perturbmask(8,10)=-2.0; 
  perturbmask(9,10)=-2.0;
}


void optimise_strategy4(Matrix& matresult, const Matrix& opt_matrixlist,
			int max_num = 3)
{
  Tracer tr("optimise_strategy4");
  // takes the output of optimise_strategy4 and reoptimises the
  //  first max_num of the optimised results, including the pre-optimised
  //  positions and perturbations of these - defined in set_perturbations()
  int dof = Min(globalopts.dof,7);
  float bestcost=0.0, costval=0.0;
  Matrix bestmat(4,4), reshaped(1,16);
  ColumnVector params_8(12);
  
  Matrix delta, perturbmask;
  set_perturbations(delta,perturbmask);

  // for each potentially best transform do an optimisation
  int verbose = globalopts.verbose;
  globalopts.verbose -= 2;
  bool no_previous_cost = true;
  int no_its = Min(max_num,opt_matrixlist.Nrows());

  for (int n=1; n<=no_its; n++) {
    for (int matno=1; matno<=2; matno++) {
      for (int pert=1; pert<=perturbmask.Ncols(); pert++) {

	if ((matno==1) && (pert>1)) {
	  break;
	} else {
	  if (matno==1) { // the previously optimised case
	    reshaped = opt_matrixlist.SubMatrix(n,n,2,17);
	    reshape(matresult,reshaped,4,4);
	  } else { // the pre-optimised case with perturbations
	    reshaped = opt_matrixlist.SubMatrix(n,n,19,34);
	    reshape(matresult,reshaped,4,4);
	    affmat2vector(matresult,12,params_8);
	    // use the elementwise product to produce the perturbation
	    params_8 += SP(perturbmask.SubMatrix(1,12,pert,pert),delta);
	    vector2affine(params_8,12,matresult);
	    reshape(reshaped,matresult,1,16);
	  }
	}
	
	optimise_strategy1(matresult,costval,dof);
	if ( (no_previous_cost) || (costval <= bestcost) ) {
	  bestcost = costval;
	  bestmat = matresult;
	  no_previous_cost = false;
	}

	if (globalopts.verbose>=3) {
	  affmat2vector(matresult,dof,params_8);
	  cout << "Current parameters [" 
	       << (n-1)*(perturbmask.Ncols()+1) + pert + matno - 1
	       << " out of " << no_its*(perturbmask.Ncols()+1)
	       << "] (cost vs bestcost = " << costval 
	       << " vs " << bestcost << ") are:\n" << params_8.t();
	  reshape(matresult,reshaped,4,4);
	  affmat2vector(matresult,dof,params_8);
	  cout << "  from : " << params_8.t();
	}
      }
    }
  }

  globalopts.verbose = verbose;
  matresult = bestmat;  // the return value
}


////////////////////////////////////////////////////////////////////////////


int get_testvol(volume& testvol)
{
  Tracer tr("get_testvol");
  read_volume(testvol,globalopts.inputfname);
  read_matrix(globalopts.initmat,globalopts.initmatfname,testvol);
  set_default_planes(testvol);

  // THE FOLLOWING IS OBSOLETE
  // use the affmat_init as an initial transform estimate for testvol
  // if (globalopts.verbose >= 1) 
  //    cout << "Applying initial transform" << endl;
  //   check to see if there is any transform to do...
  //  Matrix diffmat(4,4);
  //  identity(diffmat);
  //  diffmat = diffmat - globalopts.initmat;
  //  if (diffmat.MaximumAbsoluteValue() > 1e-5) {
  //   volume newtestvol(testvol);
  //   affine_transform(newtestvol,testvol,globalopts.initmat);
  //   testvol = newtestvol;
  //  }
  if (globalopts.verbose>=2) {
    cerr << "Testvol sampling matrix =\n" << testvol.sampling_matrix() << endl;
  }
  // calculate the plane fits
  fit_planes(testvol,globalopts.planetestfname);
  //volume validmask;
  //testvol.create_valid_mask(validmask);
  //safe_save_volume(validmask,"validmask");
  //safe_save_volume(testvol,"t2tst");
  return 0;
}  


int get_refvol(volume& refvol)
{
  Tracer tr("get_refvol");
  read_volume(refvol,globalopts.reffname);
  set_default_planes(refvol);
  return 0;
}


int resample_refvol(volume& refvol, float sampling=1.0)
{
  Tracer tr("resample_refvol");
  float sampl=sampling;
  if (sampl<1e-4) { 
    cerr << "WARNING: sampling " << sampl << " less than 0.0001" << endl
	 << "         Setting to 1.0 sampling instead" << endl;
    sampl=1.0;
  }
  if (sampl<0.1) { 
    cerr << "WARNING: under minimum sampling of 0.1" << endl;
  }
  if (globalopts.verbose >= 1) 
    cout << "Resampling refvol isotropically" << endl;
  // isotropically resample the volume
  volume tmpvol;
  blur4subsampling(tmpvol,refvol,sampl);
  isotropic_rescale(refvol,tmpvol,sampl);
  if (globalopts.verbose>=2) print_volume_info(refvol,"Refvol");      
  
  // calculate the plane fits
  fit_planes(refvol,globalopts.planereffname);
  return 0;
}


////////////////////////////////////////////////////////////////////////////


void no_optimise()
{
  Tracer tr("no_optimise");
  volume refvol, testvol;
  // set up image pair and global pointer
  
  get_refvol(refvol);
  read_volume(testvol,globalopts.inputfname);
  read_matrix(globalopts.initmat,globalopts.initmatfname,testvol);
  set_default_planes(testvol);
  fit_planes(testvol,globalopts.planetestfname);
  
  float min_sampling_ref=1.0, min_sampling_test=1.0;
  min_sampling_ref = Min(refvol.getx(),Min(refvol.gety(),refvol.getz()));
  min_sampling_test = Min(testvol.getx(),Min(testvol.gety(),testvol.getz()));
  globalopts.min_sampling = Max(min_sampling_ref,min_sampling_test);
  
  if (globalopts.iso) {
    resample_refvol(refvol,globalopts.min_sampling);
  }
  {
    volume testvol_1;
    blur4subsampling(testvol_1,testvol,globalopts.min_sampling);
    testvol = testvol_1;
  }

  imagepair global_1(refvol,testvol);
  global_1.set_no_bins(globalopts.no_bins);
  globalopts.impair = &global_1;
  
  if (globalopts.measure_cost) {
    float cost = costfn(globalopts.initmat);
    cout << "\n\nCost = " << cost << "\n\n";
  } else {
    volume outputvol = globalopts.impair->refvol;
    filled_affine_transform(outputvol,globalopts.impair->testvol,globalopts.initmat);
    safe_save_volume(outputvol,globalopts.outputfname.c_str());
    if (globalopts.verbose>=2) {
      save_matrix_data(globalopts.initmat,globalopts.impair->testvol,
		       outputvol);
    }
  }
  exit(0);
}


////////////////////////////////////////////////////////////////////////////


int firstelementless(const RowVector& r1, const RowVector& r2)
{
  Tracer tr("firstelementless");
  return (r1(1) < r2(1));
}


void stripleadingspace(string& line)
{
  Tracer tr("stripleadingspace");
  if (line.length()<1) return;
  string::size_type i2=line.find_first_not_of("	 ");
  if (i2==string::npos) {
    return;
  }
  if (i2!=0) {
    line.erase(0,i2);
  }
}


int parseline(const string& inname, std::vector<string>& words) 
{
  Tracer tr("parseline");
  string name = inname;
  string::size_type i1;
  words.clear();
  while (name.length()>=1) {
    stripleadingspace(name);
    i1=name.find_first_of("	 ");
    words.push_back(name.substr(0,i1));
    name.erase(0,i1);
  }
  return 0;
}



int setmatvariable(const string& name, MatVecPtr& namedmat)
{
  Tracer tr("setmatvariable");
  if (name.size()<1) return -2;
  if (name == "S") {
    namedmat = &globalsearchoptmat;
    return 0;
  } else if (name == "P") {
    namedmat = &globalpreoptsearchmat;
    return 0;
  } else if (name == "U") {
    namedmat = &(globalusrmat[0]);
    return 0;
  } else if (name[0] == 'U') {
    int idx = name[1] - 'A' + 1;
    if ((idx<1) || (idx>27)) {
      if ((name[1]<'0') || (name[1]>'9')) {
	cerr << "Can only accept UA to UZ, not " << name << endl;
	return -3;
      } else {
	idx = 0;
      }
    }
    namedmat = &(globalusrmat[idx]);
    return 0;
  } else {
    cerr << "Cannot interpret " << name << " as a Matrix name" << endl;
    return -1;
  }
}


int setscalarvariable(const string& name, float& namedscalar)
{
  Tracer tr("setscalarvariable");
  if (name.size()<1) return -2;
  if (name == "MAXDOF") {
    namedscalar = (float) globalopts.dof;
    return 0;
  } else if (name == "MINSAMPLING") {
    namedscalar = (float) globalopts.min_sampling;
    return 0;
  } else if (isalpha(name[0])) {
    cerr << "Cannot interpret " << name << " as a scalar variable name" << endl;
    namedscalar = 0.0;
    return -1;
  } else {
    namedscalar = atof(name.c_str());
  }
  return 0;
}


int setscalarvariable(const string& name, int& namedscalar)
{
  Tracer tr("setscalarvariable");
  float temp=0.0;
  int retval = setscalarvariable(name,temp);
  if (temp>0.0) {
    namedscalar = (int) (temp+0.5);
  } else {
    namedscalar = (int) (temp-0.5);
  }
  return retval;
}


int parsematname(const string& inname, MatVecPtr& usrdefmat, int& r1, int& r2) 
{
  Tracer tr("parsematname");
  string name = inname;
  stripleadingspace(name);
  string::size_type colon, hyphen;
  if (name.length()<1) {
    return 0;
  }
  colon=name.find_first_of(":");
  hyphen=name.find_first_of("-");
  string basename=name, row1="", row2="";
  if (colon!=string::npos) {
    basename = name.substr(0,colon);
  }
  if ((colon!=string::npos) && (hyphen!=string::npos)) {
    // there is a colon and a hyphen
    row1 = name.substr(colon+1,hyphen-colon-1);
    row2 = name.substr(hyphen+1);
  } else {
    if ((colon!=string::npos) && (hyphen==string::npos)) {
      // there is a colon but no hyphen
      row1 = name.substr(colon+1);
      row2 = row1;
    }
  }
  if (row1.length()<1) row1="1";
  if (row2.length()<1)  row2="999999";

  setmatvariable(basename, usrdefmat);
  setscalarvariable(row1,r1);
  setscalarvariable(row2,r2);
  if (r1<1) r1=1;
  if (r2<r1) r2=r1;
  return 0;
}


//////////////////////////////////////////////////////////////////////////////

void usrcopy(MatVecPtr usrsrcmat, MatVecPtr usrdestmat) 
{
  Tracer tr("usrcopy");
  // COPY src -> dest
  usrdestmat->erase(usrdestmat->begin(),usrdestmat->end());
  for (unsigned int i=0; i<usrsrcmat->size(); i++) {
    usrdestmat->push_back((*usrsrcmat)[i]);
  }
}


void usrclear(MatVecPtr usrsrcmat) 
{
  Tracer tr("usrclear");
  // CLEAR src
  usrsrcmat->erase(usrsrcmat->begin(),usrsrcmat->end());
}


void usrprint(MatVecPtr usrsrcmat, unsigned int row1, unsigned int row2) 
{
  Tracer tr("usrprint");
  // PRINT src
  for (unsigned int r=row1; r<=Min(usrsrcmat->size(),row2); r++) {
    cout << (*usrsrcmat)[r-1];
  }
}


int usrsave(string filename, MatVecPtr usrsrcmat, 
	     unsigned int row1, unsigned int row2) 
{
  Tracer tr("usrsave");
  // SAVE src

  ofstream fptr(filename.c_str());
  if (!fptr) { 
    cerr << "Could not open file " << filename << " for writing" << endl;
    return -1;
  }
  for (unsigned int r=row1; r<=Min(usrsrcmat->size(),row2); r++) {
    fptr << (*usrsrcmat)[r-1];
  }
  fptr << endl;
  fptr.close();
  return 0;
}


int usrread(string filename, MatVecPtr usrsrcmat)
{
  Tracer tr("usrread");
  // READ src

  ifstream fptr(filename.c_str());
  if (!fptr) { 
    cerr << "Could not open file " << filename << " for reading" << endl;
    return -1;
  }
  usrsrcmat->clear();
  RowVector currow(17);
  while (!fptr.eof()) {
    currow = 0.0;
    for (unsigned int c=1; c<=17; c++) {
      fptr >> currow(c);
    }
    if (!fptr.eof()) {
      usrsrcmat->push_back(currow);
    }
  }
  fptr.close();
  return 0;
}


int usrsetrow(MatVecPtr usrsrcmat,const std::vector<string> &words)
{
  Tracer tr("usrsetrow");
  // SETROW src

  RowVector currow(17);
  currow = 0.0;
  for (unsigned int c=2; c<=17; c++) {
    currow(c) = atof(words[c].c_str());
  }
  usrsrcmat->push_back(currow);
  return 0;
}


void usrsort(MatVecPtr usrsrcmat) 
{
  Tracer tr("usrsort");
  // SORT src
  sort(usrsrcmat->begin(),usrsrcmat->end(),firstelementless);
}


void usrsearch(MatVecPtr searchoptmat, MatVecPtr preoptsearchmat, int sdof=12) 
{
  Tracer tr("usrsearch");
  // SEARCH
  Matrix opt_matrixlist;
  globalopts.searchdof = sdof;
  optimise_strategy3(opt_matrixlist);
  MatVec optmatsorted(0);
  RowVector tmprow;
  for (int i=1; i<=opt_matrixlist.Nrows(); i++) {
    tmprow = opt_matrixlist.SubMatrix(i,i,1,34);
    optmatsorted.push_back(tmprow);
  }
  sort(optmatsorted.begin(),optmatsorted.end(),firstelementless);
  for (unsigned int i=0; i<optmatsorted.size(); i++) {
    tmprow = optmatsorted[i].SubMatrix(1,1,1,17);
    searchoptmat->push_back(tmprow);
    tmprow = optmatsorted[i].SubMatrix(1,1,18,34);
    preoptsearchmat->push_back(tmprow);
  }
}


int usrreadparams(string filename, MatVecPtr usrsrcmat)
{
  Tracer tr("usrreadparams");
  // READ src

  ifstream fptr(filename.c_str());
  if (!fptr) { 
    cerr << "Could not open file " << filename << " for reading" << endl;
    return -1;
  }
  usrsrcmat->clear();
  ColumnVector params(12);
  Matrix affmat(4,4), reshaped(1,16);
  RowVector rowresult(17);
  float costval;
  while (!fptr.eof()) {
    params = 0.0;
    fptr >> costval;
    for (unsigned int c=1; c<=12; c++) {
      fptr >> params(c);
    }
    if (!fptr.eof()) {
      vector2affine(params,12,affmat);
      reshape(reshaped,affmat,1,16);
      rowresult(1) = costval;
      rowresult.SubMatrix(1,1,2,17) = reshaped;
      usrsrcmat->push_back(rowresult);
    }
  }
  fptr.close();
  return 0;
}


int usrsaveparams(string filename, MatVecPtr usrmatptr, 
		    unsigned int usrrow1, unsigned int usrrow2)
{
  Tracer tr("usrsaveparams");
  // PRINTPARAMS
  ofstream fptr(filename.c_str());
  if (!fptr) { 
    cerr << "Could not open file " << filename << " for writing" << endl;
    return -1;
  }
  Matrix matresult;
  ColumnVector params(12);
  for (unsigned int crow=usrrow1; crow<=Min(usrrow2,usrmatptr->size()); crow++)
    {
      // the pre-optimised case with perturbations
      Matrix reshaped = (*usrmatptr)[crow-1].SubMatrix(1,1,2,17);
      reshape(matresult,reshaped,4,4);
      affmat2vector(matresult,12,params);
      fptr << ((*usrmatptr)[crow-1])(1) << " " << params.t() << endl;
    }
  fptr << endl;
  fptr.close();
  return 0;
}


void usrprintparams(MatVecPtr usrmatptr, 
		    unsigned int usrrow1, unsigned int usrrow2)
{
  Tracer tr("usrprintparams");
  // PRINTPARAMS
  Matrix matresult;
  ColumnVector params(12);
  for (unsigned int crow=usrrow1; crow<=Min(usrrow2,usrmatptr->size()); crow++)
    {
      // the pre-optimised case with perturbations
      Matrix reshaped = (*usrmatptr)[crow-1].SubMatrix(1,1,2,17);
      reshape(matresult,reshaped,4,4);
      affmat2vector(matresult,12,params);
      cout << ((*usrmatptr)[crow-1])(1) << " " << params.t() << endl;
    }
}


void usrmeasurecost(MatVecPtr stdresultmat, 
		    MatVecPtr usrmatptr, 
		    unsigned int usrrow1, unsigned int usrrow2, int usrdof, 
		    ColumnVector& usrperturbation, bool usrperturbrelative)
{
  Tracer tr("usrmeasurecost");
  // OPTIMISE
  Matrix delta, perturbmask, matresult;
  set_perturbations(delta,perturbmask);
  int dof = Min(globalopts.dof,usrdof);
  ColumnVector params(12);
  RowVector rowresult(17);
  for (unsigned int crow=usrrow1; crow<=Min(usrrow2,usrmatptr->size()); crow++)
    {
      // the pre-optimised case with perturbations
      Matrix reshaped = (*usrmatptr)[crow-1].SubMatrix(1,1,2,17);
      reshape(matresult,reshaped,4,4);
      affmat2vector(matresult,12,params);
      // use the elementwise product to produce the perturbation
      if (usrperturbrelative) {
	params += SP(usrperturbation,delta); // rel
      } else {
	params += usrperturbation; // abs
      }
      vector2affine(params,12,matresult);
      
      float costval=0.0;
      costval = measure_cost(matresult,dof);
      reshape(reshaped,matresult,1,16);
      rowresult(1) = costval;
      rowresult.SubMatrix(1,1,2,17) = reshaped;
      // store result
      stdresultmat->push_back(rowresult);
    }
}


void usroptimise(MatVecPtr stdresultmat, 
		 MatVecPtr usrmatptr, 
		 unsigned int usrrow1, unsigned int usrrow2, int usrdof, 
		 ColumnVector& usrperturbation, bool usrperturbrelative,
		 int usrmaxitn)
{
  Tracer tr("usroptimise");
  // OPTIMISE
  Matrix delta, perturbmask, matresult;
  set_perturbations(delta,perturbmask);
  int dof = Min(globalopts.dof,usrdof);
  ColumnVector params(12);
  RowVector rowresult(17);
  for (unsigned int crow=usrrow1; crow<=Min(usrrow2,usrmatptr->size()); crow++)
    {
      // the pre-optimised case with perturbations
      Matrix reshaped = (*usrmatptr)[crow-1].SubMatrix(1,1,2,17);
      reshape(matresult,reshaped,4,4);
      affmat2vector(matresult,12,params);
      // use the elementwise product to produce the perturbation
      if (usrperturbrelative) {
	params += SP(usrperturbation,delta); // rel
      } else {
	params += usrperturbation; // abs
      }
      vector2affine(params,12,matresult);
      
      float costval=0.0;
      optimise_strategy1(matresult,costval,dof,usrmaxitn);
      reshape(reshaped,matresult,1,16);
      rowresult(1) = costval;
      rowresult.SubMatrix(1,1,2,17) = reshaped;
      // store result
      stdresultmat->push_back(rowresult);
    }
}


void usrsetscale(int usrscale,
		 volume& testvol, volume& refvol, 
		 volume& refvol_2, volume& refvol_4, 
		 volume& refvol_8) {
  // SETSCALE (int usrscale = 8,4,2,1)
  // testvol must be passed in as a static storage is needed so that
  //  the object pointed to in globalopts.impair does not go out of scope
  Tracer tr("usrsetscale");
  int scale = usrscale;
  imagepair *globalpair=0;
  if (globalopts.min_sampling<=1.25 * ((float) scale)) {
    get_testvol(testvol);
    volume testvolnew;
    blur4subsampling(testvolnew,testvol,(float) scale);
    testvol = testvolnew;
    volume *refvolnew=0;
    if (scale==8) {
      refvolnew = &refvol_8;
    } else if (scale==4) {
      refvolnew = &refvol_4;
    } else if (scale==2) {
      refvolnew = &refvol_2;
    } else if (scale==1) {
      refvolnew = &refvol;
    } else {
      cerr << "Cannot set scale to " << scale << endl;
      cerr << "Instead using unity scale" << endl;
      scale = 1;
      refvolnew = &refvol;
    }
    globalpair = new imagepair(*refvolnew,testvol);
    globalpair->set_no_bins(globalopts.no_bins/scale);
    globalopts.impair = globalpair;
  }
}


void interpretcommand(const string& comline, bool& skip,
		      volume& testvol, volume& refvol, 
		      volume& refvol_2, volume& refvol_4, 
		      volume& refvol_8)
{
  Tracer tr("interpretcommand");
  std::vector<string> words(0);
  parseline(comline,words);
  if (words.size()<1) return;
  if ((words[0])[0] == '#') return;  // comment line

  if (skip) {
    skip=false;
    return;
  }

  if (words[0]=="copy") {
    // COPY
    if (words.size()<3) {
      cerr << "Wrong number of args to COPY" << endl;
      exit(-1);
    }
    MatVecPtr src, dest;
    int d1, d2, d3, d4;
    parsematname(words[1],src,d1,d2);
    parsematname(words[2],dest,d3,d4);
    usrcopy(src,dest);
  } else if (words[0]=="clear") {
    // CLEAR
    if (words.size()<2) {
      cerr << "Wrong number of args to CLEAR" << endl;
      exit(-1);
    }
    MatVecPtr src;
    int d1, d2;
    parsematname(words[1],src,d1,d2);
    usrclear(src);
  } else if (words[0]=="print") {
    // PRINT
    if (words.size()<2) {
      cerr << "Wrong number of args to PRINT" << endl;
      exit(-1);
    }
    MatVecPtr src;
    int d1, d2;
    parsematname(words[1],src,d1,d2);
    usrprint(src,d1,d2);
  } else if (words[0]=="save") {
    // SAVE
    if (words.size()<3) {
      cerr << "Wrong number of args to SAVE" << endl;
      exit(-1);
    }
    MatVecPtr src;
    int d1, d2;
    parsematname(words[1],src,d1,d2);
    usrsave(words[2],src,d1,d2);
  } else if (words[0]=="read") {
    // READ
    if (words.size()<3) {
      cerr << "Wrong number of args to READ" << endl;
      exit(-1);
    }
    MatVecPtr src;
    int d1, d2;
    parsematname(words[1],src,d1,d2);
    usrread(words[2],src);
  } else if (words[0]=="printparams") {
    // PRINTPARAMS
    if (words.size()<2) {
      cerr << "Wrong number of args to PRINTPARAMS" << endl;
      exit(-1);
    }
    MatVecPtr src;
    int d1, d2;
    parsematname(words[1],src,d1,d2);
    usrprintparams(src,d1,d2);
  } else if (words[0]=="saveparams") {
    // SAVEPARAMS
    if (words.size()<3) {
      cerr << "Wrong number of args to SAVEPARAMS" << endl;
      exit(-1);
    }
    MatVecPtr src;
    int d1, d2;
    parsematname(words[1],src,d1,d2);
    usrsaveparams(words[2],src,d1,d2);
  } else if (words[0]=="readparams") {
    // READPARAMS
    if (words.size()<3) {
      cerr << "Wrong number of args to READPARAMS" << endl;
      exit(-1);
    }
    MatVecPtr src;
    int d1, d2;
    parsematname(words[1],src,d1,d2);
    usrreadparams(words[2],src);
  } else if (words[0]=="sort") {
    // SORT
    if (words.size()<2) {
      cerr << "Wrong number of args to SORT" << endl;
      exit(-1);
    }
    MatVecPtr src;
    int d1, d2;
    parsematname(words[1],src,d1,d2);
    usrsort(src);
  } else if (words[0]=="search") {
    // SEARCH
    int usrdof=12;
    if (words.size()>2) {
      cerr << "Wrong number of args to SEARCH" << endl;
      exit(-1);
    }
    if (words.size()==2) {
      setscalarvariable(words[1],usrdof);
    }
    usrsearch(&globalsearchoptmat,&globalpreoptsearchmat,usrdof);
  } else if (words[0]=="optimise") {
    // OPTIMISE
    if (words.size()<5) {
      cerr << "Wrong number of args to OPTIMISE" << endl;
      exit(-1);
    }
    int usrdof=12, usrmaxitn=4, usrrow1=1, usrrow2=999999;
    ColumnVector usrperturbation(12);
    usrperturbation = 0.0;
    MatVecPtr usrdefmat;
    bool usrperturbrelative = true;
    setscalarvariable(words[1],usrdof);
    parsematname(words[2],usrdefmat,usrrow1,usrrow2);
    unsigned int wordno=3;
    float tempparam=0.0;
    while ((wordno<words.size()) 
	   && (words[wordno]!="rel") && (words[wordno]!="abs")) {
      setscalarvariable(words[wordno],tempparam);
      usrperturbation(wordno-2) = tempparam;
      wordno++;
    }
    if (wordno<words.size()) {
      if (words[wordno]=="abs")  usrperturbrelative = false;
      wordno++;
    }
    if (wordno<words.size()) {
      setscalarvariable(words[wordno],usrmaxitn);
    }
    usroptimise(&(globalusrmat[0]),usrdefmat,usrrow1,usrrow2,usrdof, 
		usrperturbation,usrperturbrelative,usrmaxitn);
  } else if (words[0]=="measurecost") {
    // MEASURECOST
    if (words.size()<5) {
      cerr << "Wrong number of args to MEASURECOST" << endl;
      exit(-1);
    }
    int usrdof=12, usrrow1=1, usrrow2=999999;
    ColumnVector usrperturbation(12);
    usrperturbation = 0.0;
    MatVecPtr usrdefmat;
    bool usrperturbrelative = true;
    setscalarvariable(words[1],usrdof);
    parsematname(words[2],usrdefmat,usrrow1,usrrow2);
    unsigned int wordno=3;
    float tempparam=0.0;
    while ((wordno<words.size()) 
	   && (words[wordno]!="rel") && (words[wordno]!="abs")) {
      setscalarvariable(words[wordno],tempparam);
      usrperturbation(wordno-2) = tempparam;
      wordno++;
    }
    if (wordno<words.size()) {
      if (words[wordno]=="abs")  usrperturbrelative = false;
      wordno++;
    }
    usrmeasurecost(&(globalusrmat[0]),usrdefmat,usrrow1,usrrow2,usrdof, 
		usrperturbation,usrperturbrelative);
  } else if (words[0]=="setscale") {
    // SETSCALE
    if (words.size()<2) {
      cerr << "Wrong number of args to SETSCALE" << endl;
      exit(-1);
    }
    int usrscale;
    setscalarvariable(words[1],usrscale);
    usrsetscale(usrscale,testvol,refvol,refvol_2,refvol_4,refvol_8);
  } else if (words[0]=="setrow") {
    // SETROW
    if (words.size()<18) {
      cerr << "Wrong number of args to SETROW" << endl;
      exit(-1);
    }
    MatVecPtr src;
    int d1, d2;
    parsematname(words[1],src,d1,d2);
    usrsetrow(src,words);
  } else if (words[0]=="if") {
    // IF
    if (words.size()<4) {
      cerr << "Wrong number of args to IF" << endl;
      exit(-1);
    }
    float arg1, arg2;
    setscalarvariable(words[1],arg1);
    setscalarvariable(words[3],arg2);
    if (words[2]==">") {
      skip = !(arg1 > arg2);
    } else if (words[2]=="<") {
      skip = !(arg1 < arg2);
    } else if (words[2]=="==") {
      skip = !(arg1 == arg2);
    } else if (words[2]=="!=") {
      skip = !(arg1 != arg2);
    } else if (words[2]=="<=") {
      skip = !(arg1 <= arg2);
    } else if (words[2]==">=") {
      skip = !(arg1 >= arg2);
    } else {
      cerr << "Cannot recognise operator " << words[2] << " in IF statement\n";
      exit(-1);
    }
  } else {
    cerr << "Unrecognised command " << words[0] << endl;
    cerr << "Quitting." << endl;
    exit(-1);
  }
}


////////////////////////////////////////////////////////////////////////////

// For random tests
void test_proc(void) {
  ColumnVector angl(3), centre(3);
  Matrix rot(4,4);
  float anglx=0.0;
  cout << "Input x angl: ";
  cin >> anglx;
  angl << anglx << 0.0 << 0.0 ;
  centre << 118.575 << 118.575 << 72.5 ;
  make_rot(angl,centre,rot);
  cout << "Rot matrix is: " << endl;
  cout << rot;
  exit(-1);
  return;
}

////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[])
{
  Tracer tr("main");
  try {

    //test_proc();  // only used for some debugging

    parse_command_line(argc, argv);

    if (!globalopts.do_optimise) {   
      no_optimise();
    }
    
    // read in the volumes
    volume refvol, testvol;
    get_refvol(refvol);
    get_testvol(testvol);

    float min_sampling_ref=1.0, min_sampling_test=1.0;
    min_sampling_ref = Min(refvol.getx(),Min(refvol.gety(),refvol.getz()));
    min_sampling_test = Min(testvol.getx(),Min(testvol.gety(),testvol.getz()));
    globalopts.min_sampling = Max(min_sampling_ref,min_sampling_test);

    
    if (globalopts.verbose>=3) {
      cout << "CoG for refvol is:  " << refvol.cog().t();
      cout << "CoG for testvol is:  " << testvol.cog().t();
    }
    
     // set up subsampled volumes by factors of 2, 4 and 8
    if (globalopts.verbose >= 1) 
      cout << "Subsampling the volumes" << endl;
    resample_refvol(refvol,1.0);
    volume refvol_2, refvol_4, refvol_8;
    subsample_by2(refvol,refvol_2);
    subsample_by2(refvol_2,refvol_4);
    subsample_by2(refvol_4,refvol_8);


    // set up image pair and global pointer
    {
      volume testvol_8;
      blur4subsampling(testvol_8,testvol,8.0);
      testvol = testvol_8;
    }
    imagepair global_8(refvol_8,testvol);
    global_8.set_no_bins(globalopts.no_bins/8);
    if (globalopts.verbose>=2) print_volume_info(testvol,"TESTVOL");

    globalopts.impair = &global_8;


    Matrix matresult(4,4);
    ColumnVector params_8(12), param_tol(12);
    float costval=0.0;

    if (globalopts.measure_cost) {
      identity(matresult);
      costval = costfn(matresult);
      cout << "\n\nCost = " << costval << "\n\n";
      return 1;
    }

    // first do the cost plot (it returns directly if not requested)
    costplot();

    // perform the optimisation

    globalopts.currentcostfn = globalopts.maincostfn;
    globalsearchoptmat.clear();
    globalpreoptsearchmat.clear();
    globalusrmat.clear();
   
    std::vector<string> schedulecoms(0);
    string comline;
    if (globalopts.schedulefname.length()<1) {
      setdefaultschedule(schedulecoms);
    } else {
      // open the schedule file
      ifstream schedulefile(globalopts.schedulefname.c_str());
      if (!schedulefile) {
	cerr << "Could not open file" << globalopts.schedulefname << endl;
	return -1;
      }
      while (!schedulefile.eof()) {
	getline(schedulefile,comline);
	schedulecoms.push_back(comline);
      }
      schedulefile.close();
    }

    // interpret each line in the schedule command vector
    bool skip=false;
    for (unsigned int i=0; i<schedulecoms.size(); i++) {
      comline = schedulecoms[i];
      if (globalopts.verbose>0) {
	cout << " >> " << comline << endl;
      }
      interpretcommand(comline,skip,testvol,refvol,refvol_2,refvol_4,refvol_8);
    }

    // re-read the initial volume, and transform it by the optimised result
    Matrix reshaped;
    if (globalusrmat[0].size()>0) {
      reshaped = (globalusrmat[0])[0].SubMatrix(1,1,2,17);
      reshape(matresult,reshaped,4,4);

      Matrix finalmat = matresult * globalopts.initmat;
      read_volume(testvol,globalopts.inputfname);
      read_volume(refvol,globalopts.reffname);
      save_matrix_data(finalmat,testvol,refvol);
      // generate the outputvolume (not safe_save st -out overrides -nosave)
      volume newtestvol = refvol;
      filled_affine_transform(newtestvol,testvol,finalmat);      
      save_volume(newtestvol,globalopts.outputfname.c_str());
      //save_global_data(finalmat,testvol,refvol);
    }

  }
  catch(Exception exc) {
    cerr << exc.what() << endl;
    throw;
  }
  catch(...) {
    cerr << "Image error" << endl;
    throw;
  } 
  return(0);
}

