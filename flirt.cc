/*  FLIRT - FMRIB's Linear Image Registration Tool

    flirt.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */

// Put current version number here:
#include <string>
const string version = "2.1.3";

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
#include "miscimfns.h"
#include "miscmaths.h"
#include "optimise.h"
#include "interpolation.h"
#include "mjimage.h"
#include "costfns.h"
#include "generalio.h"
#include "defaultschedule.h"
#include "globaloptions.h"

#ifndef NO_NAMESPACE
 using namespace MISCMATHS;
 using namespace MISCIMFNS;
 using namespace COSTFNS;
 using namespace INTERPOLATION;
 using namespace NEWMAT;
 using namespace MJIMAGE;
 using namespace GENERALIO;
#endif

////////////////////////////////////////////////////////////////////////////

void print_vector(float x, float y, float z)
{
  cerr << "(" << x << "," << y << "," << z << ")"; 
}


int round(const float val) {
  Tracer tr("round");
  if (val>0.0) {
    return (int) (val + 0.5);
  } else {
    return (int) (val - 0.5);
  }
}

//------------------------------------------------------------------------//
// Some interfaces to generalio

int safe_save_volume(const volume& source, const string& filename)
{
  if (!globaloptions::get().nosave) {
    save_volume(source,filename,globaloptions::get().datatype);
  }
  return 0;
}


void save_matrix_data(const Matrix& matresult, const volume& initvol, 
		      const volume& finalvol)
{
    write_ascii_matrix(matresult,globaloptions::get().outputmatascii);
    write_medx_matrix(matresult,globaloptions::get().outputmatmedx,
		      initvol,finalvol,"a",globaloptions::get().reffname);
}

void save_global_data(const Matrix& matresult, const volume& initvol,
		      const volume& finalvol)
{
    // save the data now
    safe_save_volume(globaloptions::get().impair->testvol,"inputvol");
    safe_save_volume(globaloptions::get().impair->refvol,"refvol");
    volume outputvol;
    outputvol = globaloptions::get().impair->refvol;
    filled_affine_transform(outputvol,globaloptions::get().impair->testvol,
		     matresult * globaloptions::get().initmat);
    safe_save_volume(outputvol,globaloptions::get().outputfname.c_str());
    save_matrix_data(matresult * globaloptions::get().initmat,initvol,finalvol);
}

void save_global_data(const Matrix& matresult) {
  save_global_data(matresult,globaloptions::get().impair->testvol,
		   globaloptions::get().impair->testvol);
}


float costfn(const Matrix& matresult);

void save_workspace_and_pause(const Matrix& matresult) {
  Tracer tr("save_workspace_and_pause");
  costfn(matresult);      
  save_global_data(matresult);
  if (globaloptions::get().interactive) {
    cerr << "Enter a number to continue..."; // forces the transform to be done
    int dummy; 
    cin >> dummy;
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

  switch (globaloptions::get().anglerep) 
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
  return vector2affine(params,n,globaloptions::get().impair->testvol.cog(),aff);
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
  switch (globaloptions::get().anglerep) 
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
  return affmat2vector(aff,n,globaloptions::get().impair->testvol.cog(),params);
}


void set_param_basis(Matrix &parambasis, int no_params)
{
  parambasis = 0.0;
  for (int i=1; i<=no_params; i++) {
    parambasis(i,i)=1.0;
  }
}


float estimate_scaling(const volume& vol) {
  Tracer tr("estimate_scaling");
  return Min(Min(vol.getx(),vol.gety()),vol.getz());
}

float estimate_scaling() {
  Tracer tr("estimate_scaling");
  return estimate_scaling(globaloptions::get().impair->refvol);
}

void set_param_tols(ColumnVector &param_tol, int no_params)
{
      // Tolerances are: 0.57 degrees (0.005 radians), 0.2mm translation
      //    0.005 scale and 0.001 skew
//    float diagonal[12]={0.005, 0.005, 0.005, 0.2, 0.2, 0.2, 0.002, 0.002, 0.002,
//    		      0.001, 0.001, 0.001};
  if (param_tol.Nrows()<no_params) {
    param_tol.ReSize(no_params);
  }
  for (int i=1; i<=no_params; i++) {
    //    param_tol(i)=diagonal[i-1];
    param_tol(i)=globaloptions::get().tolerance(i);
  }
  param_tol *= estimate_scaling();  // scale it up by the current scaling
}



void initialise_params(ColumnVector& params)
{
  Real paramsf[12] = {0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0};
  params << paramsf;
}


void optimise(ColumnVector& params, int no_params, ColumnVector& param_tol, 
	      int &no_its, float *fans, 
	      float (*costfunc)(const ColumnVector &), int itmax=4)
{
  // sets up the initial parameters and calls the optimisation routine
  if (params.MaximumAbsoluteValue() < 0.001)  initialise_params(params);
  { 
    Matrix affmattst(4,4);
    vector2affine(params,no_params,affmattst);
    if (globaloptions::get().verbose>=5) {
      cout << "Starting with : " << endl << affmattst;
    }
  }
  Matrix parambasis(no_params,no_params);
  set_param_basis(parambasis,no_params);
  float ptol[13];
  for (int i=1; i<=no_params; i++) { ptol[i] = param_tol(i); }
  
  *fans = MISCMATHS::optimise(params,no_params,param_tol,costfunc,no_its,itmax,
			      globaloptions::get().boundguess);
}


////////////////////////////////////////////////////////////////////////////

// OPTIMISATION SUPPORT (cost function interfaces)


float costfn(const Matrix& uninitaffmat)
{
  Tracer tr("costfn");
  Matrix affmat = uninitaffmat * globaloptions::get().initmat;  // apply initial matrix
  float retval = 0.0;
  switch (globaloptions::get().currentcostfn) 
    {
    case NormCorr:
      retval = 1.0 - fabs(normcorr(globaloptions::get().impair,affmat));  // MAXimise corr
      break;
    case CorrRatio:
      if (globaloptions::get().smoothsize > 0.0) {
	retval = 1.0 - corr_ratio_smoothed(globaloptions::get().impair,affmat);  // MAXimise corr
      } else {
	retval = 1.0 - corr_ratio(globaloptions::get().impair,affmat);  // MAXimise corr
      }
      break;
    case Woods:
      retval = woods_fn(globaloptions::get().impair,affmat);  // minimise variance/mean
      break;
    case MutualInfo:
      retval = -mutual_info(globaloptions::get().impair,affmat);  // MAXimise info
      break;
    case NormMI:
      retval = -normalised_mutual_info(globaloptions::get().impair,affmat);  // MAXimise
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
  vector2affine(params,globaloptions::get().no_params,affmat);
  float retval = costfn(affmat);
  if (globaloptions::get().verbose>=5) {
    cout << globaloptions::get().impair->count++ << " : ";
    cout << retval << " :: ";
    for (int i=1; i<=globaloptions::get().no_params; i++) 
      { cout << params(i) << " "; }
    cout << endl;
  }
  return retval;
}
  

//  float costfn(float params[])
//  {
//    Tracer tr("costfn");
//    Matrix affmat(4,4);
//    vector2affine(params,globaloptions::get().no_params,affmat);
//    float retval = costfn(affmat);
//    return retval;
//  }

//----------------------------------------------------------------------//


void params12toN(ColumnVector& params)
{
  // Convert the full 12 dof param vector to a small param vector
  Tracer tr("params12toN");
  ColumnVector nparams;
  nparams = pinv(globaloptions::get().parammask)*(params - globaloptions::get().refparams);
  params = nparams;
}


void paramsNto12(ColumnVector& params)
{
  // Convert small param vector to full 12 dof param vector
  Tracer tr("paramsNto12");
  ColumnVector param12;
  param12 = globaloptions::get().parammask*params + globaloptions::get().refparams;
  params = param12;
}



float subset_costfn(const ColumnVector& params)
{
  Tracer tr("subset_costfn");
  ColumnVector param12;
  param12 = params;
  paramsNto12(param12);
  float retval = costfn(param12);
  if (globaloptions::get().verbose>=7) {
    cout << globaloptions::get().impair->count++ << " : ";
    cout << retval << " :: " << param12.t() << endl;
  }
  return retval;
}


//  float subset_costfn(float params[])
//  {
//    Tracer tr("subset_costfn");
//    ColumnVector newparams(globaloptions::get().parammask.Ncols());
//    for (int n=1; n<=globaloptions::get().parammask.Ncols(); n++) {
//      newparams(n) = params[n];
//    }
//    return subset_costfn(newparams);
//  }



//------------------------------------------------------------------------//




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
	  if (cost(x,y,z) < 
	        cost(round(bestpt(1)),round(bestpt(2)),round(bestpt(3))))
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
	if ( cost.in_bounds(x,y,z) && cost.in_bounds(x+1,y+1,z+1) ) {
	  if (minv(x,y,z) < 0.5) {
	    bestpts(idx,1) = (float) x; 
	    bestpts(idx,2) = (float) y; 
	    bestpts(idx,3) = (float) z;
	    idx++;
	    if (globaloptions::get().verbose>=3)
	      cerr << "COST minima at : " << x << "," << y << "," << z << endl;
	  }
	}
      }
    }  
  }
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
  rxcoarse.ReSize(round((globaloptions::get().searchrx(2) 
			 - globaloptions::get().searchrx(1))
			/globaloptions::get().coarsedelta)+1);
  rycoarse.ReSize(round((globaloptions::get().searchry(2) 
			 - globaloptions::get().searchry(1))
			/globaloptions::get().coarsedelta)+1);
  rzcoarse.ReSize(round((globaloptions::get().searchrz(2) 
			 - globaloptions::get().searchrz(1))
			/globaloptions::get().coarsedelta)+1);
  rxfine.ReSize(round((globaloptions::get().searchrx(2) 
		       - globaloptions::get().searchrx(1))
			/globaloptions::get().finedelta)+1);
  ryfine.ReSize(round((globaloptions::get().searchry(2) 
		       - globaloptions::get().searchry(1))
			/globaloptions::get().finedelta)+1);
  rzfine.ReSize(round((globaloptions::get().searchrz(2) 
		       - globaloptions::get().searchrz(1))
			/globaloptions::get().finedelta)+1);
  // now get the appropriate angle sample values
  set_rot_sampling(rxcoarse,globaloptions::get().searchrx(1),
		   globaloptions::get().searchrx(2));
  set_rot_sampling(rycoarse,globaloptions::get().searchry(1),
		   globaloptions::get().searchry(2));
  set_rot_sampling(rzcoarse,globaloptions::get().searchrz(1),
		   globaloptions::get().searchrz(2));
  set_rot_sampling(rxfine,globaloptions::get().searchrx(1),
		   globaloptions::get().searchrx(2));
  set_rot_sampling(ryfine,globaloptions::get().searchry(1),
		   globaloptions::get().searchry(2));
  set_rot_sampling(rzfine,globaloptions::get().searchrz(1),
		   globaloptions::get().searchrz(2));
  if (globaloptions::get().verbose>=4) {
    cout << "Coarse rotation samplings are:\n" << rxcoarse.t() << rycoarse.t() 
	 << rzcoarse.t() << " and fine rotation samplings are:\n"
	 << rxfine.t() << ryfine.t() << rzfine.t() << endl;
  }
}



void search_cost(Matrix& paramlist, volume& costs, volume& tx, 
		 volume& ty, volume& tz, volume& scale) {
  Tracer tr("search_cost");
  int storedverbose = globaloptions::get().verbose;
  int storeddof = globaloptions::get().dof;
  anglereps useranglerep = globaloptions::get().anglerep;
  globaloptions::get().verbose -= 2;
  globaloptions::get().anglerep = Euler;  // a workaround hack
  globaloptions::get().currentcostfn = globaloptions::get().searchcostfn;
  globaloptions::get().dof = globaloptions::get().searchdof;

  ColumnVector coarserx, coarsery, coarserz, finerx, finery, finerz;
  set_rot_samplings(coarserx,coarsery,coarserz,finerx,finery,finerz);

  // set up the type of parameter subset (via the global mask)
  // here 3 translations and 1 (common) scaling are used
  if (globaloptions::get().dof>6) {
    globaloptions::get().parammask.ReSize(12,4);  // was 3
    globaloptions::get().parammask = 0.0;
    globaloptions::get().parammask(7,1) = 1.0;  // didn't used to exist
    globaloptions::get().parammask(8,1) = 1.0;  // didn't used to exist
    globaloptions::get().parammask(9,1) = 1.0;  // didn't used to exist
    globaloptions::get().parammask(4,2) = 1.0;
    globaloptions::get().parammask(5,3) = 1.0;
    globaloptions::get().parammask(6,4) = 1.0;
  } else {
    globaloptions::get().parammask.ReSize(12,3);
    globaloptions::get().parammask = 0.0;
    globaloptions::get().parammask(4,1) = 1.0;
    globaloptions::get().parammask(5,2) = 1.0;
    globaloptions::get().parammask(6,3) = 1.0;
  }
	
  ColumnVector param_tol, param_tol0(12), param_tol1(12), params_8(12);
  globaloptions::get().no_params = 12; // necessary for any subset_costfn call
  param_tol0 = globaloptions::get().refparams;
  params12toN(param_tol0);
  set_param_tols(param_tol1,12);
  param_tol1 = param_tol1 + globaloptions::get().refparams;
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
  globaloptions::get().refparams = 0.0;
  globaloptions::get().refparams(7) = 1.0;
  globaloptions::get().refparams(8) = 1.0;
  globaloptions::get().refparams(9) = 1.0;
  // set the initial translation (to align cog's)
  trans = globaloptions::get().impair->refvol.cog() 
            - globaloptions::get().impair->testvol.cog();
  globaloptions::get().refparams(4) = trans(1);
  globaloptions::get().refparams(5) = trans(2);
  globaloptions::get().refparams(6) = trans(3);
  for (int ix=0; ix<coarserx.Nrows(); ix++) {
    for (int iy=0; iy<coarsery.Nrows(); iy++) {
      for (int iz=0; iz<coarserz.Nrows(); iz++) {
	rx = coarserx(ix+1);
	ry = coarsery(iy+1);
	rz = coarserz(iz+1);
	globaloptions::get().refparams(1) = rx;
	globaloptions::get().refparams(2) = ry;
	globaloptions::get().refparams(3) = rz;
	params_8 = globaloptions::get().refparams;
	if (globaloptions::get().verbose>=4) {
	  cout << "Starting with " << params_8.t();
	  cout << "  and tolerance " << param_tol.t();
	}
	params12toN(params_8);
	optimise(params_8,globaloptions::get().parammask.Ncols(),
		 param_tol,no_its,&fans,subset_costfn);
	//		 param_tol,&no_its,&fans,subset_costfn);
	paramsNto12(params_8);
	tx(ix,iy,iz) = params_8(4);
	ty(ix,iy,iz) = params_8(5);
	tz(ix,iy,iz) = params_8(6);
	scale(ix,iy,iz) = params_8(7);
	
	if (globaloptions::get().verbose>=4) {
	  cout << " dearranged: " << params_8.t();
	}
      }
      if (globaloptions::get().verbose>=2) cerr << "*";
    }
  }
  if (globaloptions::get().verbose>=2) cerr << endl;

  // scale = 1.0;  // for now disallow non-unity scalings
  float medianscale = get_percentile(scale,50.0);
  scale = medianscale;  // try a constant, median scale for all
  if (globaloptions::get().verbose>=2) 
    { cerr << "Median scale = " << medianscale << endl; }

  if (globaloptions::get().verbose>=4) {
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
	globaloptions::get().refparams(1) = rx;
	globaloptions::get().refparams(2) = ry;
	globaloptions::get().refparams(3) = rz;
	xf = ((float) ix)*factorx;
	yf = ((float) iy)*factory;
	zf = ((float) iz)*factorz;
	txv = tri_interpolation(tx,xf,yf,zf);
	tyv = tri_interpolation(ty,xf,yf,zf);
	tzv = tri_interpolation(tz,xf,yf,zf);
	scv = tri_interpolation(scale,xf,yf,zf);
	if ((scv<0.5) || (scv>2.0))  scv = 1.0;
	if (globaloptions::get().dof<=6) scv = 1.0;
	globaloptions::get().refparams(4) = txv;
	globaloptions::get().refparams(5) = tyv;
	globaloptions::get().refparams(6) = tzv;
	globaloptions::get().refparams(7) = scv;
	globaloptions::get().refparams(8) = scv;
	globaloptions::get().refparams(9) = scv;
	params_8 = globaloptions::get().refparams;
	costs(ix,iy,iz) = costfn(params_8);
      }
      if (globaloptions::get().verbose>=2) cerr << "*";
    }
  }
  if (globaloptions::get().verbose>=2) cerr << endl;

  if (globaloptions::get().verbose>=4) {
    safe_save_volume(costs,"costs");
  }


  float costmin, costmax;
  get_min_max(costs,costmin,costmax);
  // the following assumes that costmin > 0
  float factor = 0.2;
  float costthresh = Min(costmin + factor*(costmax-costmin),
			 get_percentile(costs,20));
  // avoid the percentile giving the costmin (or less)
  if (costthresh <= costmin)  costthresh = Max(costmin*1.0001,costmin*0.9999);
  if (globaloptions::get().verbose>=4) {
    cout << "Cost threshold = " << costthresh << " and minimum is "
	 << costmin << endl;
  }

  // for all costs less than 150% of the best cost so far, optimise...
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
	  if (globaloptions::get().dof<=6) scv = 1.0;
	  params_8 = 0.0;
	  params_8(1) = rx;  params_8(2) = ry;  params_8(3) = rz;
	  params_8(4) = txv; params_8(5) = tyv; params_8(6) = tzv;
	  params_8(7) = scv; params_8(8) = scv; params_8(9) = scv;
	  globaloptions::get().refparams = params_8;
	  params12toN(params_8);
	  optimise(params_8,globaloptions::get().parammask.Ncols(),
		   param_tol,no_its,&fans,subset_costfn);
	           //param_tol,&no_its,&fans,subset_costfn);
	  paramsNto12(params_8);
	  costs(ix,iy,iz) = fans;
	  bestparams(n,1) = fans;
	  bestparams.SubMatrix(n,n,2,13) = params_8.t();
	  n++;
	  if (globaloptions::get().verbose>=3) {
	    cout << "(" << ix << "," << iy << "," << iz << ") => " << fans
		 << " with " << params_8.t();
	  }
	}
      }
    }
  }

  if (globaloptions::get().verbose>=3) {
    safe_save_volume(costs,"costs");
    cout << "Costs (1st column) are:\n" << bestparams << endl;
  }

  // find the cost minima and return these
  Matrix bestpts;
  find_cost_minima(bestpts,costs);
  paramlist.ReSize(bestpts.Nrows(),12);
  int ix,iy,iz;
  for (int n=1; n<=paramlist.Nrows(); n++) {
    ix = round(bestpts(n,1));
    iy = round(bestpts(n,2));
    iz = round(bestpts(n,3));
    if (globaloptions::get().verbose>=3) 
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
    if (globaloptions::get().dof<=6) scv = 1.0;
    params_8 = 0.0;
    params_8(1) = rx;  params_8(2) = ry;  params_8(3) = rz;
    params_8(4) = txv; params_8(5) = tyv; params_8(6) = tzv;
    params_8(7) = scv; params_8(8) = scv; params_8(9) = scv;
    paramlist.SubMatrix(n,n,1,12) = params_8.t();
  }

  if (globaloptions::get().verbose>=3) {
    cout << "Chosen parameters:\n" << paramlist << endl;
  }

  globaloptions::get().anglerep = useranglerep;
  globaloptions::get().verbose = storedverbose;
  globaloptions::get().dof = storeddof;
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
  globaloptions::get().no_params = dof;
  set_param_tols(param_tol,12);  // 12 used to be dof
  affmat2vector(matresult,dof,params);
  //optimise(params,dof,param_tol,&no_its,&fans,costfn,max_iterations);
  optimise(params,dof,param_tol,no_its,&fans,costfn,max_iterations);
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

  dof = Min(globaloptions::get().dof,7);
  no_its = optimise_strategy1(matresult,costval,dof);

  if ((no_its<3) || (globaloptions::get().dof>dof)) {
    dof = Min(globaloptions::get().dof,9);
    optimise_strategy1(matresult,costval,dof,1);
  }

  if (globaloptions::get().dof>=10) {
    dof = Min(globaloptions::get().dof,12);
    optimise_strategy1(matresult,costval,dof,2);
  }
}

//-------------------------------------------------------------------------//


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
	 << " rows to matrix\n";
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
  globaloptions::get().currentcostfn = globaloptions::get().searchcostfn;
  search_cost(paramlist,costs,tx,ty,tz,scale);
  
  int dof = Min(globaloptions::get().dof,7);
  ColumnVector params_8(12);
  float costval=0.0, rms_min = estimate_scaling();
  Matrix reshapedmat(1,16), matresult(4,4), premat(4,4), matrow(1,34);
  opt_matrixlist.ReSize(0,34);
  // freely optimise each member of the parameter list (allow rotns to vary)
  int verbose = globaloptions::get().verbose;
  globaloptions::get().verbose-=2;
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
      cout << costval << " ::: " << params_8.t(); 
    }
  }
  globaloptions::get().verbose = verbose;
  globaloptions::get().currentcostfn = globaloptions::get().maincostfn;

  if (globaloptions::get().verbose>=3) {
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
  // set the magnitude of the variations
  float delscale, delrx, delry, delrz;
  delrx = 3.0*param_tol(1); 
  delry = 3.0*param_tol(2);
  delrz = 3.0*param_tol(3);
  if (finerx.Nrows()>1) { delrx = 0.5*(finerx(2) - finerx(1)); }
  if (finery.Nrows()>1) { delry = 0.5*(finery(2) - finery(1)); }
  if (finerz.Nrows()>1) { delrz = 0.5*(finerz(2) - finerz(1)); }
  delscale=0.0;  
  if (globaloptions::get().dof>=7) delscale = 0.1;  // scale only set if allowed
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
  int dof = Min(globaloptions::get().dof,7);
  float bestcost=0.0, costval=0.0;
  Matrix bestmat(4,4), reshaped(1,16);
  ColumnVector params_8(12);
  
  Matrix delta, perturbmask;
  set_perturbations(delta,perturbmask);

  // for each potentially best transform do an optimisation
  int verbose = globaloptions::get().verbose;
  globaloptions::get().verbose -= 2;
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

	if (globaloptions::get().verbose>=3) {
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

  globaloptions::get().verbose = verbose;
  matresult = bestmat;  // the return value
}


////////////////////////////////////////////////////////////////////////////


int get_testvol(volume& testvol)
{
  Tracer tr("get_testvol");
  short dtype;
  read_volume(testvol,globaloptions::get().inputfname,dtype);
  if (!globaloptions::get().forcedatatype)
    globaloptions::get().datatype = dtype;
  read_matrix(globaloptions::get().initmat,globaloptions::get().initmatfname,testvol);

  ColumnVector hist;
  float minval=0.0, maxval=0.0;
  find_robust_limits(testvol,globaloptions::get().no_bins,hist,minval,maxval);
  clamp(testvol,minval,maxval);

  if (globaloptions::get().verbose>=2) {
    cerr << "Init Matrix = \n" << globaloptions::get().initmat << endl;
    cerr << "Testvol sampling matrix =\n" << testvol.sampling_matrix() << endl;
    cerr << "Testvol Data Type = " << dtype << endl;
    cerr << "Testvol intensity clamped between " 
	 << minval << " and " << maxval << endl;
  }
  return 0;
}  


int get_refvol(volume& refvol)
{
  Tracer tr("get_refvol");
  read_volume(refvol,globaloptions::get().reffname);

  ColumnVector hist;
  float minval=0.0, maxval=0.0;
  find_robust_limits(refvol,globaloptions::get().no_bins,hist,minval,maxval);
  clamp(refvol,minval,maxval);

  if (globaloptions::get().verbose>=2) {
    cerr << "Refvol intensity clamped between " 
	 << minval << " and " << maxval << endl;
  }
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

  if (globaloptions::get().verbose >= 2) 
    cout << "Resampling refvol isotropically: " << sampling << endl;

  // isotropically resample the volume
  volume tmpvol;
  blur4subsampling(tmpvol,refvol,sampl);
  isotropic_rescale(refvol,tmpvol,sampl);

  if (globaloptions::get().verbose>=2) print_volume_info(refvol,"Refvol");      
  
  return 0;
}


////////////////////////////////////////////////////////////////////////////


void no_optimise()
{
  Tracer tr("no_optimise");
  volume refvol, testvol;
  // set up image pair and global pointer
  
  get_refvol(refvol);
  short dtype;
  read_volume(testvol,globaloptions::get().inputfname,dtype);
  if (!globaloptions::get().forcedatatype)
    globaloptions::get().datatype = dtype;
  read_matrix(globaloptions::get().initmat,
	      globaloptions::get().initmatfname,testvol);
  if (globaloptions::get().verbose>=2) {
    cerr << "Init Matrix = \n" << globaloptions::get().initmat << endl;
  }
  
  float min_sampling_ref=1.0;
  min_sampling_ref = Min(refvol.getx(),Min(refvol.gety(),refvol.getz()));
  
  if (globaloptions::get().iso) {
    resample_refvol(refvol,globaloptions::get().isoscale);
  }
  {
    volume testvol_1;
    blur4subsampling(testvol_1,testvol,min_sampling_ref);
    testvol = testvol_1;
  }

  if (globaloptions::get().verbose>=2) { 
    print_volume_info(refvol,"refvol"); 
    print_volume_info(testvol,"inputvol"); 
  }

  volume outputvol = refvol;
  filled_affine_transform(outputvol,testvol,globaloptions::get().initmat);
  safe_save_volume(outputvol,globaloptions::get().outputfname.c_str());
  if (globaloptions::get().verbose>=2) {
    save_matrix_data(globaloptions::get().initmat, testvol, outputvol);
    print_volume_info(outputvol,"Resampled volume");
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
    line = line.substr(i2,string::npos);
  }
}


void striptrailingspace(string& line)
{
  Tracer tr("striptrailingspace");
  if (line.length()<1) return;
  string::size_type i2=line.find_last_not_of("	 ");
  if (i2==string::npos) {
    line="";
    return;
  }
  if (i2!=0) {
    line = line.substr(0,i2+1);
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
    striptrailingspace(name);
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
    namedmat = &globaloptions::get().searchoptmat;
    return 0;
  } else if (name == "P") {
    namedmat = &globaloptions::get().preoptsearchmat;
    return 0;
  } else if (name == "U") {
    namedmat = &(globaloptions::get().usrmat[0]);
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
    namedmat = &(globaloptions::get().usrmat[idx]);
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
    namedscalar = (float) globaloptions::get().dof;
    return 0;
  } else if (name == "MINSAMPLING") {
    namedscalar = (float) globaloptions::get().min_sampling;
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
  namedscalar = round(temp);
  return retval;
}


int parsematname(const string& inname, MatVecPtr& usrdefmat, int& r1, int& r2) 
{
  Tracer tr("parsematname");
  string name = inname;
  stripleadingspace(name);
  striptrailingspace(name);
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

void usrcopy(MatVecPtr usrsrcmat, MatVecPtr usrdestmat, 
	     unsigned int usrsrcrow1, unsigned int usrsrcrow2) 
{
  Tracer tr("usrcopy");
  // COPY src -> dest
  for (unsigned int crow=usrsrcrow1; crow<=Min(usrsrcrow2,usrsrcmat->size()); 
       crow++)
    {
      usrdestmat->push_back((*usrsrcmat)[crow-1]);
    }
}


void usrclear(MatVecPtr usrsrcmat) 
{
  Tracer tr("usrclear");
  // CLEAR src
  if (usrsrcmat->size() > 0) {
    usrsrcmat->erase(usrsrcmat->begin(),usrsrcmat->end());
  }
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



int usrsetoption(const std::vector<string> &words)
{
  Tracer tr("usrsetoption");
  // SETOPTION name values

  int len = words.size();
  if (len<3) {
    cerr << "Must specify an option and a value" << endl;
    return -2;
  }
    
  string option = words[1];

  ColumnVector fvalues(len-2);
  float fval=0.0;

  for (int i=2; i<len; i++) {
    setscalarvariable(words[i],fval);
    fvalues(i-1) = fval;
  }

  if (globaloptions::get().verbose>2) {
    cout << "setoptions: " << option << " " << fvalues.t() << endl;
  }

  if (option=="smoothing") {
    globaloptions::get().smoothsize = fvalues(1);
    return 0;
  } else if (option=="tolerance") {
    globaloptions::get().tolerance = fvalues/estimate_scaling();
    // Note: division by estimate_scaling used so that the absolute
    //       tolerance used at this scale is that specified by the user
    return 0;
  } else if (option=="boundguess") {
    globaloptions::get().boundguess = fvalues;
    return 0;
  } else if (option=="minsampling") {
    globaloptions::get().min_sampling = fvalues(1);
    return 0;
  } else {
    cerr << "Option " << option << " is unrecognised - ignoring" << endl;
  }
  return -1;
}


void usrsort(MatVecPtr usrsrcmat) 
{
  Tracer tr("usrsort");
  // SORT src
  sort(usrsrcmat->begin(),usrsrcmat->end(),firstelementless);
}


void usrdualsort(MatVecPtr usrsrcmat1, MatVecPtr usrsrcmat2) 
{
  Tracer tr("usrdualsort");
  // DUALSORT
  if (usrsrcmat1->size() != usrsrcmat2->size()) {
    cerr << "Cannot dual sort matrices of unequal size";
    return;
  }
  RowVector tmprow(34);
  MatVec optmatsorted(0);
  for (unsigned int i=0; i<usrsrcmat1->size(); i++) {
    tmprow.SubMatrix(1,1,1,17) = (*usrsrcmat1)[i];
    tmprow.SubMatrix(1,1,18,34) = (*usrsrcmat2)[i];
    optmatsorted.push_back(tmprow);
  }
  sort(optmatsorted.begin(),optmatsorted.end(),firstelementless);
  usrclear(usrsrcmat1);
  usrclear(usrsrcmat2);
  tmprow.ReSize(17);
  for (unsigned int i=0; i<optmatsorted.size(); i++) {
    tmprow = optmatsorted[i].SubMatrix(1,1,1,17);
    usrsrcmat1->push_back(tmprow);
    tmprow = optmatsorted[i].SubMatrix(1,1,18,34);
    usrsrcmat2->push_back(tmprow);
  }
}


void usrsearch(MatVecPtr searchoptmat, MatVecPtr preoptsearchmat, int sdof=12) 
{
  Tracer tr("usrsearch");
  // SEARCH
  Matrix opt_matrixlist;
  globaloptions::get().searchdof = sdof;
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
  int dof = Min(globaloptions::get().dof,usrdof);
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
  int dof = Min(globaloptions::get().dof,usrdof);
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
  //  the object pointed to in globaloptions::get().impair does not go 
  //  out of scope
  Tracer tr("usrsetscale");
  int scale = usrscale;
  imagepair *globalpair=0;
  if (globaloptions::get().min_sampling<=1.25 * ((float) scale)) {
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
    globalpair->set_no_bins(globaloptions::get().no_bins/scale);
    globalpair->smoothsize = globaloptions::get().smoothsize;
    if (globaloptions::get().verbose>=3) {
      if (globaloptions::get().impair) {
	cout << "Previous scale used " << globaloptions::get().impair->count
	     << " cost function evaluations" << endl;
      }
    }
    globaloptions::get().impair = globalpair;
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
    usrcopy(src,dest,d1,d2);
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
  } else if (words[0]=="exit") {
    // EXIT
    exit(0);
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
  } else if (words[0]=="dualsort") {
    // DUALSORT
    if (words.size()<3) {
      cerr << "Wrong number of args to DUALSORT" << endl;
      exit(-1);
    }
    MatVecPtr mat1, mat2;
    int d1, d2, d3, d4;
    parsematname(words[1],mat1,d1,d2);
    parsematname(words[2],mat2,d3,d4);
    usrdualsort(mat1,mat2);
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
    usrsearch(&globaloptions::get().searchoptmat,
	      &globaloptions::get().preoptsearchmat,usrdof);
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
    usroptimise(&(globaloptions::get().usrmat[0]),
		usrdefmat,usrrow1,usrrow2,usrdof, 
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
    usrmeasurecost(&(globaloptions::get().usrmat[0]),
		   usrdefmat,usrrow1,usrrow2,usrdof, 
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
  } else if (words[0]=="setoption") {
    // SETOPTION
    if (words.size()<3) {
      cerr << "Wrong number of args to SETOPTION" << endl;
      exit(-1);
    }
    usrsetoption(words);
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
    cerr << " ... ignoring" << endl;
    //cerr << "Quitting." << endl;
    //exit(-1);
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

  //test_proc();  // only used for some debugging

  globaloptions::get().parse_command_line(argc, argv,version);

  if (!globaloptions::get().do_optimise) {   
    no_optimise();
  }
    
  // read in the volumes
  volume refvol, testvol;
  get_refvol(refvol);
  get_testvol(testvol);

  float min_sampling_ref=1.0, min_sampling_test=1.0, min_sampling=1.0;
  min_sampling_ref = Min(refvol.getx(),Min(refvol.gety(),refvol.getz()));
  min_sampling_test = Min(testvol.getx(),Min(testvol.gety(),testvol.getz()));
  min_sampling = (float) ceil(Max(min_sampling_ref,min_sampling_test));
  if (globaloptions::get().min_sampling < min_sampling) 
    globaloptions::get().min_sampling = min_sampling;

    
  if (globaloptions::get().verbose>=3) {
    cout << "CoG for refvol is:  " << refvol.cog().t();
    cout << "CoG for testvol is:  " << testvol.cog().t();
  }
  
  volume refvol_2, refvol_4, refvol_8;
  if (globaloptions::get().resample) {
    // set up subsampled volumes by factors of 2, 4 and 8
    if (globaloptions::get().verbose >= 2) 
      cout << "Subsampling the volumes" << endl;
      
    resample_refvol(refvol,globaloptions::get().min_sampling);
    if (globaloptions::get().min_sampling < 1.9) {
      subsample_by2(refvol,refvol_2);
    } else {
      refvol_2 = refvol;
    }
    if (globaloptions::get().min_sampling < 3.9) {
      subsample_by2(refvol_2,refvol_4);
    } else {
      refvol_4 = refvol_2;
    }
    if (globaloptions::get().min_sampling < 7.9) {
      subsample_by2(refvol_4,refvol_8);
    } else {
      refvol_8 = refvol_4;
    }
      
    {
      volume testvol_8;
      blur4subsampling(testvol_8,testvol,8.0);
      testvol = testvol_8;
    }
  } else {
    // if no resampling chosen, then refvol is simply copied
    refvol_8 = refvol;
    refvol_4 = refvol;
    refvol_2 = refvol;
  }

  // set up image pair and global pointer
  imagepair global_8(refvol_8,testvol);
  global_8.set_no_bins(globaloptions::get().no_bins/8);
  global_8.smoothsize = globaloptions::get().smoothsize;
  if (globaloptions::get().verbose>=2) print_volume_info(testvol,"TESTVOL");

  globaloptions::get().impair = &global_8;


  Matrix matresult(4,4);
  ColumnVector params_8(12), param_tol(12);

    // perform the optimisation

  globaloptions::get().currentcostfn = globaloptions::get().maincostfn;
   
  std::vector<string> schedulecoms(0);
  string comline;
  if (globaloptions::get().schedulefname.length()<1) {
    setdefaultschedule(schedulecoms);
  } else {
    // open the schedule file
    ifstream schedulefile(globaloptions::get().schedulefname.c_str());
    if (!schedulefile) {
      cerr << "Could not open file" << globaloptions::get().schedulefname << endl;
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
    if (globaloptions::get().verbose>=1) {
      cout << " >> " << comline << endl;
    }
    interpretcommand(comline,skip,testvol,refvol,refvol_2,refvol_4,refvol_8);
  }

  // re-read the initial volume, and transform it by the optimised result
  Matrix reshaped;
  if (globaloptions::get().usrmat[0].size()>0) {
    reshaped = (globaloptions::get().usrmat[0])[0].SubMatrix(1,1,2,17);
    reshape(matresult,reshaped,4,4);

    Matrix finalmat = matresult * globaloptions::get().initmat;
    read_volume(testvol,globaloptions::get().inputfname);
    read_volume(refvol,globaloptions::get().reffname);
    save_matrix_data(finalmat,testvol,refvol);
    // generate the outputvolume (not safe_save st -out overrides -nosave)
    if (globaloptions::get().outputfname.size()>0) {
      volume newtestvol = refvol;
      filled_affine_transform(newtestvol,testvol,finalmat);      
      save_volume(newtestvol,globaloptions::get().outputfname.c_str(),
		  globaloptions::get().datatype);
    }
    if ( (globaloptions::get().outputmatascii.size()<=0) && 
	 (globaloptions::get().outputmatmedx.size()<=0) ) {
      cout << endl << "Final result: " << endl << finalmat << endl;
    }
  }

  return(0);
}

