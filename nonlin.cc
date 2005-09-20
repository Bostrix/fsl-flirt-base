/*  nonlin.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2005 University of Oxford  */

/*  CCOPYRIGHT  */

// Initial effort at generic non-linear registration

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

#include "newimage/newimageall.h"
#include "newimage/costfns.h"
#include "newimage/warpfns.h"
#include "miscmaths/miscmaths.h"
#include "utils/options.h"

using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace Utilities;

// The two strings below specify the title and example usage that is
//  printed out as the help or usage message

string title="nonlin (Version 1.0)\nCopyright(c) 2005, University of Oxford (Mark Jenkinson)";
string examples="nonlin [options] -i <input image> -r <reference image> -o <output image>";

// Each (global) object below specificies as option and can be accessed
//  anywhere in this file (since they are global).  The order of the
//  arguments needed is: name(s) of option, default value, help message,
//       whether it is compulsory, whether it requires arguments
// Note that they must also be included in the main() function or they
//  will not be active.

Option<bool> verbose(string("-v,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);
Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);
Option<bool> debug(string("-d,--debug"), false, 
		     string("switch on debugging output"), 
		     false, no_argument);
Option<bool> nullbc(string("--nullbc"), false, 
		     string("set null boundary condition (zero gradient)"), 
		     false, no_argument);
Option<int> maxiter(string("--maxiter"), 2, 
		     string("maximum number of iterations"), 
		     false, requires_argument);
Option<int> nbins(string("--nbins"), 80, 
		  string("number of bins (default=80)"), 
		  false, requires_argument);
Option<float> blursize(string("--blursize"), 20, 
		     string("amount of blurring of gradient (in mm)"), 
		     false, requires_argument);
Option<string> inname(string("-i,--in"), string(""),
		      string("input image filename"),
		      true, requires_argument);
Option<string> refname(string("-r,--ref"), string(""),
		       string("reference image filename"),
		       true, requires_argument);
Option<string> outname(string("-o,--out"), string(""),
		       string("output image filename"),
		       true, requires_argument);
Option<string> warpname(string("-w,--warp"), string(""),
			string("output warp-field filename"),
			false, requires_argument);
Option<string> affname(string("--initaff"), string(""),
			string("filename for initial affine matrix"),
			false, requires_argument);
Option<string> initwarp(string("--initwarp"), string(""),
			string("filename for initial warp"),
			false, requires_argument);
int nonoptarg;

Costfn* globalcostfnptr=0;

////////////////////////////////////////////////////////////////////////////

// Optimisation functions

  // A temporary fix of including the std:: in front of all abs() etc
  //  has been done for now
  using std::abs;

  bool estquadmin(float &xnew, float x1, float xmid, float x2, 
		   float y1, float ymid, float y2)
  {
    // Finds the estimated quadratic minimum's position
    float ad=0.0, bd=0.0, det=0.0;
    ad = (xmid - x2)*(ymid - y1) - (xmid - x1)*(ymid - y2);
    bd = -(xmid*xmid - x2*x2)*(ymid - y1) + (xmid*xmid - x1*x1)*(ymid - y2);
    det = (xmid - x2)*(x2 -x1)*(x1 - xmid);
    if ((fabs(det)>1e-15) && (ad/det < 0)) {  // quadratic only has a maxima
      xnew = 0.0;
      return false;
    }
    if (fabs(ad)>1e-15) {
      xnew = -bd/(2*ad);
      return true;
    } else {  // near linear condition -> get closer to an end point
      xnew = 0.0;
      return false;
    }
    return false;
  }


  float extrapolatept(float x1, float xmid, float x2)
  {
    // xmid must be between x1 and x2
    // use the golden ratio (scale similar result)
    const float extensionratio = 0.3819660;
    float xnew;
    if (fabs(x2-xmid)>fabs(x1-xmid)) {
      xnew = extensionratio * x2 + (1 - extensionratio) * xmid;
    } else {
      xnew = extensionratio * x1 + (1 - extensionratio) * xmid;
    }
    return xnew;
  }
  


  float nextpt(float x1, float xmid, float x2, float y1, float ymid, float y2)
  {
    // x1 and x2 are the bounds, xmid is between them

    float xnew;
    bool quadok=false;
    quadok = estquadmin(xnew,x1,xmid,x2,y1,ymid,y2);

    // check to see that the quadratic result is in the range
    if ((!quadok) || (xnew < Min(x1,x2)) || (xnew > Max(x1,x2))) {
      xnew = extrapolatept(x1,xmid,x2);
    }
    return xnew;
  }

      

  void findinitialbound(float &x1, float &xmid, float &x2, 
			float &y1, float &ymid, float &y2, 
			float (*func)(const volume4D<float> &),
			const volume4D<float> &unitdir, 
			const volume4D<float> &pt)
  {
    const float extrapolationfactor = 1.6;
    const float maxextrap = extrapolationfactor*2;
    if (y1==0)  y1 = (*func)(x1*unitdir + pt);
    if (ymid==0)  ymid = (*func)(xmid*unitdir + pt);
    if (y1<ymid) {   // swap a and b if this is the case
      float tempx = x1, tempy = y1;
      x1 = xmid;     y1 = ymid;
      xmid = tempx;  ymid = tempy;
    }

    float newx2 = 0.0, newy2=0.0, maxx2=0.0;
    float dir=1.0;
    if (xmid<x1) dir=-1.0;

    bool quadok;

    x2 = xmid + extrapolationfactor*(xmid - x1);
    y2 = (*func)(x2*unitdir + pt);

    while (ymid > y2) {  // note: must maintain y1 >= ymid
	
      // cout << "    <" << Min(x1,x2) << "," << xmid 
      //   << "," << Max(x1,x2) << ">" << endl;
      maxx2 = xmid + maxextrap*(x2 - xmid);
      quadok = estquadmin(newx2,x1,xmid,x2,y1,ymid,y2);
      if ((!quadok) || ((newx2 - x1)*dir<0) || ((newx2 - maxx2)*dir>0)) {
	newx2 = xmid + extrapolationfactor*(x2-x1);
      }
      
      newy2 = (*func)(newx2*unitdir + pt);

      if ((newx2 - xmid)*(newx2 - x1)<0) {  // newx2 is between x1 and xmid
	if (newy2 < ymid) {  // found a bracket!
	  x2 = xmid;  y2 = ymid;
	  xmid = newx2;  ymid = newy2;
	  break;
	} else {  // can use newx2 as a new value for x1 (as newy2 >= ymid)
	  x1 = newx2;  y1 = newy2;
	}
      } else {  // newx2 is between xmid and maxx2
	if (newy2 > ymid) { // found a bracket!
	  x2 = newx2;  y2 = newy2;
	  break;
	} else if ((newx2 - x2)*dir<0) {  // newx2 closer to xmid than old x2
	  x1 = xmid;  y1 = ymid;
	  xmid = newx2;  ymid = newy2;
	} else {
	  x1 = xmid;  y1 = ymid;
	  xmid = x2;  ymid = y2;
	  x2 = newx2;  y2 = newy2;
	}
      }
	
    }

    if ( (y2<ymid) || (y1<ymid) ) {
      cerr << "findinitialbound failed to bracket: current triplet is" << endl;
    }
  }
  

  float optimise1d(volume4D<float> &pt, const volume4D<float>& unitdir, 
		  float unittol, int &iterations_done, 
		  float (*func)(const volume4D<float>&), int max_iter,
		  float init_value, float boundguess) 
  {
    // Golden Search Routine
    // Must pass in the direction vector in N-space (dir), the initial
    //  N-dim point (pt), the acceptable tolerance (tol) and other
    //  stuff
    // Note that the length of the direction vector is unimportant
    // Pass in previous costfn value as init_value, if known, otherwise
    //  pass in 0.0 and it will force the calculation
    // Unlike the version in optimise.cc the boundguess is in absolute
    //  units, not in units of unittol

    float y1,y2,ymid;
    float x1,x2,xmid;

    // set up initial points
    xmid = 0.0;
    x1 = boundguess;  // initial guess (bound)
    if (init_value==0.0) ymid = (*func)(xmid*unitdir + pt);
    else ymid = init_value;
    y1 = (*func)(x1*unitdir + pt);
    findinitialbound(x1,xmid,x2,y1,ymid,y2,func,unitdir,pt);

    // cout << "(" << x1 << "," << y1 << ")  ";
    // cout << "(" << xmid << "," << ymid << ")  ";
    // cout << "(" << x2 << "," << y2 << ")" << endl;

    float min_dist = 0.1 * unittol;
    float xnew, ynew;
    int it=0;
    while ( ((++it)<=max_iter) && (fabs((x2-x1)/unittol)>1.0) )
      {
	// cout << "  [" << Min(x1,x2) << "," << Max(x1,x2) << "]" << endl;

	if (it>0) {
	  xnew = nextpt(x1,xmid,x2,y1,ymid,y2);
	} else {
	  xnew = extrapolatept(x1,xmid,x2);
	}

	float dirn=1.0;
	if (x2<x1) dirn=-1.0;

	if (fabs(xnew - x1)<min_dist) {
	  xnew = x1 + dirn*min_dist;
	}

	if (fabs(xnew - x2)<min_dist) {
	  xnew = x2 - dirn*min_dist;
	}

	if (fabs(xnew - xmid)<min_dist) {
	  xnew = extrapolatept(x1,xmid,x2);
	}

	if (fabs(xmid - x1)<0.4*unittol) {
	  xnew = xmid + dirn*0.5*unittol;
	}

	if (fabs(xmid - x2)<0.4*unittol) {
	  xnew = xmid - dirn*0.5*unittol;
	}

	if (verbose.value()) { cout << "xnew = " << xnew << endl; }
	ynew = (*func)(xnew*unitdir + pt);

	if ((xnew - xmid)*(x2 - xmid) > 0) {  // is xnew between x2 and xmid ?
	  // swap x1 and x2 so that xnew is between x1 and xmid
	  float xtemp = x1;  x1 = x2;  x2 = xtemp;
	  float ytemp = y1;  y1 = y2;  y2 = ytemp;
	}
	if (ynew < ymid) {
	  // new interval is [xmid,x1] with xnew as best point in the middle
	  x2 = xmid;  y2 = ymid;
	  xmid = xnew;  ymid = ynew;
	} else {
	  // new interval is  [x2,xnew] with xmid as best point still
	  x1 = xnew;  y1 = ynew;
	}
      }
    iterations_done = it;
    pt = xmid*unitdir + pt;
    return ymid;
  }


float localcostfn(const volume4D<float>& warp)
{
  float cost = globalcostfnptr->cost(warp);
  if (verbose.value()) { cout << "Cost = " << cost << endl; }
  return cost;
}

// 1D minimisation function
float line_minimise(volume4D<float>& warp,
		    const volume4D<float>& gradvec, 
		    const Costfn& costfnobj, float scalefac, 
		    float current_cost)
{
  int niter=1;
  float mincost=1.0;
  mincost = optimise1d(warp, gradvec*scalefac, 0.1, niter, localcostfn,
		       20,  0.0, -1.0); 
  return mincost;
}

////////////////////////////////////////////////////////////////////////////

// this is the non-linear registration bit!
int do_work(int argc, char* argv[]) 
{
  volume<float> vin, vref;
  read_volume(vin,inname.value());
  read_volume(vref,refname.value());

  volume<float> refweight, inweight, outvol;
  refweight = vref * 0.0f + 1.0f;
  inweight = vin * 0.0f + 1.0f;

  Matrix affmat;
  if (affname.set()) {
    affmat = read_ascii_matrix(affname.value());
  } else {
    affmat = Identity(4);
  }

  if (verbose.value()) { cerr << "Setting up costfn object" << endl; }
  Costfn costfnobj(vref,vin,refweight,inweight);
  costfnobj.set_no_bins(nbins.value());
  costfnobj.set_costfn(CorrRatio);
  globalcostfnptr = &costfnobj;

  volume4D<float> warp, bestwarp;
  affine2warp(affmat,warp,vref);
  bestwarp = warp;
  if (initwarp.set()) {
    volume4D<float> initwarpvol;
    read_volume4D(initwarpvol,initwarp.value());
    concat_warps(initwarpvol,warp,bestwarp); 
  }
  warp = bestwarp;
  if (verbose.value()) { print_volume_info(warp,"warp"); }

  if (verbose.value()) { cerr << "Calling cost()" << endl; }
  float cost;
  cost = costfnobj.cost(warp);
  cout << "Non-grad Cost = " << cost << endl;
  volume4D<float> gradvec;

  // start iteration
  for (int iter=1; iter<=maxiter.value(); iter++) { 
    cost = costfnobj.cost_gradient(gradvec,warp,nullbc.value());
    cout << "Grad Cost = " << cost << endl;
    
    // project gradient onto constraint manifold
    gradvec = blur(gradvec,blursize.value());
    if (debug.value()) {
      save_volume4D(gradvec,fslbasename(outname.value())+"_grad");
    }
    
    // calculate scale factor for gradient move
    float scalefac=1.0;
    volume<float> dummy;
    dummy = sumsquaresvol(gradvec);
    dummy = sqrt(dummy);
    scalefac = 1.0 / dummy.percentile(0.95);
    if (verbose.value()) { 
      cout << "scaled gradient motion percentiles (90,95,99,max) are: "
	   << scalefac * dummy.percentile(0.9) << " " 
	   << scalefac * dummy.percentile(0.95) << " "
	   << scalefac * dummy.percentile(0.99) << " " << scalefac * dummy.max() 
	   << endl;
    }

    // clamp the gradvec to avoid very large deformations
    clamp(gradvec,-2.0f/scalefac,+2.0f/scalefac);
    
    cost = line_minimise(warp,gradvec,costfnobj,scalefac,cost);
    bestwarp = warp;

  } // end iteration
  
  warp = bestwarp;

  if (warpname.set()) {
    save_volume4D(warp,warpname.value());
  }

  if (outname.set()) {
    volume<float> vout(vref);
    apply_warp(vin,vout,warp);
    save_volume(vout,outname.value());
  }

  return 0;
}

////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[])
{

  Tracer tr("main");
  OptionParser options(title, examples);

  try {
    // must include all wanted options here (the order determines how
    //  the help message is printed)
    options.add(inname);
    options.add(refname);
    options.add(outname);
    options.add(warpname);
    options.add(affname);
    options.add(initwarp);
    options.add(blursize);
    options.add(maxiter);
    options.add(nbins);
    options.add(nullbc);
    options.add(verbose);
    options.add(debug);
    options.add(help);
    
    nonoptarg = options.parse_command_line(argc, argv);

    // line below stops the program if the help was requested or 
    //  a compulsory option was not set
    if ( (help.value()) || (!options.check_compulsory_arguments(true)) )
      {
	options.usage();
	exit(EXIT_FAILURE);
      }
    
  }  catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  } catch(std::exception &e) {
    cerr << e.what() << endl;
  } 

  // Call the local functions

  return do_work(argc,argv);
}

