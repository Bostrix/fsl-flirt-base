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
#include "newimage/imfft.h"
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
Matrix initaffmat;

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

// topology preservation code

void jacobian_check(volume4D<float>& jvol,
		    ColumnVector& jacobian_stats, 
		    const volume4D<float>& warp,
		    float minJ, float maxJ, bool use_vol=true)
{
  // set up jacobian stats to contain: min, max, num < minJ, num > maxJ
  if (jacobian_stats.Nrows()!=4) { jacobian_stats.ReSize(4); }
  jacobian_stats = 0.0;
  jacobian_stats(1)=1.0; jacobian_stats(2)=1.0; 
  if (use_vol) {
    if (!samesize(jvol[0],warp[0]) || (jvol.tsize()!=8)) {
      jvol = warp;  // set up all the right properties
      jvol=0.0f; 
      for (int n=1; n<=5; n++) { jvol.addvolume(jvol[0]); }
    }
  }
  float Jfff, Jbff, Jfbf, Jffb, Jbbf, Jbfb, Jfbb, Jbbb;
  float wx000=0,wx001=0,wx010=0,wx011=0,wx100=0,wx101=0,wx110=0,wx111=0;
  float wy000=0,wy001=0,wy010=0,wy011=0,wy100=0,wy101=0,wy110=0,wy111=0;
  float wz000=0,wz001=0,wz010=0,wz011=0,wz100=0,wz101=0,wz110=0,wz111=0;
  float volscale=1.0/(warp.xdim() * warp.ydim() * warp.zdim());
  for (int z=warp.minz(); z<=warp.maxz()-1; z++) {
    for (int y=warp.miny(); y<=warp.maxy()-1; y++) {
      for (int x=warp.minx(); x<=warp.maxx()-1; x++) {
	warp[0].getneighbours(x,y,z,wx000,wx001,wx010,wx011,
			      wx100,wx101,wx110,wx111);
	warp[1].getneighbours(x,y,z,wy000,wy001,wy010,wy011,
			      wy100,wy101,wy110,wy111);
	warp[2].getneighbours(x,y,z,wz000,wz001,wz010,wz011,
			      wz100,wz101,wz110,wz111);
	Jfff = (wx100-wx000) *
	  ((wy010-wy000) * (wz001-wz000) - (wy001-wy000) * (wz010-wz000)) 
	  - (wx010-wx000) *
	  ((wy100-wy000) * (wz001-wz000) - (wy001-wy000) * (wz100-wz000)) 
	  + (wx001-wx000) *
	  ((wy100-wy000) * (wz010-wz000) - (wy010-wy000) * (wz100-wz000));
	Jfff *= volscale;
	if (Jfff<jacobian_stats(1)) { jacobian_stats(1)=Jfff; }
	if (Jfff>jacobian_stats(2)) { jacobian_stats(2)=Jfff; }
	if (Jfff<minJ) { jacobian_stats(3)+=1.0; }
	if (Jfff>maxJ) { jacobian_stats(4)+=1.0; }
	Jbff = (wx100-wx000) *
	  ((wy110-wy100) * (wz101-wz100) - (wy101-wy100) * (wz110-wz100)) 
	  - (wx110-wx100) *
	  ((wy100-wy000) * (wz101-wz100) - (wy101-wy100) * (wz100-wz000)) 
	  + (wx101-wx100) *
	  ((wy100-wy000) * (wz110-wz100) - (wy110-wy100) * (wz100-wz000));
	Jbff *= volscale;
	if (Jbff<jacobian_stats(1)) { jacobian_stats(1)=Jbff; }
	if (Jbff>jacobian_stats(2)) { jacobian_stats(2)=Jbff; }
	if (Jbff<minJ) { jacobian_stats(3)+=1.0; }
	if (Jbff>maxJ) { jacobian_stats(4)+=1.0; }
	Jfbf = (wx110-wx010) *
	  ((wy010-wy000) * (wz011-wz010) - (wy011-wy010) * (wz010-wz000)) 
	  - (wx010-wx000) *
	  ((wy110-wy010) * (wz011-wz010) - (wy011-wy010) * (wz110-wz010)) 
	  + (wx011-wx010) *
	  ((wy110-wy010) * (wz010-wz000) - (wy010-wy000) * (wz110-wz010));
	Jfbf *= volscale;
	if (Jfbf<jacobian_stats(1)) { jacobian_stats(1)=Jfbf; }
	if (Jfbf>jacobian_stats(2)) { jacobian_stats(2)=Jfbf; }
	if (Jfbf<minJ) { jacobian_stats(3)+=1.0; }
	if (Jfbf>maxJ) { jacobian_stats(4)+=1.0; }
	Jffb = (wx101-wx001) *
	  ((wy011-wy001) * (wz001-wz000) - (wy001-wy000) * (wz011-wz001)) 
	  - (wx011-wx001) *
	  ((wy101-wy001) * (wz001-wz000) - (wy001-wy000) * (wz101-wz001)) 
	  + (wx001-wx000) *
	  ((wy101-wy001) * (wz011-wz001) - (wy011-wy001) * (wz101-wz001));
	Jffb *= volscale;
	if (Jffb<jacobian_stats(1)) { jacobian_stats(1)=Jffb; }
	if (Jffb>jacobian_stats(2)) { jacobian_stats(2)=Jffb; }
	if (Jffb<minJ) { jacobian_stats(3)+=1.0; }
	if (Jffb>maxJ) { jacobian_stats(4)+=1.0; }
	Jfbb = (wx111-wx011) *
	  ((wy011-wy001) * (wz011-wz010) - (wy011-wy010) * (wz011-wz001)) 
	  - (wx011-wx001) *
	  ((wy111-wy011) * (wz011-wz010) - (wy011-wy010) * (wz111-wz011)) 
	  + (wx011-wx010) *
	  ((wy111-wy011) * (wz011-wz001) - (wy011-wy001) * (wz111-wz011));
	Jfbb *= volscale;
	if (Jfbb<jacobian_stats(1)) { jacobian_stats(1)=Jfbb; }
	if (Jfbb>jacobian_stats(2)) { jacobian_stats(2)=Jfbb; }
	if (Jfbb<minJ) { jacobian_stats(3)+=1.0; }
	if (Jfbb>maxJ) { jacobian_stats(4)+=1.0; }
	Jbfb = (wx101-wx001) *
	  ((wy111-wy101) * (wz101-wz100) - (wy101-wy100) * (wz111-wz101)) 
	  - (wx111-wx101) *
	  ((wy101-wy001) * (wz101-wz100) - (wy101-wy100) * (wz101-wz001)) 
	  + (wx101-wx100) *
	  ((wy101-wy001) * (wz111-wz101) - (wy111-wy101) * (wz101-wz001));
	Jbfb *= volscale;
	if (Jbfb<jacobian_stats(1)) { jacobian_stats(1)=Jbfb; }
	if (Jbfb>jacobian_stats(2)) { jacobian_stats(2)=Jbfb; }
	if (Jbfb<minJ) { jacobian_stats(3)+=1.0; }
	if (Jbfb>maxJ) { jacobian_stats(4)+=1.0; }
	Jbbf = (wx110-wx010) *
	  ((wy110-wy100) * (wz111-wz110) - (wy111-wy110) * (wz110-wz100)) 
	  - (wx110-wx100) *
	  ((wy110-wy010) * (wz111-wz110) - (wy111-wy110) * (wz110-wz010)) 
	  + (wx111-wx110) *
	  ((wy110-wy010) * (wz110-wz100) - (wy110-wy100) * (wz110-wz010));
	Jbbf *= volscale;
	if (Jbbf<jacobian_stats(1)) { jacobian_stats(1)=Jbbf; }
	if (Jbbf>jacobian_stats(2)) { jacobian_stats(2)=Jbbf; }
	if (Jbbf<minJ) { jacobian_stats(3)+=1.0; }
	if (Jbbf>maxJ) { jacobian_stats(4)+=1.0; }
	Jbbb = (wx111-wx011) *
	  ((wy111-wy101) * (wz111-wz110) - (wy111-wy110) * (wz111-wz101)) 
	  - (wx111-wx101) *
	  ((wy111-wy011) * (wz111-wz110) - (wy111-wy110) * (wz111-wz011)) 
	  + (wx111-wx110) *
	  ((wy111-wy011) * (wz111-wz101) - (wy111-wy101) * (wz111-wz011));
	Jbbb *= volscale;
	if (Jbbb<jacobian_stats(1)) { jacobian_stats(1)=Jbbb; }
	if (Jbbb>jacobian_stats(2)) { jacobian_stats(2)=Jbbb; }
	if (Jbbb<minJ) { jacobian_stats(3)+=1.0; }
	if (Jbbb>maxJ) { jacobian_stats(4)+=1.0; }
// 	if (debug.value()) {
// 	  cout << "Jacobian values are: " << Jfff << ", " << Jffb << ", " 
// 	       << Jfbf << ", " << Jbff << ", " << Jfbb << ", " << Jbfb 
// 	       << ", " << Jbbf << "," << Jbbb  << endl;
// 	}
	if (use_vol) {
	  // the following must be consistent with get_jac_offset()
	  jvol(x,y,z,0) = Jfff;
	  jvol(x,y,z,1) = Jbff;
	  jvol(x,y,z,2) = Jfbf;
	  jvol(x,y,z,3) = Jffb;
	  jvol(x,y,z,4) = Jfbb;
	  jvol(x,y,z,5) = Jbfb;
	  jvol(x,y,z,6) = Jbbf;
	  jvol(x,y,z,7) = Jbbb;
	}
      }
    }
  }
}


volume4D<float> jacobian_check(ColumnVector& jacobian_stats, 
			       const volume4D<float>& warp,
			       float minJ, float maxJ)
{
  volume4D<float> jvol;
  jacobian_check(jvol,jacobian_stats,warp,minJ,maxJ,true);
  return jvol;
}


void jacobian_check_quick(ColumnVector& jacobian_stats, 
			  const volume4D<float>& warp,
			  float minJ, float maxJ)
{
  volume4D<float> dummy;
  jacobian_check(dummy,jacobian_stats,warp,minJ,maxJ,false);
}

void grad_calc(volume4D<float>& gradvols, const volume4D<float>& warp)
{
  // returns gradients in the order: dx'/dx, dx'/dy, dx'/dz, dy'/dx, etc
  if (!samesize(gradvols[0],warp[0]) || (gradvols.tsize()!=9)) {
    gradvols = warp;  // set up all the right properties
    gradvols=0.0f; 
    for (int n=1; n<=6; n++) { gradvols.addvolume(gradvols[0]); }
  }
  float wx000=0,wx001=0,wx010=0,wx011=0,wx100=0,wx101=0,wx110=0,wx111=0;
  float wy000=0,wy001=0,wy010=0,wy011=0,wy100=0,wy101=0,wy110=0,wy111=0;
  float wz000=0,wz001=0,wz010=0,wz011=0,wz100=0,wz101=0,wz110=0,wz111=0;
  for (int z=warp.minz(); z<=warp.maxz()-1; z++) {
    for (int y=warp.miny(); y<=warp.maxy()-1; y++) {
      for (int x=warp.minx(); x<=warp.maxx()-1; x++) {
	warp[0].getneighbours(x,y,z,wx000,wx001,wx010,wx011,
			      wx100,wx101,wx110,wx111);
	warp[1].getneighbours(x,y,z,wy000,wy001,wy010,wy011,
			      wy100,wy101,wy110,wy111);
	warp[2].getneighbours(x,y,z,wz000,wz001,wz010,wz011,
			      wz100,wz101,wz110,wz111);
	gradvols[0](x,y,z) = wx100-wx000;
	gradvols[1](x,y,z) = wx010-wx000;
	gradvols[2](x,y,z) = wx001-wx000;
	gradvols[3](x,y,z) = wy100-wy000;
	gradvols[4](x,y,z) = wy010-wy000;
	gradvols[5](x,y,z) = wy001-wy000;
	gradvols[6](x,y,z) = wz100-wz000;
	gradvols[7](x,y,z) = wz010-wz000;
	gradvols[8](x,y,z) = wz001-wz000;
      }
    }
  }
}


void integrate_gradient_field(volume4D<float>& newwarp, 
			      const volume4D<float>& grad,
			      float warpmeanx, float warpmeany, float warpmeanz)
{
  // enforces integrability constraints and returns the integrated grad field
  // Note that the mean of the newwarp will be equal to warpmean{x,y,z}
  //  pass in: oldwarp[0].mean(), oldwarp[1].mean(), oldwarp[2].mean()
  
  int Nx, Ny, Nz;
  Nx = grad.maxx();
  Ny = grad.maxy();
  Nz = grad.maxz();
  if (!samesize(newwarp[0],grad[0]) || (newwarp.tsize()!=3) ) {
    newwarp = grad;
    for (int n=8; n>2; n++) { newwarp.deletevolume(n); }
    newwarp = 0.0f;
  }
  volume4D<float> gradkre(newwarp), gradkim(newwarp);
  float dotprodre, dotprodim, norm, argx, argy, argz;
  float gradkrealx, gradkimagx, gradkrealy, gradkimagy, gradkrealz, gradkimagz;
  // enforce things separately for gradients of warp[0], warp[1] and warp[2]
  for (int n=0; n<3; n++) {
    // take FFT of the x,y,z gradient fields (of warp[n])
    fft3(grad[n*3+0],grad[n*3+0]*0.0f,gradkre[0],gradkim[0]);
    fft3(grad[n*3+1],grad[n*3+1]*0.0f,gradkre[1],gradkim[1]);
    fft3(grad[n*3+2],grad[n*3+2]*0.0f,gradkre[2],gradkim[2]);
    // take normalised dot product of gradient vector and "A" vector
    for (int z=grad.minz(); z<=grad.maxz()-1; z++) {
      for (int y=grad.miny(); y<=grad.maxy()-1; y++) {
	for (int x=grad.minx(); x<=grad.maxx()-1; x++) {
	  argx=2.0*M_PI*x/Nx;  argy=2.0*M_PI*y/Ny;  argz=2.0*M_PI*z/Nz;  
	  norm = 6.0 - 2.0*cos(argx) - 2.0*cos(argy) - 2.0*cos(argz);
	  gradkrealx = gradkre[0](x,y,z);
	  gradkimagx = gradkim[0](x,y,z);
	  gradkrealy = gradkre[1](x,y,z);
	  gradkimagy = gradkim[1](x,y,z);
	  gradkrealz = gradkre[2](x,y,z);
	  gradkimagz = gradkim[2](x,y,z);
	  dotprodre = 0.0;  dotprodim = 0.0;
	  dotprodre += gradkrealx * (cos(argx)-1) + gradkimagx * sin(argx);
	  dotprodim += gradkimagx * (cos(argx)-1) - gradkrealx * sin(argx);
	  dotprodre += gradkrealy * (cos(argy)-1) + gradkimagy * sin(argy);
	  dotprodim += gradkimagy * (cos(argy)-1) - gradkrealy * sin(argy);
	  dotprodre += gradkrealz * (cos(argz)-1) + gradkimagz * sin(argz);
	  dotprodim += gradkimagz * (cos(argz)-1) - gradkrealz * sin(argz);
	  // write back values into gradkre[0] and gradkim[0]
	  gradkre[0](x,y,z) = dotprodre / norm;
	  gradkim[0](x,y,z) = dotprodim / norm;
	}
      }
    }
    // take IFFT to get the integrated gradient field
    ifft3(gradkre[0],gradkim[0]);
    newwarp[n] = gradkre[0];
  }
  // adjust the mean values
  newwarp[0] += warpmeanx;
  newwarp[1] += warpmeany;
  newwarp[2] += warpmeanz;
}

void get_jac_offset(int jacnum, int& xoff, int& yoff, int& zoff)
{
  xoff=0; yoff=0; zoff=0;
  if (jacnum==0) { xoff=0; yoff=0; zoff=0; } // Jfff
  if (jacnum==1) { xoff=1; yoff=0; zoff=0; } // Jbff
  if (jacnum==2) { xoff=0; yoff=1; zoff=0; } // Jfbf
  if (jacnum==3) { xoff=0; yoff=0; zoff=1; } // Jffb
  if (jacnum==4) { xoff=0; yoff=1; zoff=1; } // Jfbb
  if (jacnum==5) { xoff=1; yoff=0; zoff=1; } // Jbfb
  if (jacnum==6) { xoff=1; yoff=1; zoff=0; } // Jbbf
  if (jacnum==7) { xoff=1; yoff=1; zoff=1; } // Jbbb
}

void limit_grad(volume4D<float>& grad, const volume4D<float>& jvol, 
	   float minJ, float maxJ)
{
  // use initaffmat for the default config (needs to be modified if
  //  this becomes a library function)
  Matrix J, Jnew, J0;
  J0 = initaffmat;
  float alpha, detJ;
  int xoff, yoff, zoff;
  for (int z=grad.minz(); z<=grad.maxz()-1; z++) {
    for (int y=grad.miny(); y<=grad.maxy()-1; y++) {
      for (int x=grad.minx(); x<=grad.maxx()-1; x++) {
	for (jacnum=0; jacnum<8; jacnum++) {
	  if ((jvol[jacnum](x,y,z)<minJ) || (jvol[jacnum](x,y,z)>maxJ)) {
	    get_jac_offset(jacnum,xoff,yoff,zoff);
	    // interpolate between current J matrix and initaffmat
	    for (int n1=1; n1<=3; n1++) { for (int n2=1; n2<=3; n2++) {
	      J(n1,n2) = grad[(n1-1)*3 + (n2-1)](x+xoff,y+yoff,z+zoff);
	    } }
	    alpha = 0.0;
	    Jnew = J;
	    detJ = Jnew.Determinant();
	    while ( (detJ > maxJ) || (detJ < minJ) ) {
	      alpha += 0.1;
	      if (alpha>1.0) alpha=1.0;
	      Jnew = (1 - alpha ) * J + alpha * J0;
	      detJ = Jnew.Determinant();
	    }
	    // rescale gradients as required
	    for (int n1=1; n1<=3; n1++) { for (int n2=1; n2<=3; n2++) {
	      grad[(n1-1)*3 + (n2-1)](x+xoff,y+yoff,z+zoff) = 
		(1 - alpha) * J(n1,n2) + alpha * J0(n1,n2);
	    } }
	  }
	}
      }
    }
  }
}


void constrain_topology(volume4D<float>& warp, float minJ, float maxJ)
{
  ColumnVector jstats(4);
  jacobian_check_quick(jstats,warp,minJ,maxJ);
  volume4D<float> grad, jvol;
  int n=1, maxit=10;
  while ( (n++<maxit) && ( (jstats(3)>0.5) || (jstats(4)>0.5) ) ) {
    grad_calc(grad,warp);
    jacobian_check(jvol,jstats,warp,minJ,maxJ);
    limit_grad(grad,jvol,minJ,maxJ);
    integrate_gradient_field(warp, grad, warp[0].mean(), warp[1].mean(), 
			     warp[2].mean());
    jacobian_check_quick(jstats,warp,minJ,maxJ);
  }
}

void constrain_topology(volume4D<float>& warp)
{
  constrain_topology(warp,0.01,100.0);  // mainly just enforcing positivity
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

  if (affname.set()) {
    initaffmat = read_ascii_matrix(affname.value());
  } else {
    initaffmat = Identity(4);
  }

  if (verbose.value()) { cerr << "Setting up costfn object" << endl; }
  Costfn costfnobj(vref,vin,refweight,inweight);
  costfnobj.set_no_bins(nbins.value());
  costfnobj.set_costfn(CorrRatio);
  globalcostfnptr = &costfnobj;

  volume4D<float> warp, bestwarp;
  affine2warp(initaffmat,warp,vref);
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

  volume4D<float> jvol;
  ColumnVector jstats(4);
  jstats=0;
  jvol = jacobian_check(jstats,bestwarp,0.1,7.0);
  if (verbose.value()) { 
    cout << "Jacobian stats: min J = " << jstats(1) << " , max J = " <<
      jstats(2) << " , num < 0.1 = " << jstats(3) << " , num > 7.0 = " <<
      jstats(4) << endl;
  }
  if (debug.value()) {
    save_volume4D(jvol,fslbasename(outname.value())+"_jac");
  }

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
    jvol = jacobian_check(jstats,bestwarp,0.1,7.0);
    if (verbose.value()) { 
      cout << "Jacobian stats: min J = " << jstats(1) << " , max J = " <<
	jstats(2) << " , num < 0.1 = " << jstats(3) << " , num > 7.0 = " <<
	jstats(4) << endl;
    }
    if (debug.value()) {
      save_volume4D(jvol,fslbasename(outname.value())+"_jac");
    }

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

