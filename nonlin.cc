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
int nonoptarg;

////////////////////////////////////////////////////////////////////////////

// Local functions

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

  volume4D<float> warp;
  affine2warp(affmat,warp,vref);

  if (verbose.value()) { cerr << "Setting up costfn object" << endl; }
  Costfn costfnobj(vref,vin,refweight,inweight);
  costfnobj.set_no_bins(256);
  costfnobj.set_costfn(CorrRatio);

  float cost;
  if (verbose.value()) { cerr << "Calling cost()" << endl; }
  if (verbose.value()) { print_volume_info(warp,"warp"); }

  cerr << "Pre-calling affine cost" << endl;
  cost = costfnobj.cost(affmat);
  cerr << "Returned affine cost = " << cost << endl;

  cost = costfnobj.cost(warp);
  cout << "Cost = " << cost << endl;

  for (int z=vref.minz(); z<=vref.maxz(); z++) {
    for (int y=vref.miny(); y<=vref.maxy(); y++) {
      for (int x=vref.minx(); x<=vref.maxx(); x++) {
	warp[1](x,y,z) += 0.001*x*x + z/200;
      }
    }
  }
  cost = costfnobj.cost(warp);
  cout << "Cost = " << cost << endl;

  for (int z=vref.minz(); z<=vref.maxz(); z++) {
    for (int y=vref.miny(); y<=vref.maxy(); y++) {
      for (int x=vref.minx(); x<=vref.maxx(); x++) {
	warp[1](x,y,z) += 0.01*x*x + z/20;
      }
    }
  }
  
  cout << "Returning to original warp" << endl;
  affine2warp(affmat,warp,vref);
  cost = costfnobj.cost(warp);
  cout << "Non-grad Cost = " << cost << endl;
  volume4D<float> gradvec;
  cost = costfnobj.cost_gradient(gradvec,warp);
  cout << "Grad Cost = " << cost << endl;

  // project gradient onto constraint manifold
  float blursize=20.0; // mm
  gradvec[0] = blur(gradvec[0],blursize);
  gradvec[1] = blur(gradvec[1],blursize);
  gradvec[2] = blur(gradvec[2],blursize);

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

  ColumnVector sfacs(7);
  sfacs << -0.5 << -1.0 << -1.0 << +2.5 << 0.5 << 1.0 << 1.0;
  sfacs *= scalefac;

  float cumfac = 0.0, bestcumfac=0, mincost=cost;
  for (int idx=1; idx<=sfacs.Nrows(); idx++) {
    cumfac += sfacs(idx);
    warp += ((float) sfacs(idx))*gradvec;
    cost = costfnobj.cost(warp);
    if (cost<mincost) { mincost=cost; bestcumfac=cumfac; }
    cout << "Post-grad (" << cumfac/scalefac << ") Cost = " << cost << endl;
  }

  // reset warp to best one found so far
  affine2warp(affmat,warp,vref);
  warp += bestcumfac*gradvec;

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

