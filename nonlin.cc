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

  volume4D<float> warp;
  affine2warp(Identity(4),warp,vref);

  if (verbose.value()) { cerr << "Setting up costfn object" << endl; }
  Costfn costfnobj(vref,vin,refweight,inweight);
  costfnobj.set_no_bins(256);
  costfnobj.set_costfn(CorrRatio);

  float cost;
  if (verbose.value()) { cerr << "Calling cost()" << endl; }
  if (verbose.value()) { print_volume_info(warp,"warp"); }

  cerr << "Pre-calling affine cost" << endl;
  cost = costfnobj.cost(Identity(4));
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
  cost = costfnobj.cost(warp);
  cout << "Cost = " << cost << endl;

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
    options.add(verbose);
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

