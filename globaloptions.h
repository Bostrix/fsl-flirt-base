/*  globaloptions.h

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */

#ifndef __GLOBALOPTIONS_
#define __GLOBALOPTIONS_

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "costfns.h"

using namespace COSTFNS;

  enum anglereps { Euler, Quaternion };
  enum interps { TriLinear, NearestNeighbour, Sinc };
  enum windowtype { Rectangular, Hanning, Blackman };

  typedef std::vector<RowVector> MatVec;
  typedef MatVec* MatVecPtr;




class globaloptions {
 public:
  static globaloptions& get();
  ~globaloptions() { delete gopt; }
  
  string version;
  
  std::vector<MatVec> usrmat;
  MatVec searchoptmat;
  MatVec preoptsearchmat;

  string inputfname;
  string outputfname;
  string reffname;
  string outputmatascii;
  string outputmatmedx;
  string initmatfname;
  string refweightfname;
  string testweightfname;
  Matrix initmat;

  string schedulefname;

  imagepair *impair;
  ColumnVector refparams;
  Matrix parammask;
  int no_params;
  int dof;
  bool usrsubset;
  int searchdof;
  int no_bins;
  costfns maincostfn;
  costfns searchcostfn;
  costfns currentcostfn;
  anglereps anglerep;
  float isoscale;
  float min_sampling;
  int lastsampling;
  float requestedscale;
  bool force_scaling;
  float smoothsize;
  float fuzzyfrac;
  ColumnVector tolerance;
  ColumnVector boundguess;

  ColumnVector searchrx;
  ColumnVector searchry;
  ColumnVector searchrz;
  float coarsedelta;
  float finedelta;

  short datatype;
  bool forcedatatype;
  int verbose;
  bool interactive;
  bool do_optimise;
  bool nosave;
  bool iso;
  bool resample;
  bool useweights;
  bool mode2D;
  interps interpmethod;
  int sincwidth;
  windowtype sincwindow;

  int single_param;


  void parse_command_line(int argc, char** argv, const string &);

 private:
  globaloptions();
  
  const globaloptions& operator=(globaloptions&);
  globaloptions(globaloptions&);
      
  static globaloptions* gopt;

  void print_usage(int argc, char *argv[]);  
  void print_version();
  
};



//--------------------------------------------------------------------------//


inline globaloptions& globaloptions::get(){
  if(gopt == NULL)
    gopt = new globaloptions();
  
  return *gopt;
}

inline globaloptions::globaloptions()
{
  // set up defaults

  version = "";

  searchoptmat.clear();
  preoptsearchmat.clear();
  usrmat.resize(27);
  for (unsigned int i=0; i<usrmat.size(); i++) {
    usrmat[i].clear();
  }

  outputfname = "";
  reffname = "";

  inputfname = "";
  outputmatascii = "";
  outputmatmedx = "";
  initmatfname = "";
  refweightfname = "";
  testweightfname = "";
  initmat = Identity(4);

  schedulefname = "";
  
  impair = 0;
  refparams.ReSize(12);
  refparams << 0.0 << 0.0 << 0.0 << 0.0 << 0.0 << 0.0 << 1.0 << 1.0 << 1.0
	    << 0.0 << 0.0 << 0.0;
  parammask.ReSize(12,12);
  Identity(parammask);
  no_params = 12;
  dof = 12;
  usrsubset = false;
  searchdof = 12;
  no_bins = 256;
  maincostfn = CorrRatio;
  searchcostfn = CorrRatio;
  currentcostfn = CorrRatio;
  anglerep = Euler;
  isoscale = 1.0;
  min_sampling = 1.0;
  lastsampling = 8;
  requestedscale = 1.0;
  force_scaling = false;
  smoothsize = 1.0;
  fuzzyfrac = 0.5;
  tolerance.ReSize(12);
  tolerance << 0.005 << 0.005 << 0.005 << 0.2 << 0.2 << 0.2 << 0.002 
	    << 0.002 << 0.002 << 0.001 << 0.001 << 0.001; 
  boundguess.ReSize(2);
  boundguess << 10.0 << 1.0;

  searchrx.ReSize(2);
  searchry.ReSize(2);
  searchrz.ReSize(2);
  searchrx << -M_PI/2.0 << M_PI/2.0;
  searchry << -M_PI/2.0 << M_PI/2.0;
  searchrz << -M_PI/2.0 << M_PI/2.0;
  coarsedelta = 60.0*M_PI/180.0;
  finedelta = 18.0*M_PI/180.0;

  datatype = -1;
  forcedatatype = false;
  verbose = 0;
  interactive = false;
  do_optimise = true;
  nosave = true;
  iso = true;
  resample = true;
  useweights = false;
  mode2D = false;
  interpmethod = TriLinear;
  sincwidth = 3;
  sincwindow = Hanning;

  single_param = -1;
}

#endif

