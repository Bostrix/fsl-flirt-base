/*  defaultschedule.h

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */

// Set default schedule file
//  Written by Mark Jenkinson  11/10/99

#if !defined(__defaultschedule_h)
#define __defaultschedule_h

#include <vector>
#include <string>

void setdefaultschedule(std::vector<string>& comms)
{
  comms.clear();
  comms.push_back("# 8mm scale");
  comms.push_back("setscale 8");
  comms.push_back("setoption smoothing 8");
  comms.push_back("clear S");
  comms.push_back("clear P");
  comms.push_back("search");

  comms.push_back("# 4mm scale");
  comms.push_back("setscale 4");
  comms.push_back("setoption smoothing 4");
  comms.push_back("clear U");
  comms.push_back("clear UA ");
  comms.push_back("clear UB");
  comms.push_back("clear US");
  comms.push_back("clear UP");

  comms.push_back("# remeasure costs at this scale");
  comms.push_back("measurecost 7 S 0 0 0 0 0 0 rel");
  comms.push_back("copy U US");
  comms.push_back("clear U");
  comms.push_back("measurecost 7 P 0 0 0 0 0 0 rel");
  comms.push_back("copy U UP");
  comms.push_back("dualsort US UP");

  comms.push_back("# optimise best 3 candidates (pre and post 8mm optimisations)");
  comms.push_back("clear U");
  comms.push_back("optimise 7 US:1-3  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4");
  comms.push_back("optimise 7 UP:1-3  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4");
  comms.push_back("sort U");
  comms.push_back("copy U UA");

  comms.push_back("# select best 3 optimised solutions and try perturbations of these");
  comms.push_back("clear U");
  comms.push_back("copy UA:1-3 U");
  comms.push_back("optimise 7 UA:1-3  1.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4");
  comms.push_back("optimise 7 UA:1-3 -1.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4");
  comms.push_back("optimise 7 UA:1-3  0.0   1.0   0.0   0.0   0.0   0.0   0.0  rel 4");
  comms.push_back("optimise 7 UA:1-3  0.0  -1.0   0.0   0.0   0.0   0.0   0.0  rel 4");
  comms.push_back("optimise 7 UA:1-3  0.0   0.0   1.0   0.0   0.0   0.0   0.0  rel 4");
  comms.push_back("optimise 7 UA:1-3  0.0   0.0  -1.0   0.0   0.0   0.0   0.0  rel 4");
  comms.push_back("optimise 7 UA:1-3  0.0   0.0   0.0   0.0   0.0   0.0   0.1  abs 4");
  comms.push_back("optimise 7 UA:1-3  0.0   0.0   0.0   0.0   0.0   0.0  -0.1  abs 4");
  comms.push_back("optimise 7 UA:1-3  0.0   0.0   0.0   0.0   0.0   0.0   0.2  abs 4");
  comms.push_back("optimise 7 UA:1-3  0.0   0.0   0.0   0.0   0.0   0.0  -0.2  abs 4");
  comms.push_back("sort U");
  comms.push_back("copy U UB");

  comms.push_back("# 2mm scale");
  comms.push_back("setscale 2");
  comms.push_back("setoption smoothing 2");
  comms.push_back("clear U");
  comms.push_back("clear UC");
  comms.push_back("clear UD");
  comms.push_back("clear UE");
  comms.push_back("clear UF");

  comms.push_back("# remeasure costs at this scale");
  comms.push_back("measurecost 7 UB 0 0 0 0 0 0 rel");
  comms.push_back("sort U");
  comms.push_back("copy U UC");

  comms.push_back("clear U");
  comms.push_back("optimise 7  UC:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4");
  comms.push_back("copy U UD");
  comms.push_back("setoption boundguess 1");
  comms.push_back("if MAXDOF > 7");
  comms.push_back(" clear U");
  comms.push_back("if MAXDOF > 7");
  comms.push_back(" optimise 9  UD:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 1");
  comms.push_back("copy U UE");
  comms.push_back("if MAXDOF > 9");
  comms.push_back(" clear U");
  comms.push_back("if MAXDOF > 9");
  comms.push_back(" optimise 12 UE:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 2");
  comms.push_back("sort U");
  comms.push_back("copy U UF");

  comms.push_back("# 1mm scale");
  comms.push_back("setscale 1");
  comms.push_back("setoption smoothing 1");
  comms.push_back("setoption boundguess 1");
  comms.push_back("clear U");
  comms.push_back("optimise 12 UF:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 1");
  comms.push_back("sort U");
}

#endif
