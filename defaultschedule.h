// Set default schedule file
//  Written by Mark Jenkinson  11/10/99

#if !defined(__defaultschedule_h)
#define __defaultschedule_h

#include <vector>
#include <string>

void setdefaultschedule(std::vector<string>& comms)
{
  comms.clear();
  comms.push_back("setscale 8");
  comms.push_back("clear S");
  comms.push_back("clear P");
  comms.push_back("search");

  comms.push_back("setscale 4");
  comms.push_back("clear U");
  comms.push_back("optimise 7 S:1-3  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4");
  comms.push_back("optimise 7 P:1-3  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4");
  comms.push_back("optimise 7 P:1-3  1.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4");
  comms.push_back("optimise 7 P:1-3 -1.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4");
  comms.push_back("optimise 7 P:1-3  0.0   1.0   0.0   0.0   0.0   0.0   0.0  rel 4");
  comms.push_back("optimise 7 P:1-3  0.0  -1.0   0.0   0.0   0.0   0.0   0.0  rel 4");
  comms.push_back("optimise 7 P:1-3  0.0   0.0   1.0   0.0   0.0   0.0   0.0  rel 4");
  comms.push_back("optimise 7 P:1-3  0.0   0.0  -1.0   0.0   0.0   0.0   0.0  rel 4");
  comms.push_back("optimise 7 P:1-3  0.0   0.0   0.0   0.0   0.0   0.0   1.0  rel 4");
  comms.push_back("optimise 7 P:1-3  0.0   0.0   0.0   0.0   0.0   0.0  -1.0  rel 4");
  comms.push_back("optimise 7 P:1-3  0.0   0.0   0.0   0.0   0.0   0.0   2.0  rel 4");
  comms.push_back("optimise 7 P:1-3  0.0   0.0   0.0   0.0   0.0   0.0  -2.0  rel 4");
  comms.push_back("sort U");
  comms.push_back("copy U UA");

  comms.push_back("setscale 2");
  comms.push_back("clear U");
  comms.push_back("optimise 7  UA:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4");
  comms.push_back("if MAXDOF > 7");
  comms.push_back("optimise 9  UA:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 1");
  comms.push_back("if MAXDOF > 9");
  comms.push_back("optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 2");
  comms.push_back("sort U");
  comms.push_back("copy U UB");

  comms.push_back("setscale 1");
  comms.push_back("clear U");
  comms.push_back("optimise 12 UB:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 1");
  comms.push_back("sort U");

}

#endif
