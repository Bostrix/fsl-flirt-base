#include <string>
#include <iostream>
#include <unistd.h>
#define WANT_STREAM
#define WANT_MATH

#include "newmatap.h"
#include "newmatio.h"
#include "miscmaths.h"
#include "generalio.h"

#ifndef NO_NAMESPACE
 using namespace MISCMATHS;
 using namespace NEWMAT;
 using namespace GENERALIO;
#endif


////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[])
{

  if (argc!=3) { 
    cerr << "Usage: " << argv[0] << " matrixfile1 matrixfile2\n"; 
    return -1; 
  }
  
  Matrix affmat1(4,4), affmat2(4,4);
  if (read_ascii_matrix(affmat1,argv[1])<0)   return -2;
  if (read_ascii_matrix(affmat2,argv[2])<0)   return -2;

  float rms = rms_deviation(affmat1,affmat2);

  cout << rms << endl;
  
  return 0;

}








