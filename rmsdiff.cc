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

  float rmax=80.0;

  if ((argc<3) || (argc>4)) { 
    cerr << "Usage: " << argv[0] << " matrixfile1 matrixfile2 [refvol]\n"; 
    return -1; 
  }
  
  Matrix affmat1(4,4), affmat2(4,4);
  if (read_ascii_matrix(affmat1,argv[1])<0)   return -2;
  if (read_ascii_matrix(affmat2,argv[2])<0)   return -2;

  ColumnVector centre(3);
  centre = 0;

  if (argc==4) {
    volume refvol;
    read_volume_hdr_only(refvol,argv[3]);
    // compute the centre of volume (in world coordinates)
    centre(1) = 0.5*(refvol.xsize() - 1.0)*refvol.getx();
    centre(2) = 0.5*(refvol.ysize() - 1.0)*refvol.gety();
    centre(3) = 0.5*(refvol.zsize() - 1.0)*refvol.getz();
  }

  float rms = rms_deviation(affmat1,affmat2,centre,rmax);

  cout << rms << endl;
  
  return 0;

}





