/*  rmsdiff.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */


#include <string>
#include <iostream>
#include <unistd.h>
#define WANT_STREAM
#define WANT_MATH

#include "newmatap.h"
#include "newmatio.h"
#include "miscmaths/miscmaths.h"
#include "newimage/newimageall.h"

#ifndef NO_NAMESPACE
 using namespace MISCMATHS;
 using namespace NEWMAT;
 using namespace NEWIMAGE;
#endif


////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[])
{

  float rmax=80.0;

  if (argc!=4) { 
    cerr << "Usage: " << argv[0] << " matrixfile1 matrixfile2 refvol" << endl; 
    cerr << "        Outputs rms deviation between matrices (in mm)" << endl;
    return -1; 
  }
  
  Matrix affmat1(4,4), affmat2(4,4);
  affmat1 = read_ascii_matrix(argv[1]);
  if (affmat1.Nrows()<4) {
    cerr << "Could not read matrix " << argv[1] << endl;
    return -2;
  }
  affmat2 = read_ascii_matrix(argv[2]);
  if (affmat2.Nrows()<4) {
    cerr << "Could not read matrix " << argv[2] << endl;
    return -2;
  }

  if (fabs(affmat1.Determinant())<0.1) {
    cerr << "WARNING:: matrix 1 has low determinant" << endl;
    cerr << affmat1 << endl;
  }

  if (fabs(affmat2.Determinant())<0.1) {
    cerr << "WARNING:: matrix 2 has low determinant" << endl;
    cerr << affmat2 << endl;
  }

  ColumnVector centre(3);
  centre = 0;

  volume<float> refvol;
  read_volume_hdr_only(refvol,argv[3]);
  // compute the centre of volume (in world coordinates)
  centre(1) = 0.5*(refvol.xsize() - 1.0)*refvol.xdim();
  centre(2) = 0.5*(refvol.ysize() - 1.0)*refvol.ydim();
  centre(3) = 0.5*(refvol.zsize() - 1.0)*refvol.zdim();

  float rms = rms_deviation(affmat1,affmat2,centre,rmax);

  cout << rms << endl;
  
  return 0;

}





