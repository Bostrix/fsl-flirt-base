/*  avscale.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */

#include <string>
#include <iostream>
#include <fstream>
#include <unistd.h>
#define WANT_STREAM
#define WANT_MATH

#include "newmatap.h"
#include "newmatio.h"
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"

#ifndef NO_NAMESPACE
 using namespace MISCMATHS;
 using namespace NEWMAT;
 using namespace NEWIMAGE;
#endif


////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[])
{

  try {

    if (argc<3) { 
      cerr << "Usage: " << argv[0] << " matrixfile reslice_volume\n"; 
      return -1; 
    }

    volume<float> testvol;
    if (read_volume_hdr_only(testvol,argv[2])<0)  return -1;
    Matrix affmat(4,4);
    ColumnVector params(12);
    if (read_matrix(affmat,argv[1],testvol)<0)   return -2;
    //cout << affmat << endl;
    decompose_aff(params,affmat,rotmat2euler);

    Matrix rotmat(4,4);
    construct_rotmat_euler(params,6,rotmat);
    
    cout << "Rotation & Translation Matrix:\n" << rotmat << endl;
    cout << "Scales (x,y,z) = " << params.SubMatrix(7,9,1,1).t() << endl;
    cout << "Skews (xy,xz,yz) = " << params.SubMatrix(10,12,1,1).t() << endl;
    float avscale = (params(7) + params(8) + params(9))/3.0;
    cout << "Average scaling = " << avscale << endl << endl;
    Matrix m2 = sqrtaff(affmat);
    Matrix m0 = m2*affmat.i();
    cout << "Forward half transform =\n" << m2 << endl;
    cout << "Backward half transform =\n" << m0 << endl;

    Matrix scale(4,4), skew(4,4);
    Identity(scale);
    Identity(skew);
    scale(1,1) = params(7);
    scale(2,2) = params(8);
    scale(3,3) = params(9);
    skew(1,2) = params(10);
    skew(1,3) = params(11);
    skew(2,3) = params(12);
    Matrix ans;
    ans = affmat.i() * rotmat * skew * scale;
    //cout << endl << "Matrix check: mat^-1 * rotmat * skew * scale =\n";
    //cout << ans << endl;

    return 0;
  }
  catch(Exception exc) {
    cerr << exc.what() << endl;
    throw;
  }
  catch(...) {
    cerr << "Image error" << endl;
    throw;
  } 
  return(0);
}








