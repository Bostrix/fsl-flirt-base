/*  applyxfm.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2000 University of Oxford  */

/*  CCOPYRIGHT  */

#include "newimageall.h"

using namespace NEWIMAGE;

int main(int argc,char *argv[])
{

  Tracer tr("main");
  if (argc<5) { 
    cerr << "Usage: " << argv[0] << " <input volume> <ref volume>"
	 << " <output volume> <transformation matrix>\n"; 
    return -1; 
  }

  string oname = argv[3], iname = argv[1], 
    transname = argv[4], refname = argv[2];

  volume<float> invol, outvol;

  read_volume(invol,iname);
  read_volume(outvol,refname);
  invol.setextrapolationmethod(extraslice);

  Matrix affmat(4,4);
  read_matrix(affmat,transname,invol);

  affine_transform(invol,outvol,affmat);
  save_volume(outvol,oname);

  return 0;
}




