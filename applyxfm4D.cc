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
	 << " <output volume> <transformation matrix file/[dir]> [-4D [-singlematrix]]\n"; 
    return -1; 
  }

  string oname = argv[3], iname = argv[1], 
    transname = argv[4], refname = argv[2];
  bool fourd = false;
  if (argc>=6) {
    string option = argv[5];
    if (option == "-4D" )  fourd = true;
  }
  bool singlematrix = false;
  if (argc>=7) {
    string option = argv[6];
    if (option == "-singlematrix" )  singlematrix = true;
  }

  if (fourd) {
    volume4D<float> invol, outvol;
    volume<float> refvol, dummy;
    volumeinfo vinfo;
    read_volume4D(invol,iname,vinfo);
    invol.setextrapolationmethod(extraslice);

    refvol = invol[atoi(refname.c_str())];

    Matrix affmat(4,4);
    string matname;
    for (int m=invol.mint(); m<=invol.maxt(); m++) {

      if (singlematrix)
	{
	  read_matrix(affmat,transname,invol[0]);
	  cout << affmat << endl;
	}
      else
	{
	  matname = transname + "/MAT_0";
	  char nc='0';
	  int n = m;
	  matname += (nc+(n / 100));
	  n -= (n/100)*100;
	  matname += (nc+(n / 10));
	  n -= (n/10)*10;
	  matname += (nc+n);
	  cerr << matname << endl;
	  read_matrix(affmat,matname,invol[0]);
	}
      
      dummy = refvol;
      affine_transform(invol[m],dummy,affmat);
      outvol.addvolume(dummy);
    }
    save_volume4D(outvol,oname,vinfo);

  } else {
    volume<float> invol, outvol;
    volumeinfo vinfo;
    
    read_volume(invol,iname,vinfo);
    read_volume(outvol,refname);
    invol.setextrapolationmethod(extraslice);
    
    Matrix affmat(4,4);
    read_matrix(affmat,transname,invol);
    
    affine_transform(invol,outvol,affmat);
    save_volume(outvol,oname,vinfo);
  }

  return 0;
}




