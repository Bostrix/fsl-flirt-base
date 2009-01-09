/*  applyxfm4D.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2000 University of Oxford  */

/*  CCOPYRIGHT  */

#include "newimage/newimageall.h"
#include "newimage/fmribmain.h"

using namespace NEWIMAGE;

// Globals - needed by fmrib_main

string oname, iname, transname, refname, matprefix="/MAT_0";
bool singlematrix, fourd;


//////////////////////////////////////////////////////////////////////

template <class T>
int fmrib_main(int argc, char* argv[])
{
  if (fourd) {
    // 4D mode
    volume4D<T> invol, outvol;
    volume<T> refvol, dummy;
    read_volume4D(invol,iname);
    for (int t=0; t<invol.tsize(); t++) {
      invol[t].setpadvalue(invol[t].backgroundval());
    }
    invol.setextrapolationmethod(extraslice);
    invol.setinterpolationmethod(sinc);
    invol.definesincinterpolation("b",7);

    // old form used a volume number
    //    refvol = invol[atoi(refname.c_str())];  

    read_volume(refvol,refname);

    Matrix affmat(4,4);
    string matname;
    if (singlematrix) { affmat = read_ascii_matrix(transname); }

    if (invol.maxt() - invol.mint() > 10000) {
      cerr << "WARNING:: More than 10000 volumes - only doing first 10000" << endl;
    }

    for (int m=invol.mint(); m<=Min(invol.maxt(),invol.mint()+10000); m++) {

      if (!singlematrix) {
	matname = transname + matprefix;
	char nc='0';
	int n = m;
	matname += (nc + (char) (n / 1000));
	n -= (n/1000)*1000;
	matname += (nc + (char) (n / 100));
	n -= (n/100)*100;
	matname += (nc + (char) (n / 10));
	n -= (n/10)*10;
	matname += (nc + (char) n);
	cout << matname << endl;
	affmat = read_ascii_matrix(matname);
      }
      
      dummy = refvol;
      affine_transform(invol[m],dummy,affmat);
      outvol.addvolume(dummy);
    }
    outvol.setDisplayMaximumMinimum(0,0);
    save_volume4D(outvol,oname);

  } else {
    // 3D mode
    volume<T> invol, outvol;
    
    read_volume(invol,iname);
    read_volume(outvol,refname);
    invol.setextrapolationmethod(extraslice);
    invol.setinterpolationmethod(sinc);
    invol.definesincinterpolation("r",9);

    Matrix affmat(4,4);
    affmat = read_ascii_matrix(transname);
    
    affine_transform(invol,outvol,affmat);
    outvol.setDisplayMaximumMinimum(0,0);
    save_volume(outvol,oname);
  }

  return 0;
}


int main(int argc,char *argv[])
{

  Tracer tr("main");
  if (argc<5) { 
    cerr << "Usage: " << argv[0] << " <input volume> <ref volume>"
	 << " <output volume> <transformation matrix file/[dir]> [-singlematrix/-fourdigit/-userprefix <prefix>]]\n"; 
    return -1; 
  }
  
  // NB: a hidden option (-3D) exists (must appear after singlematrix)


  // parse the command line
  oname = argv[3];
  iname = argv[1];
  transname = argv[4]; 
  refname = argv[2];
  singlematrix = false;
  if (argc>=6) {
    string option = argv[5];
    if (option == "-singlematrix" )  singlematrix = true;
    if (option == "-fourdigit" )  matprefix = "/MAT_";
    if ( (option == "-userprefix" ) && (argc>=7) ) {
      string uprefix = argv[6];
      matprefix = "/" + uprefix;
    }
  }
  fourd = true;
  if (argc>=7) {
    string option = argv[6];
    if (option == "-3D" )  fourd = false;
  }

  // call the templated main
  return call_fmrib_main(dtype(iname),argc,argv);

}



