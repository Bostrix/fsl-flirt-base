#include <string>
#include <iostream>
#include <fstream>
#include <unistd.h>
#define WANT_STREAM
#define WANT_MATH

#include "newmatap.h"
#include "newmatio.h"
#include "mjimage.h"
#include "miscmaths.h"
#include "miscimfns.h"
#include "generalio.h"

#ifndef NO_NAMESPACE
 using namespace MISCMATHS;
 using namespace MISCIMFNS;
 using namespace MJIMAGE;
 using namespace NEWMAT;
 using namespace GENERALIO;
#endif

////////////////////////////////////////////////////////////////////////////
// the real defaults are provided in the function parse_command_line

class globaloptions {
public:
  string testfname;
  string reffname;
  string outputmatascii;
  string outputmatmedx;
  string initmatfname;
  string xfm_type;
  int verbose;
  bool inverse;
public:
  globaloptions();
  ~globaloptions() {};
};

globaloptions globalopts;


globaloptions::globaloptions()
{
  // set up defaults
  reffname = "/usr/people/steve/reg/stroke_stpendle/average_305";

  testfname = "";
  outputmatascii = "";
  outputmatmedx = "";
  initmatfname = "";
  xfm_type = "u";
  inverse = false;
}


////////////////////////////////////////////////////////////////////////////

// Parsing functions for command line parameters

void print_usage(int argc, char *argv[])
{
  cout << "Usage: " << argv[0] << " [options] <input-matrix-filename>\n\n"
       << "  Available options are:\n"
       << "        -ref <refvol>                      (default is "
                                        << globalopts.reffname << ")\n"
       << "        -in <inputvol>                     (no default)\n"
       << "        -omat <matrix-filename>            (4x4 ascii format)\n"
       << "        -omedx <matrix-filename>           (MEDx format)\n"
    //       << "        -ominc <matrix-filename>           (MINC format)\n"
       << "        -xfmtype {u,a,m,g,s}              (Specify MEDx xfm format)\n"
       << "                      u   = UserTransformation (default)\n"
       << "                      a,m = AlignLinearReslice\n"
       << "                      g   = GenericReslice\n"
       << "                      s   = ShadowTransform\n"
       << "        -inverse\n"
       << "        -help\n";
}


void parse_command_line(int argc, char* argv[])
{
  if(argc<2){
    print_usage(argc,argv);
    exit(1);
  }


  int n=1;
  string arg;
  char first;

  while (n<argc) {
    arg=argv[n];
    if (arg.size()<1) { n++; continue; }
    first = arg[0];
    if (first!='-') {
      globalopts.initmatfname = arg;
      n++;
      continue;
    }
    
    // put options without arguments here
    if ( arg == "-help" ) {
      print_usage(argc,argv);
      exit(0);
    } else if ( arg == "-inverse" ) {
      globalopts.inverse = true;
      n++;
      continue;
    } else if ( arg == "-v" ) {
      globalopts.verbose = 5;
      n++;
      continue;
    }

    if (n+1>=argc) 
      { 
	cerr << "Lacking argument to option " << arg << endl;
	break; 
      }

    // put options with 1 argument here
    if ( arg == "-ref") {
      globalopts.reffname = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-in") {
      globalopts.testfname = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-omat") {
      globalopts.outputmatascii = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-omedx") {
      globalopts.outputmatmedx = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-xfmtype") {
      globalopts.xfm_type = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-verbose") {
      globalopts.verbose = atoi(argv[n+1]);
      n+=2;
      continue;
    }

    if (n+2>=argc) 
      { 
	cerr << "Lacking argument to option " << arg << endl;
	break; 
      }

  }  // while (n<argc)

  if (globalopts.initmatfname.size()<1) {
    cerr << "Input matrix filename not found\n\n";
    print_usage(argc,argv);
    exit(2);
  }
  if (globalopts.testfname.size()<1) {
    cerr << "WARNING:: Inputvol filename not found\n\n";
  }
}

////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[])
{
  try {

    parse_command_line(argc,argv);

    volume testvol, refvol;
    if (read_volume_hdr_only(testvol,globalopts.testfname)<0)  return -1;
    if (read_volume_hdr_only(refvol,globalopts.reffname)<0)  return -1;

    if (globalopts.verbose>3) {
      print_volume_info(refvol,"Reference Volume");
      cout << " origin = " << refvol.avw_origin.t() << endl << endl;
      print_volume_info(testvol,"Input Volume");
      cout << " origin = " << testvol.avw_origin.t() << endl;
    }

    Matrix affmat(4,4);
    ColumnVector params(12);
    if (read_matrix(affmat,globalopts.initmatfname,testvol,refvol)<0) 
      return -2;
    if (globalopts.inverse) {
      affmat = affmat.i();
    }
    if (globalopts.outputmatmedx.size() >= 1) {
      if (globalopts.inverse) {
	write_medx_matrix(affmat,globalopts.outputmatmedx,refvol,testvol,
			  globalopts.xfm_type,globalopts.testfname);
      } else {
	write_medx_matrix(affmat,globalopts.outputmatmedx,testvol,refvol,
			  globalopts.xfm_type,globalopts.reffname);
      }
    }
    if (globalopts.outputmatascii.size() >= 1) {
      write_ascii_matrix(affmat,globalopts.outputmatascii);
    }

    if (globalopts.verbose>0) {
      cout << affmat << endl;
    }

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




