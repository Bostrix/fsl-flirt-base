// Interpolation functions
//  Written by Mark Jenkinson  18/2/99

#if !defined(__costfns_h)
#define __costfns_h

#include "mjimage.h"
#include "miscimfns.h"

#ifndef NO_NAMESPACE
 using namespace MJIMAGE;
 using namespace MISCIMFNS;

 namespace COSTFNS {
#endif

  enum costfns { Woods, CorrRatio, MutualInfo, NormCorr, NormMI };

  float normcorr(const imagepair* ims, const Matrix& aff); 

  float woods_fn(const imagepair* ims, const Matrix& aff); 

  float corr_ratio(const imagepair* ims, const Matrix& aff); 

  float mutual_info(imagepair* ims, const Matrix& aff);

  float normalised_mutual_info(imagepair* ims, const Matrix& aff);


#ifndef NO_NAMESPACE
 }
#endif

#endif
