/*  costfns.h

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */

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

  enum costfns { Woods, CorrRatio, MutualInfo, NormCorr, NormMI, LeastSq };

  float normcorr(const imagepair* ims, const Matrix& aff); 

  float normcorr_smoothed(const imagepair* ims, const Matrix& aff);

  float normcorr_fully_weighted(const imagepair* ims, const Matrix& aff,
				const volume& refweight, 
				const volume& testweight);
 
  float leastsquares(const imagepair* ims, const Matrix& aff);
 
  float leastsquares_smoothed(const imagepair* ims, const Matrix& aff);

  float leastsquares_fully_weighted(const imagepair* ims, const Matrix& aff, 
				    const volume& refweight, 
				    const volume& testweight);

  float woods_fn(const imagepair* ims, const Matrix& aff); 

  float woods_fn_smoothed(const imagepair* ims, const Matrix& aff); 

  float corr_ratio(const imagepair* ims, const Matrix& aff); 

  float corr_ratio_smoothed(const imagepair* ims, const Matrix& aff); 

  float corr_ratio_fully_weighted(const imagepair* ims, const Matrix& aff,
				  const volume& refweight, 
				  const volume& testweight);

  float mutual_info(imagepair* ims, const Matrix& aff);

  float mutual_info_smoothed(imagepair* ims, const Matrix& aff);

  float mutual_info_fully_weighted(imagepair* ims, const Matrix& aff,
				   const volume& refweight, 
				   const volume& testweight);

  float normalised_mutual_info(imagepair* ims, const Matrix& aff);

  float normalised_mutual_info_smoothed(imagepair* ims, const Matrix& aff);

  float normalised_mutual_info_fully_weighted(imagepair* ims, const Matrix& aff,
					      const volume& refweight, 
					      const volume& testweight);

#ifndef NO_NAMESPACE
 }
#endif

#endif







