// Interpolation functions
//  Written by Mark Jenkinson  18/2/99

#include <iostream>
#include <cassert>
#include <values.h>
#include "costfns.h"
#include "miscmaths.h"
#include "interpolation.h"

#ifndef NO_NAMESPACE
 using namespace MISCMATHS;
 using namespace INTERPOLATION;

 namespace COSTFNS {
#endif

   inline bool quick_in_plane_bounds(const float x, const float y, 
				     const float z, const float xb, 
				     const float yb, const float zb, 
				     const volume& vol) {
     return ( (y>=0.0) && (z>=0.0) && (y<yb) && (z<zb)
	      && (x>=vol.validxmin(y,z)) && (x<=vol.validxmax(y,z)) );
   }
   
   //--------------------------------------------------------------------//
   
   float corr_ratio(const volume& vref, const volume& vtest,
		    const Matrix& aff,
		    const float min, const float max,
		    const int no_bins)
    {
      // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1
      Matrix iaffbig = vtest.sampling_matrix().i() * aff.i() *
	                     vref.sampling_matrix();  
      Matrix iaff=iaffbig.SubMatrix(1,3,1,3);
      unsigned int xb1=vref.columns()-1, yb1=vref.rows()-1, zb1=vref.slices()-1;
      float  xb2 = ((float) vtest.columns())-1.0001,
	yb2=((float) vtest.rows())-1.0001, zb2=((float) vtest.slices())-1.0001;
      int io1, io2, io3;
      unsigned int x1, x2;

      float a11=iaff(1,1), a12=iaff(1,2), a13=iaff(1,3), a14=iaffbig(1,4),
	a21=iaff(2,1), a22=iaff(2,2), a23=iaff(2,3), a24=iaffbig(2,4),
	a31=iaff(3,1), a32=iaff(3,2), a33=iaff(3,3), a34=iaffbig(3,4), o1,o2,o3;

      float corr_ratio=0.0, var=0.0, totsumy=0.0, totsumy2=0.0;
      float *sumy, *sumy2;
      sumy = new float[no_bins+1];
      sumy2 = new float[no_bins+1];
      int *numy, *numx;
      numx = new int[no_bins+1];
      numy = new int[no_bins+1];
      int b=0, numtotx=0, numtoty=0;
      float a1=no_bins/(max-min), a0=-min*no_bins/(max-min), val;

      for (int i=0; i<=no_bins; i++) {
	numx[i]=0; numy[i]=0; sumy[i]=0.0;  sumy2[i]=0.0;
      }
  
  // cerr << "Bounds are " << xb << " " << yb << " " << zb << endl;
  // The matrix algebra below has been hand-optimized from
  //  [o1 o2 o3] = a * [x y z]  at each iteration
      for (unsigned int z=0; z<=zb1; z++) { 
	for (unsigned int y=0; y<=yb1; y++) { 
	  x1 = Max(vref.validxmin(y,z),(unsigned int) 0);
	  x2 = Min(vref.validxmax(y,z),xb1);
	  o1=x1*a11 + y*a12 + z*a13 + a14;  // x=x1
	  o2=x1*a21 + y*a22 + z*a23 + a24;  // x=x1
	  o3=x1*a31 + y*a32 + z*a33 + a34;  // x=x1
	  for (unsigned int x=x1; x<=x2; x++) {
	    if (quick_in_plane_bounds(o1,o2,o3,xb2,yb2,zb2,vtest)) {
	      io1=(int) o1;
	      io2=(int) o2;
	      io3=(int) o3;
	      val = q_tri_interpolation(vtest,o1,o2,o3,io1,io2,io3);
	      // do the cost function record keeping...
	      b=(int) (vref(x,y,z)*a1 + a0);
	      numx[b]++;
	      numy[b]++;
	      sumy[b]+=val;
	      sumy2[b]+=val*val;
	    }
	    o1 += a11;
	    o2 += a21;
	    o3 += a31;
	  }
	}
	//cerr << ".";
      }


      // now calculate the individual variances for each iso-set
      //  weighting them by the number of pixels from Image x that contribute
      for (b=0; b<=no_bins; b++) {
	numtotx += numx[b];
	if (numy[b]>2) {
	  numtoty += numy[b];
	  totsumy += sumy[b];
	  totsumy2 += sumy2[b];
	  // the following should be the variance of the bth iso-subset
	  var = (sumy2[b] - sumy[b]*sumy[b]/((float) numy[b]) ) /
	    ((float) (numy[b]-1));
	  // cerr << "Set #" << b << " has " << numx[b] << " elements and " 
	  //   << var << " variance" << endl;
	  corr_ratio += var * ((float) numx[b]);
	}
      }
      delete [] numx; delete [] numy; delete [] sumy; delete [] sumy2;

      // normalise the weighting of numx[]
      if (numtotx>0)  corr_ratio/=((float) numtotx);
      // calculate the total variance of Image y and then normalise by this
      if (numtoty>1)
	var = ( totsumy2 - totsumy*totsumy/((float) numtoty) ) /
	  ((float) (numtoty - 1));
      //cerr << "TOTALS are:" << endl 
      //   << " numerator variance is : " << corr_ratio << endl
      //   << " and denominator variance is: " << var << " from " << numtoty 
      //   << " valid elements" << endl;
      if (var>0.0)  corr_ratio/=var;
      // the above is actually 1 - correlation ratio, so correct this now
      if ( (numtotx<=0) || (numtoty<=1) || (var<=0.0) )
	return 0.0;   // the totally uncorrelated condition
      else
	return (1.0 - corr_ratio);

      // an alternative is to return 1.0/corr_ratio (=1/(1-correlation ratio))
      //  which may be better at rewarding gains near the best solution

      return 0;

    }


  ///////////////////////////////////////////////////////////////////////


  float woods_fn(const volume& vref, const volume& vtest,
		 const Matrix& aff,
		 const float min, const float max,
		 const int no_bins)
    {
      // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1
      Matrix iaffbig = vtest.sampling_matrix().i() * aff.i() *
	                     vref.sampling_matrix();  
      Matrix iaff=iaffbig.SubMatrix(1,3,1,3);
      unsigned int xb1=vref.columns()-1, yb1=vref.rows()-1, zb1=vref.slices()-1;
      float  xb2 = ((float) vtest.columns())-1.0001,
	yb2=((float) vtest.rows())-1.0001, zb2=((float) vtest.slices())-1.0001;
      int io1, io2, io3;
      unsigned int x1, x2;

      float a11=iaff(1,1), a12=iaff(1,2), a13=iaff(1,3), a14=iaffbig(1,4),
	a21=iaff(2,1), a22=iaff(2,2), a23=iaff(2,3), a24=iaffbig(2,4),
	a31=iaff(3,1), a32=iaff(3,2), a33=iaff(3,3), a34=iaffbig(3,4), o1,o2,o3;

      float *sum, *sum2;
      sum = new float[no_bins+1];
      sum2 = new float[no_bins+1];
      int *num;
      num = new int[no_bins+1];
      int b=0;
      float a1=no_bins/(max-min), a0=-min*no_bins/(max-min), val;

      for (int i=0; i<=no_bins; i++) {
	num[i]=0; sum[i]=0.0;  sum2[i]=0.0;
      }
  
  // cerr << "Bounds are " << xb << " " << yb << " " << zb << endl;
  // The matrix algebra below has been hand-optimized from
  //  [o1 o2 o3] = a * [x y z]  at each iteration
      for (unsigned int z=0; z<=zb1; z++) { 
	for (unsigned int y=0; y<=yb1; y++) { 
	  x1 = Max(vref.validxmin(y,z),(unsigned int) 0);
	  x2 = Min(vref.validxmax(y,z),xb1);
	  o1=x1*a11 + y*a12 + z*a13 + a14;  // x=x1
	  o2=x1*a21 + y*a22 + z*a23 + a24;  // x=x1
	  o3=x1*a31 + y*a32 + z*a33 + a34;  // x=x1
	  for (unsigned int x=x1; x<=x2; x++) {
	    if (quick_in_plane_bounds(o1,o2,o3,xb2,yb2,zb2,vtest)) {
	      io1=(int) o1;
	      io2=(int) o2;
	      io3=(int) o3;
	      val = q_tri_interpolation(vtest,o1,o2,o3,io1,io2,io3);
	      // do the cost function record keeping...
	      b=(int) (vref(x,y,z)*a1 + a0);
	      num[b]++;
	      sum[b]+=val;  // val = vtest(x',y',z')
	      sum2[b]+=val*val;
	    }
	    o1 += a11;
	    o2 += a21;
	    o3 += a31;
	  }
	}
	//cerr << ".";
      }


      // now calculate  W = sum_j (n_j/N)*(sigma_j / mu_j)
      //  where n_j = num[j], N = sum_j n_j, mu_j = sum[j]/num[j]
      //        sigma_j^2 = sum_j 1/(n_j - 1) * (sum2[j] - mu_j^2 * n_j)
      float woods=0.0, stdev=0.0, var=0.0;
      int numtot=0;
      for (b=0; b<=no_bins; b++) {
	if (num[b]>2) {
	  numtot += num[b];
	  // the following should be the variance of the bth subset
	  var = (sum2[b] - sum[b]*sum[b]/((float) num[b]) ) /
	    ((float) (num[b]-1));
	  if (var>0.0)
	    stdev = sqrt(var);
	  else
	    stdev = 0.0;
	  if (sum[b]>0)
	    woods += Sqr((float) num[b])*stdev/sum[b];
	  else
	    woods += Sqr((float) num[b])*stdev;
	}
      }
      delete [] num; delete [] sum; delete [] sum2;
      if (numtot>0) {
	woods/=((float) numtot);
	return woods;
      } else {
	return MAXFLOAT;
      }
    }


  ///////////////////////////////////////////////////////////////////////


  float normcorr(const volume& vref, const volume& vtest,
		 const Matrix& aff)
    {
      // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1
      Matrix iaffbig = vtest.sampling_matrix().i() * aff.i() *
	                     vref.sampling_matrix();  
      Matrix iaff=iaffbig.SubMatrix(1,3,1,3);
      unsigned int xb1=vref.columns()-1, yb1=vref.rows()-1, zb1=vref.slices()-1;
      float  xb2 = ((float) vtest.columns())-1.0001,
	yb2=((float) vtest.rows())-1.0001, zb2=((float) vtest.slices())-1.0001;
      int io1, io2, io3;
      unsigned int x1, x2;

      float a11=iaff(1,1), a12=iaff(1,2), a13=iaff(1,3), a14=iaffbig(1,4),
	a21=iaff(2,1), a22=iaff(2,2), a23=iaff(2,3), a24=iaffbig(2,4),
	a31=iaff(3,1), a32=iaff(3,2), a33=iaff(3,3), a34=iaffbig(3,4), o1,o2,o3;

      float corr=0.0;
      float sumx=0.0, sumy=0.0, sumx2=0.0, sumy2=0.0, sumxy=0.0;
      float varx=0.0, vary=0.0, varxy=0.0;
      float valx=0.0, valy=0.0, val=0.0;
      int num=0;

      // cerr << "Bounds are " << xb << " " << yb << " " << zb << endl;
      // The matrix algebra below has been hand-optimized from
      //  [o1 o2 o3] = a * [x y z]  at each iteration
      for (unsigned int z=0; z<=zb1; z++) { 
	for (unsigned int y=0; y<=yb1; y++) { 
	  x1 = Max(vref.validxmin(y,z),(unsigned int) 0);
	  x2 = Min(vref.validxmax(y,z),xb1);
	  o1=x1*a11 + y*a12 + z*a13 + a14;  // x=x1
	  o2=x1*a21 + y*a22 + z*a23 + a24;  // x=x1
	  o3=x1*a31 + y*a32 + z*a33 + a34;  // x=x1
	  for (unsigned int x=x1; x<=x2; x++) {
	    if (quick_in_plane_bounds(o1,o2,o3,xb2,yb2,zb2,vtest)) {
	      io1=(int) o1;
	      io2=(int) o2;
	      io3=(int) o3;
	      val = q_tri_interpolation(vtest,o1,o2,o3,io1,io2,io3);
	      // do the cost function record keeping...
	      num++;
	      valx = vref(x,y,z);
	      valy = val;  // vtest(x',y',z')
	      sumx += valx;
	      sumx2 += valx*valx;
	      sumy += valy;
	      sumy2 += valy*valy;
	      sumxy += valx*valy;
	    }
	    o1 += a11;
	    o2 += a21;
	    o3 += a31;
	  }
	}
	//cerr << ".";
      }

  
      corr = 0.0;  // uncorrelated (worst) case
      if (num>2) {
	varxy = sumxy/((float) num-1) - (sumx*sumy)/((float) num*num);
	varx = sumx2/((float) num-1) - (sumx*sumx)/((float) num*num);
	vary = sumy2/((float) num-1) - (sumy*sumy)/((float) num*num);
	if ((varx>0.0) && (vary>0.0)) {
	  corr = varxy/sqrt(varx)/sqrt(vary);
	} 
      }
      return corr;
    }

  
  ///////////////////////////////////////////////////////////////////////

  void calc_entropy(const volume& vref, const volume& vtest,
		    const Matrix& aff,
		    const float minref, const float maxref,
		    const float mintest, const float maxtest,
		    const int no_bins, const ColumnVector& plnp, 
		    int *jointhist, int *marghist1, int *marghist2,
		    float& jointentropy, float& margentropy1,
		    float& margentropy2)
    {
      // actually the joint entropy between the two images is returned
      //  since minimising the joint entropy is equivalent to maximising
      //  the mutual information
      // the last parameter, plnp, is a vector containing values of -p*log(p)
      //  which makes the calculation significantly more efficient

      // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1
      Matrix iaffbig = vtest.sampling_matrix().i() * aff.i() *
	                     vref.sampling_matrix();  
      Matrix iaff=iaffbig.SubMatrix(1,3,1,3);
      unsigned int xb1=vref.columns()-1, yb1=vref.rows()-1, zb1=vref.slices()-1;
      float  xb2 = ((float) vtest.columns())-1.0001,
	yb2=((float) vtest.rows())-1.0001, zb2=((float) vtest.slices())-1.0001;
      int io1, io2, io3;
      unsigned int x1, x2;

      float a11=iaff(1,1), a12=iaff(1,2), a13=iaff(1,3), a14=iaffbig(1,4),
	a21=iaff(2,1), a22=iaff(2,2), a23=iaff(2,3), a24=iaffbig(2,4),
	a31=iaff(3,1), a32=iaff(3,2), a33=iaff(3,3), a34=iaffbig(3,4), o1,o2,o3;

      for (long int i=0; i<=(no_bins*no_bins+2*no_bins-1); i++) {
	  jointhist[i]=0;
      }
      for (int i=0; i<=no_bins; i++) {
	  marghist1[i]=0;
	  marghist2[i]=0;
      }

      long int a,b;
      float a1=no_bins/(maxref-minref), a0=-minref*no_bins/(maxref-minref);
      float b1=no_bins/(maxtest-mintest), b0=-mintest*no_bins/(maxtest-mintest);
      float val=0.0;
 
      // cerr << "Bounds are " << xb2 << " " << yb2 << " " << zb2 << endl;
      // The matrix algebra below has been hand-optimized from
      //  [o1 o2 o3] = a * [x y z]  at each iteration
      for (unsigned int z=0; z<=zb1; z++) { 
	for (unsigned int y=0; y<=yb1; y++) { 
	  x1 = Max(vref.validxmin(y,z),(unsigned int) 0);
	  x2 = Min(vref.validxmax(y,z),xb1);
	  o1=x1*a11 + y*a12 + z*a13 + a14;  // x=x1
	  o2=x1*a21 + y*a22 + z*a23 + a24;  // x=x1
	  o3=x1*a31 + y*a32 + z*a33 + a34;  // x=x1
	  for (unsigned int x=x1; x<=x2; x++) {
	    if (quick_in_plane_bounds(o1,o2,o3,xb2,yb2,zb2,vtest)) {
	      io1=(int) o1;
	      io2=(int) o2;
	      io3=(int) o3;
	      val = q_tri_interpolation(vtest,o1,o2,o3,io1,io2,io3);
	      // do the cost function record keeping...
	      a=(long int) (vref(x,y,z)*a1 + a0) + 1;
	      b=(long int) (val*b1 + b0) + 1;  // val = vtest(x',y',z')
	      (jointhist[(a-1)*(no_bins+1) + b-1])++;
	      (marghist1[a-1])++;
	      (marghist2[b-1])++;
	    }
	    o1 += a11;
	    o2 += a21;
	    o3 += a31;
	  }
	}
      }


      // note that the plnp values indexed by integers such that:
      //   plnp(n) = n/N * log(n/N)
      float p=0.0;
      int n=0, psize=plnp.Nrows();
      long int nvoxels = (long int) (vref.rows() * vref.columns() * 
				     vref.slices());
      for (long int i=0; i<(no_bins*no_bins+2*no_bins); i++) {
	n = jointhist[i];
	if (n>0) {
	  if (n<=psize)
	    jointentropy+=plnp(n);
	  else {
	    p = ((float) n) / ((float) nvoxels);
	    jointentropy+= -((float) p)*(log((float) p));
	  }
	}
      }
      for (int i=0; i<=no_bins; i++) {
	n = marghist1[i];
	if (n>0) {
	  if (n<=psize)
	    margentropy1+=plnp(n);
	  else {
	    p = ((float) n) / ((float) nvoxels);
	    margentropy1+= -((float) p)*(log((float) p));
	  }
	}
      }
      long int noverlap=0;
      for (int i=0; i<=no_bins; i++) {
	n = marghist2[i];
	if (n>0) {
	  noverlap += n;
	  if (n<=psize)
	    margentropy2+=plnp(n);
	  else {
	    //cerr << ":";
	    p = ((float) n) / ((float) nvoxels);
	    margentropy2+= -((float) p)*(log((float) p));
	  }
	}
      }

      // correct for difference in total histogram size
      //  that is: noverlap vs nvoxels
      // H_1 = N_0/N_1 * H_0 + log(N_1/N_0) + (1 - N_0/N_1)*log(N_0)
      if (noverlap > 0) {
	float nratio = ((float) nvoxels) / ((float) noverlap);
	jointentropy = nratio * jointentropy + log(1.0/nratio)
	                + (1.0 - nratio)*log(nvoxels);
	margentropy1 = nratio * margentropy1 + log(1.0/nratio)
	                + (1.0 - nratio)*log(nvoxels);
	margentropy2 = nratio * margentropy2 + log(1.0/nratio)
	                + (1.0 - nratio)*log(nvoxels);
      } else {
	// Put in maximum entropy values as base cases = BAD registration
	jointentropy = 2.0*log(no_bins);
	margentropy1 = log(no_bins);
	margentropy2 = log(no_bins);
      }
      return;
    }



  float mutual_info(const volume& vref, const volume& vtest,
		    const Matrix& aff,
		    const float minref, const float maxref,
		    const float mintest, const float maxtest,
		    const int no_bins, const ColumnVector& plnp, 
		    int *jointhist, int *marghist1, int *marghist2)
    {
      float jointentropy=0.0, margentropy1=0.0, margentropy2=0.0;
      calc_entropy(vref,vtest,aff,minref,maxref,mintest,maxtest,no_bins,
		   plnp,jointhist,marghist1,marghist2,
		   jointentropy,margentropy1,margentropy2);
      float mutualinformation = margentropy1 + margentropy2 - jointentropy;
      return mutualinformation;
    }



  float normalised_mutual_info(const volume& vref, const volume& vtest,
			       const Matrix& aff,
			       const float minref, const float maxref,
			       const float mintest, const float maxtest,
			       const int no_bins, const ColumnVector& plnp, 
			       int *jointhist, int *marghist1, int *marghist2)
    {
      float jointentropy=0.0, margentropy1=0.0, margentropy2=0.0;
      calc_entropy(vref,vtest,aff,minref,maxref,mintest,maxtest,no_bins,
		   plnp,jointhist,marghist1,marghist2,
		   jointentropy,margentropy1,margentropy2);
      float normmi;
      if (fabs(jointentropy)<1e-9) {
	normmi = 0.0;  // BAD registration result
      } else {
	normmi = (margentropy1 + margentropy2)/jointentropy;
      }
      return normmi;
    }


  ///////////////////////////////////////////////////////////////////////

  // Supporting interfaces


  float normcorr(const imagepair* ims, const Matrix& aff) 
    {
      return normcorr(ims->refvol,ims->testvol,aff);
    }


  float woods_fn(const imagepair* ims, const Matrix& aff) 
    {
      return woods_fn(ims->refvol,ims->testvol,aff,
			ims->refmin,ims->refmax,ims->no_bins);
    }


  float corr_ratio(const imagepair* ims, const Matrix& aff) 
    {
      return corr_ratio(ims->refvol,ims->testvol,aff,
			ims->refmin,ims->refmax,ims->no_bins);
    }
  

  float mutual_info(imagepair* ims, const Matrix& aff)
    {
      return mutual_info(ims->refvol,ims->testvol,aff,
			 ims->refmin,ims->refmax,ims->testmin,ims->testmax,
			 ims->no_bins,ims->plnp,ims->jointhist,
			 ims->marghist1,ims->marghist2);
    }


  float normalised_mutual_info(imagepair* ims, const Matrix& aff)
    {
      return mutual_info(ims->refvol,ims->testvol,aff,
			 ims->refmin,ims->refmax,ims->testmin,ims->testmax,
			 ims->no_bins,ims->plnp,ims->jointhist,
			 ims->marghist1,ims->marghist2);
    }


  ///////////////////////////////////////////////////////////////////////////

#ifndef NO_NAMESPACE
 }
#endif



