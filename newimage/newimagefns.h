/*  newimagefns.h

    Mark Jenkinson et al, FMRIB Image Analysis Group

    Copyright (C) 2000-2006 University of Oxford  */

/*  CCOPYRIGHT  */

// General image processing functions

#if !defined(__newimagefns_h)
#define __newimagefns_h

#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <queue>
#include "armawrap/newmat.h"
#include "miscmaths/miscmaths.h"
#include "newimage.h"
#include "complexvolume.h"
#include "imfft.h"


#ifndef MAX
#define MAX(a,b) (((a)>(b))?(a):(b))
#endif

#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif

namespace NEWIMAGE {


  // The following lines are ignored by the current SGI compiler
  //  (version egcs-2.91.57)
  // A temporary fix of including the std:: in front of all abs() etc
  //  has been done for now
  using std::abs;
  using std::sqrt;

  ///////////////////////////////////////////////////////////////////////////

  // BASIC IMAGE SUPPORT FUNCTIONS

  // Complex volume support
  volume<float> abs(const volume<float>& realvol,
		    const volume<float>& imagvol);

  volume<float> phase(const volume<float>& realvol,
		      const volume<float>& imagvol);

  volume<float> real(const volume<float>& absvol,
		     const volume<float>& phasevol);

  volume<float> imag(const volume<float>& absvol,
		     const volume<float>& phasevol);

  // Basic Arithmetic Operations
  template <class T>
  volume<T> abs(const volume<T>& vol);

  template <class T, class S>
  volume<T> divide(const volume<T>& numervol, const volume<T>& denomvol,
		   const volume<S>& mask);

  template <class T, class S>
  volume<T> mask_volume( const volume<T>& invol,const volume<S>& mask );

  // General Mathematical Operations

  template <class T>
  void clamp(volume<T>& vol, T minval, T maxval);


  template <class T>
  volume<T> binarise(const volume<T>& vol, T lowerth, T upperth, threshtype tt=inclusive, bool invert=false);
  template <class T>
  volume<T> binarise(const volume<T>& vol, T thresh, bool invert=false);


  template <class T>
  volume<T> threshold(const volume<T>& vol, T lowerth, T upperth, threshtype tt=inclusive);
  template <class T>
  volume<T> threshold(const volume<T>& vol, T thresh);


  template <class T>
  volume<T> min(const volume<T>& v1, const volume<T>& v2);
  template <class T>
  volume<T> max(const volume<T>& v1, const volume<T>& v2);

  volume<float> sqrt(const volume<char>& vol);
  volume<float> sqrt(const volume<short>& vol);
  volume<float> sqrt(const volume<int>& vol);
  volume<float> sqrt(const volume<float>& vol);
  volume<double> sqrt(const volume<double>& vol);
  template <class T>
  volume<float> sqrt_float(const volume<T>& vol);

  template <class T>
  volume<float> meanvol(const volume4D<T>& vol4);
  template <class T>
  volume<float> stddevvol(const volume4D<T>& vol4);
  template <class T>
  volume<float> variancevol(const volume4D<T>& vol4);
  template <class T>
  volume<float> sumvol(const volume4D<T>& vol4);
  template <class T>
  volume<float> sumsquaresvol(const volume4D<T>& vol4);
  template <class T>
  volume<float> dotproductvol(const volume4D<T>& vol4,
                              const NEWMAT::ColumnVector& vec);

  template <class T>
  void pad(const volume<T>& vol, volume<T>& paddedvol);
  template <class T>
  void pad(const volume<T>& vol, volume<T>& paddedvol,
	     int offsetx, int offsety, int offsetz);

  // Considers each volume as a vector and returns the dotproduct
  template <class T>
  double dotproduct(const volume<T>&  vol1,
                    const volume<T>&  vol2);
  template <class T, class S>
  double dotproduct(const volume<T>&  vol1,
                    const volume<T>&  vol2,
                    const volume<S>&  mask);
  template <class T, class S>
  double dotproduct(const volume<T>& vol1,
                    const volume<T>& vol2,
                    const volume<S>  *mask);

  // This is a bit of a special-needs funtion. If we consider vol1 and vol2
  // as column vectors v1 and v2 then powerdotproduct returns (v1.^n1)' * (v2.^n2)
  template <class T>
  double powerdotproduct(const volume<T>&  vol1,
                         unsigned int      n1,
                         const volume<T>&  vol2,
                         unsigned int      n2);
  template <class T, class S>
  double powerdotproduct(const volume<T>&  vol1,
                         unsigned int      n1,
                         const volume<T>&  vol2,
                         unsigned int      n2,
                         const volume<S>&  mask);
  template <class T, class S>
  double powerdotproduct(const volume<T>&  vol1,
                         unsigned int      n1,
                         const volume<T>&  vol2,
                         unsigned int      n2,
                         const volume<S>   *mask);

  // These are global functions that duplicate some of the member functions.
  // The reason is that we want to be able to double template them in a
  // convenient manner (e.g. use a <char> mask for <float> data).
  template <class T>
  double mean(const volume<T>&  vol);
  template <class T, class S>
  double mean(const volume<T>&  vol,
              const volume<S>&  mask);
  template <class T, class S>
  double mean(const volume<T>& vol,
              const volume<S> *mask);

  template <class T>
  double sum(const volume<T>&  vol);
  template <class T, class S>
  double sum(const volume<T>&  vol,
             const volume<S>&  mask);
  template <class T, class S>
  double sum(const volume<T>& vol,
             const volume<S> *mask);


  ///////////////////////////////////////////////////////////////////////////

  // IMAGE PROCESSING ROUTINES

  template <class T>
  void affine_transform(const volume<T>& vin, volume<T>& vout,
			     const NEWMAT::Matrix& aff, float paddingsize=0.0);
  template <class T>
  void affine_transform(const volume<T>& vin, volume<T>& vout,
			const NEWMAT::Matrix& aff, float paddingsize, bool set_backgnd);

  template <class T>
  void affine_transform(const volume<T>& vin, volume<T>& vout,
			const NEWMAT::Matrix& aff, interpolation interptype,
			float paddingsize=0.0);
  template <class T>
  volume<T> affine_transform_mask(const volume<T>& vin, const volume<T>& vout,
				  const NEWMAT::Matrix& aff, float padding=0.0);

  template <class T>
  void get_axis_orientations(const volume<T>& inp1,
			     int& icode, int& jcode, int& kcode);


  // the following convolve function do not attempt to normalise the kernel
  template <class T, class S>
  volume<T> convolve(const volume<T>& source, const volume<S>& kernel);
  template <class T, class S>
  volume<T> convolve_separable(const volume<T>& source,
			       const NEWMAT::ColumnVector& kernelx,
			       const NEWMAT::ColumnVector& kernely,
			       const NEWMAT::ColumnVector& kernelz);

  // the following convolve functions take in a mask and also renormalise
  //  the result according to the overlap of kernel and mask at each point
  // NB: these functions should NOT be used with zero-sum kernels (eg.Laplacian)
  template <class T, class S, class M>
  volume<T> convolve(const volume<T>& source, const volume<S>& kernel,
		     const volume<M>& mask, bool ignoremask=false, bool renormalise=true);
  template <class T, class M>
  volume<T> convolve_separable(const volume<T>& source,
			       const NEWMAT::ColumnVector& kernelx,
			       const NEWMAT::ColumnVector& kernely,
			       const NEWMAT::ColumnVector& kernelz,
			       const volume<M>& mask, bool ignoremask=false, bool renormalise=true);

  //This implements Professor Smith's SUSAN convolve algorithm, note the number of optional parameters
  template <class T, class S>
    volume<T> susan_convolve(volume<T> source, const volume<S>& kernel, const float sigmabsq, const bool use_median, int num_usan,volume<T>* usan_area = new volume<T>(1,1,1),volume<T> usan_vol1=volume<T>(1,1,1),const float sigmab1sq=0,volume<T> usan_vol2 = volume<T>(1,1,1),const float sigmab2sq=0);

  template <class T, class S>
    volume4D<T> generic_convolve(const volume4D<T>& source, const volume<S>& kernel, bool seperable=false, bool renormalise=true);
  template <class T, class S>
    volume<T> efficient_convolve(const volume<T>& vin, const volume<S>& vker);
  template <class T, class S>
  int insertpart(volume<T>& v1, const volume<S>& v2);
  template <class T, class S, class U>
  volume<S> extractpart(const volume<T>& v1, const volume<S>& v2, const volume<U>& kernel) ;
   float fsllog2(float x);


   template <class T>
     volume4D<T> bandpass_temporal_filter(volume4D<T>& source,double hp_sigma, double lp_sigma);


  template <class T, class S>
  volume<T> morphfilter(const volume<T>& source, const volume<S>& kernel,
			const std::string& filtertype);

  template <class T>
  volume<T> dilall(const volume<T>& im, volume<T>& mask);
  template <class T>
  volume<T> fill_holes(const volume<T>& im, int connectivity=6);

  template <class T>
  volume<T> isotropic_resample(const volume<T>& aniso, float scale);

  template <class T>
  volume<T> subsample_by_2(const volume<T>& refvol, bool centred=true);
  template <class T>
  int upsample_by_2(volume<T>& highresvol, const volume<T>& lowresvol,
		    bool centred=true);
  template <class T>
  volume<T> upsample_by_2(const volume<T>& lowresvol, bool centred=true);

  // for all blur functions the size of blurring is in mm
  void make_blur_mask(NEWMAT::ColumnVector& bmask, const float final_vox_dim,
		     const float init_vox_dim);
  template <class T>
  volume<T> blur(const volume<T>& source, const NEWMAT::ColumnVector& resel_size);
  template <class T>
  volume<T> blur(const volume<T>& source, float iso_resel_size);

    /*template <class T>
  volume<T> smooth(const volume<T>& source, float sigma_mm);
  template <class T>
  volume<T> smooth2D(const volume<T>& source, float sigma_mm, int nulldir=3);*/


  NEWMAT::ColumnVector gaussian_kernel1D(float sigma, int radius);
  volume<float> gaussian_kernel2D(float sigma, int radius);
  volume<float> gaussian_kernel3D(float sigma, int radius);
  volume<float> gaussian_kernel3D(float sigma,float xdim,float ydim,float zdim,float cutoff=4.0);
  volume<float> spherical_kernel(float radius, float xdim, float ydim, float zdim);
  volume<float> box_kernel(float length,float xdim,float ydim,float zdim);  //mm dimensions
  volume<float> box_kernel(int x,int y, int z);                        //voxel dimensions

  void make_grad_masks(volume<float>& maskx, volume<float>& masky,
		       volume<float>& maskz);

  template <class T>
  volume<float> gradient(const volume<T>& source);  // in voxel coords

  template <class T>
  void gradient(const volume<T>& source,volume<float>& grad);  // in voxel coords

  // separate left and right gradients (changes at voxel mid-point)
  template <class T>
  volume4D<float> lrxgrad(const volume<float>& im, const volume<T>& mask);
  template <class T>
  volume4D<float> lrygrad(const volume<float>& im, const volume<T>& mask);
  template <class T>
  volume4D<float> lrzgrad(const volume<float>& im, const volume<T>& mask);


  template <class T>
  volume<T> log_edge_detect(const volume<T>& source,
			    float sigma1, float sigma2,
			    int mode=0); //mode==0: 3D, mode==1: 2D x, mode==2: 2D y, mode==3: 2D z
  template <class T>
  volume<T> fixed_edge_detect(const volume<T>& source, float threshold,
			      bool twodimensional=false);
  template <class T>
  volume4D<T> edge_strengthen(const volume4D<T>& source);



  template <class T>
    volume<int> connected_components(const volume<T>& vol,  NEWMAT::ColumnVector& clustersize, int numconnected=26);

  template <class T>
    volume<int> connected_components(const volume<T>& vol,
                                   const volume<T>& mask,
                                   bool (*binaryrelation)(T , T), NEWMAT::ColumnVector& clustersize);

  template <class T>
  volume<int> connected_components(const volume<T>& vol,
				   int numconnected=26);
  template <class T>
  volume<int> connected_components(const volume<T>& vol,
                                   const volume<T>& mask,
                                   bool (*binaryrelation)(T , T));

  template <class T>
  volume<float> distancemap(const volume<T>& binaryvol);
  template <class T>
  volume<float> distancemap(const volume<T>& binaryvol, const volume<T>& mask);

  template <class T>
  volume4D<float> sparseinterpolate(const volume4D<T>& sparsesamps,
				    const volume<T>& mask,
				    const std::string& interpmethod="general");
  // can have "general" or "nearestneighbour" (or "nn") for interpmethod


 template <class T>
 NEWMAT::Matrix NewimageVox2NewimageVoxMatrix(const NEWMAT::Matrix& flirt_in2ref,
		      const volume<T>& invol, const volume<T>& refvol);



  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  // DEBUG
  template <class T>
  void print_volume_info(const volume<T>& source, const std::string& name, std::ostream& output=std::cout)
  {
    output << name << ":: Size = (" << source.xsize() << ","
	 << source.ysize() << "," << source.zsize() << "," << source.tsize() << ")" << std::endl;
    output << name << ":: Dims = (" << source.xdim() << ","
	 << source.ydim() << "," << source.zdim() << "," << source.tdim() << ")" << std::endl;
    output << name << ":: ROI Size = (" << source.maxx() - source.minx() + 1
	 << "," << source.maxy() - source.miny() + 1
	 << "," << source.maxz() - source.minz() + 1
	 << "," << source.maxt() - source.mint() + 1 << ")" << std::endl;
    output << name << ":: Minimum and maximum intensities are: "
         << source.min() << " and " << source.max() << std::endl;
  }

  ///////////////////////////////////////////////////////////////////////////
  // BASIC IMAGE SUPPORT FUNCTIONS
  ///////////////////////////////////////////////////////////////////////////

  template <class T>
  void clamp(volume<T>& vol, T minval, T maxval)
  {
      if (maxval < minval) return;
      for ( typename volume<T>::nonsafe_fast_iterator it=vol.nsfbegin(), itEnd=vol.nsfend(); it < itEnd; it++ ) {
	if ( *it > maxval )
	  *it = maxval;
	else if ( *it < minval )
	  *it = minval;
      }
  }

  template <class T>
  volume<T> abs(const volume<T>& vol)
  {
      volume<T> newvol(vol);
      for ( typename volume<T>::nonsafe_fast_iterator it=newvol.nsfbegin(), itEnd=newvol.nsfend(); it < itEnd; it++ ) {
	*it = (T)fabs( (double)*it );
      }
      return newvol;
  }

  template <class T>
    volume<T> binarise(const volume<T>& vol, T lowerth, T upperth, threshtype tt, bool invert)
    {
      volume<T> newvol(vol);
      newvol.binarise(lowerth,upperth,tt,invert);
      return newvol;
    }

  template <class T>
    volume<T> binarise(const volume<T>& vol, T lowthresh, bool invert)
    {
      return binarise(vol,lowthresh,vol.max(),inclusive, invert);
    }

  template <class T>
  volume<T> threshold(const volume<T>& vol, T lowerth, T upperth, threshtype tt)
    {
      volume<T> newvol(vol);
      newvol.threshold(lowerth,upperth,tt);
      return newvol;
    }

  template <class T>
  volume<T> threshold(const volume<T>& vol, T thresh)
    {
      return threshold(vol,thresh,vol.max(),inclusive);
    }

  template <class T>
  volume<T> min(const volume<T>& v1, const volume<T>& v2)
  {
      if (!samesize(v1,v2)) {
	imthrow("Must use volumes of same size in min(v1,v2)",3);
      }
      volume<T> newvol(v1, false);
      typename volume<T>::nonsafe_fast_iterator it=newvol.nsfbegin();
      for ( typename volume<T>::fast_const_iterator v1it=v1.fbegin(), v1itEnd=v1.fend(), v2it=v2.fbegin(); v1it < v1itEnd; it++, v1it++, v2it++ )
        *it=MISCMATHS::Min( *v1it, *v2it );
      return newvol;
  }

  template <class T>
  volume<T> max(const volume<T>& v1, const volume<T>& v2)
    {
      if (!samesize(v1,v2)) {
	imthrow("Must use volumes of same size in max(v1,v2)",3);
      }
      volume<T> newvol;
      copyconvert(v1, newvol, false);
      typename volume<T>::nonsafe_fast_iterator it=newvol.nsfbegin();
      for ( typename volume<T>::fast_const_iterator v1it=v1.fbegin(), v1itEnd=v1.fend(), v2it=v2.fbegin(); v1it < v1itEnd; it++, v1it++, v2it++ )
	*it=MISCMATHS::Max( *v1it, *v2it );
      return newvol;
    }

  template <class T>
  volume<float> sqrt_float(const volume<T>& vol)
  {
    volume<float> retvol;
    copyconvert(vol, retvol, false);
    typename volume<float>::nonsafe_fast_iterator it=retvol.nsfbegin();
    for ( typename volume<T>::fast_const_iterator vit=vol.fbegin(), vend=vol.fend(); vit < vend; it++, vit++ ) {
      if ( *vit > 0 )
	*it=sqrt( (double) *vit);
      else
	*it=0;
    }
    return retvol;
  }

  template <class T>
    volume<float> sumvol(const volume<T>& vol4) //4D only
  {
    if (vol4.tsize()<1) { volume<float> newvol; return newvol; }
    volume<float> SumVol, dummy;
    copyconvert(vol4.constSubVolume(0),SumVol);
    for (int ctr= 1; ctr < vol4.tsize(); ctr++) {
      copyconvert(vol4.constSubVolume(ctr),dummy);
      SumVol += dummy;
    }
    return SumVol;
  }


  template <class T>
    volume<float> meanvol(const volume<T>& vol4)//4D only
  {
    if (vol4.tsize()<1) { volume<float> newvol; return newvol; }
    volume<float> MeanVol = sumvol(vol4) / (float)vol4.tsize();
    return MeanVol;
  }


  template <class T>
    volume<float> sumsquaresvol(const volume<T>& vol4) //4D only
  {
    if (vol4.tsize()<1) { volume<float> newvol; return newvol; }
    volume<float> SumSq, dummy;
    copyconvert( vol4.constSubVolume(0) * vol4.constSubVolum(0) , SumSq);
    for (int ctr=1; ctr < vol4.tsize(); ctr++) {
      copyconvert( vol4.constSubVolume(ctr) * vol4.constSubVolume(ctr), dummy);
      SumSq += dummy;
    }
    return SumSq;
  }

  //This rewrite of variancevol gives the same output as doing the sumsquaresvol etc in double precision internally
  template <class T>
    volume<float> variancevol(const volume<T>& vol4) //4D only
  {
     volume<float> variance;
     if (vol4.tsize()<1)
       return variance;
     volume<float> Mean = meanvol(vol4);
     variance = Mean*0.0f;

     for (int z=0; z<vol4.zsize(); z++)
       for (int y=0; y<vol4.ysize(); y++)
	 for (int x=0; x<vol4.xsize(); x++) {
	   double total(0);
	   for (int t=0; t<vol4.tsize(); t++)
	     total+=pow(vol4(x,y,z,t)-Mean(x,y,z),2.0);
	   variance(x,y,z)=(float)total;
	 }
     variance /= (float) (vol4.tsize()-1.0);
     return variance;
  }

  template <class T>
    volume<float> stddevvol(const volume<T>& vol4) //4D only
  {
    if (vol4.tsize()<1) { volume<float> newvol; return newvol; }
    volume<float> StdVol = variancevol(vol4);
    for (int z=0; z<StdVol.zsize(); z++)
      for (int y=0; y<StdVol.ysize(); y++)
	for (int x=0; x<StdVol.xsize(); x++)
	  StdVol(x,y,z) = sqrt(StdVol(x,y,z));

    return StdVol;
  }

  template <class T>
    volume<float> dotproductvol(const volume4D<T>& vol4, //4D only
			      const NEWMAT::ColumnVector& vec)
  {
    if (vol4.mint()<0) { volume<float> newvol; return newvol; }
    if ( vol4.tsize() != vec.Nrows() )
      {
	std::cerr << "ERROR::Time series length differs from vector length in"
	     << " dotproductvol()" << std::endl;
	volume<float> newvol; return newvol;
      }
    volume<float> vol4copy;
    copyconvert(vol4[vol4.mint()],vol4copy);
    volume<float> DotVol(vol4copy);
    DotVol *= (float) vec(1);
    for (int n=vol4.mint() + 1; n <= vol4.maxt(); n++) {
      copyconvert(vol4[n],vol4copy);
      DotVol += (vol4copy * (float) vec(1 + n - vol4.mint()));
    }
    return DotVol;
  }


  template <class T>
  void pad(const volume<T>& vol, volume<T>& paddedvol)
    {
      // The default type of padding is central padding
      if ( ( paddedvol.xsize() < vol.xsize() ) || ( paddedvol.ysize() < vol.ysize() ) || ( paddedvol.zsize() < vol.zsize() ) ) {
	imthrow("Cannot pad when target volume is smaller than original",7);
      }
      int64_t offx = ( paddedvol.xsize() - vol.xsize() ) / 2;
      int64_t offy = ( paddedvol.ysize() - vol.ysize() ) / 2;
      int64_t offz = ( paddedvol.zsize() - vol.zsize() ) / 2;
      pad(vol,paddedvol,offx,offy,offz);
    }

  template <class T>
  double dotproduct(const volume<T>& vol1,
                    const volume<T>& vol2)
  {
    if (!samesize(vol1,vol2,true)) imthrow("dotproduct: Image dimension mismatch",99);

    double rval = 0.0;
    for (typename volume<T>::fast_const_iterator it1=vol1.fbegin(), it_end=vol1.fend(), it2=vol2.fbegin(); it1 != it_end; ++it1, ++it2) {
      rval += static_cast<double>((*it1)*(*it2));
    }
    return(rval);
  }

  template <class T, class S>
  double dotproduct(const volume<T>& vol1,
                    const volume<T>& vol2,
                    const volume<S>& mask)
  {
    if (!samesize(vol1,vol2,true) || !samesize(vol1,mask,true)) imthrow("dotproduct: Image dimension mismatch",99);

    double rval = 0.0;
    typename volume<S>::fast_const_iterator  itm = mask.fbegin();
    for (typename volume<T>::fast_const_iterator it1=vol1.fbegin(), it_end=vol1.fend(), it2=vol2.fbegin(); it1 != it_end; ++it1, ++it2, ++itm) {
      if (*itm > 0.5) {
        rval += static_cast<double>((*it1)*(*it2));
      }
    }
    return(rval);
  }

  template <class T, class S>
  double dotproduct(const volume<T>& vol1,
                    const volume<T>& vol2,
                    const volume<S>  *mask)
  {
    if (mask) return(dotproduct(vol1,vol2,*mask));
    else return(dotproduct(vol1,vol2));
  }

  template <class T>
  double powerdotproduct(const volume<T>&  vol1,
                         unsigned int      n1,
                         const volume<T>&  vol2,
                         unsigned int      n2)
  {
    if (!samesize(vol1,vol2,true)) imthrow("powerdotproduct: Image dimension mismatch",99);

    double rval = 0.0;
    for (typename volume<T>::fast_const_iterator it1=vol1.fbegin(), it_end=vol1.fend(), it2=vol2.fbegin(); it1 != it_end; ++it1, ++it2) {
      double val1 = 1.0;
      for (unsigned int i=0; i<n1; i++) val1 *= *it1;
      double val2 = 1.0;
      for (unsigned int i=0; i<n2; i++) val2 *= *it2;
      rval += static_cast<double>(val1*val2);
    }
    return(rval);
  }

  template <class T, class S>
  double powerdotproduct(const volume<T>&  vol1,
                         unsigned int      n1,
                         const volume<T>&  vol2,
                         unsigned int      n2,
                         const volume<S>&  mask)
  {
    if (!samesize(vol1,vol2,true) || !samesize(vol1,mask,true)) imthrow("powerdotproduct: Image dimension mismatch",99);

    double rval = 0.0;
    typename volume<S>::fast_const_iterator itm = mask.fbegin();
    for (typename volume<T>::fast_const_iterator it1=vol1.fbegin(), it_end=vol1.fend(), it2=vol2.fbegin(); it1 != it_end; ++it1, ++it2, ++itm) {
      if (*itm > 0.5) {
        double val1 = 1.0;
        for (unsigned int i=0; i<n1; i++) val1 *= *it1;
        double val2 = 1.0;
        for (unsigned int i=0; i<n2; i++) val2 *= *it2;
        rval += static_cast<double>(val1*val2);
      }
    }
    return(rval);
  }

  template <class T, class S>
  double powerdotproduct(const volume<T>&  vol1,
                         unsigned int      n1,
                         const volume<T>&  vol2,
                         unsigned int      n2,
                         const volume<S>   *mask)
  {
    if (mask) return(powerdotproduct(vol1,n1,vol2,n2,*mask));
    else return(powerdotproduct(vol1,n1,vol2,n2));
  }

  template <class T>
  double mean(const volume<T>& vol)
  {
    double rval = sum(vol);
    rval /= static_cast<double>( vol.totalElements() );
    return(rval);
  }

  template <class T, class S>
  double mean(const volume<T>& vol,
              const volume<S>& mask)
  {
    if (!samesize(vol,mask,true)) imthrow("mean: Image-Mask dimension mismatch",99);

    double rval(0);
    int64_t n(0);
    typename volume<S>::fast_const_iterator  itm = mask.fbegin();
    for (typename volume<T>::fast_const_iterator it=vol.fbegin(), it_end=vol.fend(); it != it_end; ++it, ++itm) {
      if (*itm > 0.5) {
        n++;
        rval += static_cast<double>(*it);
      }
    }
    rval /= static_cast<double>(n);
    return(rval);
  }

  template <class T, class S>
  double mean(const volume<T>& vol,
              const volume<S> *mask)
  {
    if (mask) return(mean(vol,*mask));
    else return(mean(vol));
  }

  template <class T>
  double sum(const volume<T>& vol)
  {
    double rval(0);
    for (typename volume<T>::fast_const_iterator it=vol.fbegin(), it_end=vol.fend(); it != it_end; ++it) {
      rval += static_cast<double>(*it);
    }
    return(rval);
  }

  template <class T, class S>
  double sum(const volume<T>& vol,
             const volume<S>& mask)
  {
    if (!samesize(vol,mask,true)) imthrow("sum: Image-Mask dimension mismatch",99);

    double rval(0);
    typename volume<S>::fast_const_iterator  itm = mask.fbegin();
    for (typename volume<T>::fast_const_iterator it=vol.fbegin(), it_end=vol.fend(); it != it_end; ++it, ++itm) {
      if (*itm > 0.5) {
        rval += static_cast<double>(*it);
      }
    }
    return(rval);
  }

  template <class T, class S>
  double sum(const volume<T>& vol,
              const volume<S> *mask)
  {
    if (mask) return(sum(vol,*mask));
    else return(sum(vol));
  }

  template <class T, class S>
  volume<T> divide(const volume<T>& numervol, const volume<T>& denomvol,
		   const volume<S>& mask)
  {
    if ((!samesize(numervol,denomvol)) || (!samesize(mask,denomvol,SUBSET))) {
      imthrow("Attempted to divide images of different sizes",3);
    }
    volume<T> resvol(numervol);

    typename volume<S>::fast_const_iterator maskIt = mask.fbegin(), maskEnd=mask.fend();
    typename volume<T>::nonsafe_fast_iterator it=resvol.nsfbegin();
    for (typename volume<T>::fast_const_iterator itd=denomvol.fbegin(), itd_end=denomvol.fend(); itd != itd_end; ++it, ++itd, maskIt++) {
      if ( *maskIt != 0 ) {
	*it /= *itd;
      } else {
	*it = 0;
      }
      if ( maskIt == maskEnd )
    	maskIt=mask.fbegin();
    }
    return resvol;
  }



template <class T, class S>
  volume<T> mask_volume( const volume<T>& invol, const volume<S>& mask)
  {
    if ( !samesize(invol,mask,SUBSET)) {
      imthrow("Attempted to mask_volume with wrong sized mask",3);
    }
    volume<T> resvol(invol, false);

    typename volume<S>::fast_const_iterator maskIt = mask.fbegin(), maskEnd=mask.fend();
    typename volume<T>::nonsafe_fast_iterator itr=resvol.nsfbegin();
    for (typename volume<T>::fast_const_iterator it=invol.fbegin(), it_end=invol.fend(); it != it_end; ++it, ++itr, maskIt++) {
      if ( *maskIt != 0 ) {
	*itr = *it;
      } else {
	*itr = 0;
      }
      if ( maskIt == maskEnd )
    	maskIt=mask.fbegin();
    }
    return resvol;
  }


  // AFFINE TRANSFORM
 template <class T>
 void raw_affine_transform(const volume<T>& vin, volume<T>& vout,
			   const NEWMAT::Matrix& aff);

  template <class T>
  void affine_transform_mask(const volume<T>& vin, volume<T>& vout,
			     const NEWMAT::Matrix& aff, float padding, const T padval);

  template <class T>
  volume<T> affine_transform_mask(const volume<T>& vin, const volume<T>& vout,
				  const NEWMAT::Matrix& aff, float padding)
    {
      volume<T> affmask;
      affmask = vout;
      affmask = (T) 1;
      affine_transform_mask(vin,affmask,aff,padding,(T) 0);
      return affmask;
    }


  template <class T>
  void affine_transform(const volume<T>& vin, volume<T>& vout,
			const NEWMAT::Matrix& aff, float paddingsize, bool set_backgnd)
  {
    T padval=0;
    extrapolation oldex;
    if (set_backgnd) {
      padval = vin.getpadvalue();
      oldex = vin.getextrapolationmethod();
      vin.setpadvalue(vin.backgroundval());
      vin.setextrapolationmethod(extraslice);
    }


    raw_affine_transform(vin,vout,aff);
    // now mask the output to eliminate streaks formed by the sinc interp...
    affine_transform_mask(vin,vout,aff,paddingsize,vin.getpadvalue());
    if (set_backgnd) {
      vin.setpadvalue(padval);
      vin.setextrapolationmethod(oldex);
    }
  }

  template <class T>
  void affine_transform(const volume<T>& vin, volume<T>& vout,
 	                         const NEWMAT::Matrix& aff, float paddingsize)
  {
    affine_transform(vin,vout,aff,paddingsize,true);
  }

  template <class T>
  void affine_transform(const volume<T>& vin, volume<T>& vout,
			const NEWMAT::Matrix& aff, interpolation interptype,
			float paddingsize)
    {
      interpolation oldinterp;
      oldinterp = vin.getinterpolationmethod();
      vin.setinterpolationmethod(interptype);
      affine_transform(vin,vout,aff,paddingsize);
      vin.setinterpolationmethod(oldinterp);
    }



  ///////////////////////////////////////////////////////////////////////////
  // CONVOLVE
    template <class T, class S>
    volume<T> susan_convolve(const volume<T> source, const volume<S>& kernel, const float sigmabsq, const bool use_median, int num_usan,volume<T>* usan_area,volume<T> usan_vol1,const float sigmab1sq,volume<T> usan_vol2,const float sigmab2sq)
    //template <class T, class S, class U, class V, class W>
    //volume<T> susan_convolve(const volume<T>& source, const volume<S>& kernel, const float sigmabsq, const bool use_median, int num_usan,volume<U>* usan_area = new volume<T>(1,1,1),const volume<V>& usan_vol1=volume<T>(1,1,1),const float sigmab1sq=0,const volume<W>& usan_vol2 = volume<T>(1,1,1),const float sigmab2sq=0)
    //Note that the commented out declaration won't work with the optional arguements (since U,V,W need to be defined in call...). Code is provided for possible
    //future improvments, as it is all usans are templated as input
{
//need to use a pointer for usan_area as creating a default parameter for a pass-by-reference gives
//a "assignment to tempory memory" warning in gcc
//default values for lut1 etc
  if (    (( (kernel.maxz() - kernel.minz()) % 2)==1) ||
	  (( (kernel.maxy() - kernel.miny()) % 2)==1) ||
	  (( (kernel.maxx() - kernel.minx()) % 2)==1) )
	  std::cerr << "WARNING:: Off-centre convolution being performed as kernel has even dimensions" << std::endl;
  if ((num_usan>=1 && !samesize(source,usan_vol1)) || (num_usan>=2 && !samesize(source,usan_vol2)) || (num_usan>=1 && !samesize(source,*usan_area))      )
  {
    std::cerr << "Warning: an external usan or output usan is not the same size as the source image, reverting to num_usans=0 mode" << std::endl;
    num_usan=0;
  }
  int midx, midy, midz,lz,uz,lx,ux,ly,uy,lutsize=16384;
  lz=source.minz();
  uz=source.maxz();
  lx=source.minx();
  ux=source.maxx();
  ly=source.miny();
  uy=source.maxy();
  volume<T> result(source);
  midz=(kernel.maxz() - kernel.minz())/2;
  midy=(kernel.maxy() - kernel.miny())/2;
  midx=(kernel.maxx() - kernel.minx())/2;
  //generate look up table
  float range1=1,range2=1,range = (source.max() - source.min())/(float)lutsize;
  float **lut=new float *[3];
  for(int i=0;i<=num_usan;i++) lut[i]=new float[2*lutsize+1];
  for(int i=0;i<=num_usan;i++) lut[i]+=lutsize;
  for (int i=0;i<=lutsize;i++) lut[0][-i]=lut[0][i]= exp(-pow(i*range,2.0)/sigmabsq);
  if (num_usan>=1)
  {
     range1= (usan_vol1.max() - usan_vol1.min())/(float)lutsize;
     for (int i=0;i<=lutsize;i++)   lut[1][-i]=lut[1][i]= exp(-pow(i*range1,2.0)/sigmab1sq);
  }
  if (num_usan>=2)
  {
     range2= (usan_vol2.max() - usan_vol2.min())/(float)lutsize;
     for (int i=0;i<=lutsize;i++)   lut[2][-i]=lut[2][i]= exp(-pow(i*range2,2.0)/sigmab2sq);
  }
  NEWMAT::ColumnVector mediankernel((kernel.zsize()>1)?26:8); // cube or square, minus central voxel
  int medoffst=((kernel.zsize()>1)?1:0);
  for (int z=lz; z<=uz; z++)
    for (int y=ly; y<=uy; y++)
      for (int x=lx; x<=ux; x++)
      {
	 int xmin=x-midx,ymin=y-midy,zmin=z-midz;
	 int xmax=MIN(x+midx,ux),ymax=MIN(y+midy,uy),zmax=MIN(z+midz,uz);
         float num=0, denom=0,center_val1=0,center_val2=0,factor;
         float center_val=source.value(x,y,z);
         if (num_usan>=1) center_val1=usan_vol1.value(x,y,z);
         if (num_usan>=2) center_val2=usan_vol2.value(x,y,z);
	 for(int mz=MAX(zmin,lz); mz<=zmax; mz++)
	   for(int my=MAX(ymin,ly); my<=ymax; my++)
	     for(int mx=MAX(xmin,lx); mx<=xmax; mx++)
	       if ((factor=(float)kernel.value(mx-xmin,my-ymin,mz-zmin)))
	       {
		 if (num_usan==0) factor*= lut[0][(int)((source.value(mx,my,mz)-center_val)/range)];
		 else factor*= lut[1][(int)((usan_vol1.value(mx,my,mz)-center_val1)/range1)];
		 if (num_usan>=2) factor*=lut[2][(int)((usan_vol2.value(mx,my,mz)-center_val2)/range2)];
		 num+=source.value(mx,my,mz) * factor;
		 denom+=factor;
	       }
	 if (num_usan>=1) usan_area->value(x,y,z)=(T) denom;
	     if (use_median && denom<1.5)
             {
               int count=1;
               for(int x2=MAX(x-1,lx);x2<=MIN(x+1,ux);x2++)
		 for(int y2=MAX(y-1,ly);y2<=MIN(y+1,uy);y2++)
		   for(int z2=MAX(z-medoffst,lz);z2<=MIN(z+medoffst,uz);z2++)
		     if ( (x2-x) || (y2-y) || (z2-z) ) mediankernel(count++)=source.value(x2,y2,z2);
	       NEWMAT::ColumnVector subkernel = mediankernel.SubMatrix(1,count-1,1,1);
	       SortAscending(subkernel);
	       result(x,y,z) = (T)((subkernel(count/2)+subkernel((count+1)/2))/2.0);
	     }
             else result.value(x,y,z)=(T) (num/denom);
      }
  for(int i=0;i<=num_usan;i++) delete[] (lut[i]-lutsize);
  delete[] lut;
  return result;
}

  template <class T, class S>
  volume<T> convolve(const volume<T>& source, const volume<S>& kernel)
    {
      extrapolation oldex = source.getextrapolationmethod();
      int offset=0;
      if ((oldex==boundsassert) || (oldex==boundsexception))
	{ source.setextrapolationmethod(constpad); }
      volume<T> result(source);
      if ( (kernel.zsize()%2==0) || (kernel.ysize()%2==0) || (kernel.xsize()%2==0) )
	{
	  std::cerr << "WARNING:: Off-centre convolution being performed as kernel"
	       << " has even dimensions" << std::endl;
          //offset=2;
	  //offset gives correction to convolve to match results with fft for even kernel
          //not even kernel with normalise (e.g. -fmean in fslmaths) still has edge problems
	}
      int midx, midy, midz;
      midz=(kernel.zsize()-1)/2 + offset;
      midy=(kernel.ysize()-1)/2 + offset;
      midx=(kernel.xsize()-1)/2 + offset;

      float val;
      for (int z=0; z<result.zsize(); z++)
	for (int y=0; y<result.ysize(); y++)
	  for (int x=0; x<result.xsize(); x++) {
	    val=0.0;
	    for (int mz=0; mz<kernel.zsize(); mz++)
	      for (int my=0; my<kernel.ysize(); my++)
		for (int mx=0; mx<kernel.xsize(); mx++) {
		  val+=source(x+mx-midx,y+my-midy,z+mz-midz) * kernel(mx,my,mz);
		}
	    result(x,y,z)=(T) val;
	  }


      source.setextrapolationmethod(oldex);
      return result;
    }

  template <class T, class S, class M>
  volume<T> convolve(const volume<T>& source, const volume<S>& kernel,
		     const volume<M>& mask, bool ignoremask, bool renormalise)
    {
      if (!ignoremask && !samesize(mask, source))
	imthrow("convolve: mask and source are not the same size",10);

      if (    (( (kernel.maxz() - kernel.minz()) % 2)==1) ||
	      (( (kernel.maxy() - kernel.miny()) % 2)==1) ||
	      (( (kernel.maxx() - kernel.minx()) % 2)==1) )
	{
	  std::cerr << "WARNING:: Off-centre convolution being performed as kernel"
	       << " has even dimensions" << std::endl;
	}
      volume<T> result(source);
      int midx, midy, midz;
      int lz,uz,lx,ux,ly,uy;
      lz=result.minz();
      uz=result.maxz();
      lx=result.minx();
      ux=result.maxx();
      ly=result.miny();
      uy=result.maxy();
      midz=(kernel.maxz() - kernel.minz())/2;
      midy=(kernel.maxy() - kernel.miny())/2;
      midx=(kernel.maxx() - kernel.minx())/2;
      for (int z=lz; z<=uz; z++)
	for (int y=ly; y<=uy; y++)
	  for (int x=lx; x<=ux; x++)
	    if (ignoremask || mask(x,y,z)>0.5)
            {
	      float val(0),norm(0);
              int x3,y3,z3;
              x3=x-midx;
              y3=y-midy;
              z3=z-midz;
	      for (int mz=kernel.minz(); mz<=kernel.maxz(); mz++)
		for (int my=kernel.miny(); my<=kernel.maxy(); my++)
		  for (int mx=kernel.minx(); mx<=kernel.maxx(); mx++)
                  {
                    int x2,y2,z2;
                    x2=x3+mx;
                    y2=y3+my;
                    z2=z3+mz;
 		    if ((ignoremask && (x2<=ux && x2>=lx && y2<=uy && y2>=ly && z2<=uz && z2>=lz))  || (!ignoremask && mask(x2,y2,z2)>0.5))
                    {
		      val+=source.value(x2,y2,z2) * kernel.value(mx,my,mz);
		      norm+=kernel.value(mx,my,mz);
		    }
		  }
	      if (renormalise && fabs(norm)>1e-12) result.value(x,y,z)=(T) (val/norm);
	      else result.value(x,y,z)=(T) val;
	    }
      return result;
    }

  template <class T>
  volume<T> convolve_separable(const volume<T>& source,
			       const NEWMAT::ColumnVector& kernelx,
			       const NEWMAT::ColumnVector& kernely,
			       const NEWMAT::ColumnVector& kernelz)
    {
      volume<T> result(source);
      volume<double> kerx(kernelx.Nrows(),1,1);
      volume<double> kery(1,kernely.Nrows(),1);
      volume<double> kerz(1,1,kernelz.Nrows());
      for (int n=1; n<=kernelx.Nrows(); n++)  kerx.value(n-1,0,0) = kernelx(n);
      for (int n=1; n<=kernely.Nrows(); n++)  kery.value(0,n-1,0) = kernely(n);
      for (int n=1; n<=kernelz.Nrows(); n++)  kerz.value(0,0,n-1) = kernelz(n);
      result = convolve(result,kerx);
      result = convolve(result,kery);
      result = convolve(result,kerz);
      return result;
    }

  template <class T, class M>
  volume<T> convolve_separable(const volume<T>& source,
			       const NEWMAT::ColumnVector& kernelx,
			       const NEWMAT::ColumnVector& kernely,
			       const NEWMAT::ColumnVector& kernelz,
			       const volume<M>& mask, bool ignoremask, bool renormalise)
    {
      volume<T> result(source);
      volume<double> kerx(kernelx.Nrows(),1,1);
      volume<double> kery(1,kernely.Nrows(),1);
      volume<double> kerz(1,1,kernelz.Nrows());
      for (int n=1; n<=kernelx.Nrows(); n++)  kerx.value(n-1,0,0) = kernelx(n);
      for (int n=1; n<=kernely.Nrows(); n++)  kery.value(0,n-1,0) = kernely(n);
      for (int n=1; n<=kernelz.Nrows(); n++)  kerz.value(0,0,n-1) = kernelz(n);
      result = convolve(result,kerx,mask,ignoremask,renormalise);
      result = convolve(result,kery,mask,ignoremask,renormalise);
      result = convolve(result,kerz,mask,ignoremask,renormalise);
      return result;
    }

 ///////////////////////////////////////////////////////////////////////////
  // GENERAL DILATION INCLUDING MEDIAN FILTERING

  template <class T, class S>
  volume<T> morphfilter(const volume<T>& source, const volume<S>& kernel,
			const std::string& filtertype)
    {
      // implements a whole range of filtering, set by the string filtertype:
      //  dilateM (mean), dilateD (mode), median, max or dilate, min or erode,
      //  erodeS (set to zero if any neighbour is zero)
      extrapolation oldex = source.getextrapolationmethod();
      if ((oldex==boundsassert) || (oldex==boundsexception))
	{ source.setextrapolationmethod(constpad); }
      volume<T> result(source);
      result = 0;

      int nker;
      {
	volume<S> dummy(kernel);
	dummy.binarise((S)0.5);  //new cast to avoid warnings for int-type templates when compiling
	nker = (int) dummy.sum();
      }

      int midx, midy, midz;
      if (    (( (kernel.maxz() - kernel.minz()) % 2)==1) ||
	      (( (kernel.maxy() - kernel.miny()) % 2)==1) ||
	      (( (kernel.maxx() - kernel.minx()) % 2)==1) )
	{
	  std::cerr << "WARNING:: Off-centre morphfilter being performed as kernel"
	       << " has even dimensions" << std::endl;
	}
      midz=(kernel.maxz() - kernel.minz())/2;
      midy=(kernel.maxy() - kernel.miny())/2;
      midx=(kernel.maxx() - kernel.minx())/2;
      int count=1;

      for (int z=result.minz(); z<=result.maxz(); z++) {
	for (int y=result.miny(); y<=result.maxy(); y++) {
	  for (int x=result.minx(); x<=result.maxx(); x++) {
        NEWMAT::ColumnVector vals(nker);
	    count=1;
	    for (int mz=MISCMATHS::Max(kernel.minz(),result.minz()-z+midz);
		 mz<=MISCMATHS::Min(kernel.maxz(),result.maxz()-z+midz); mz++) {
	      for (int my=MISCMATHS::Max(kernel.miny(),result.miny()-y+midy);
		   my<=MISCMATHS::Min(kernel.maxy(),result.maxy()-y+midy); my++) {
		for (int mx=MISCMATHS::Max(kernel.minx(),result.minx()-x+midx);
		     mx<=MISCMATHS::Min(kernel.maxx(),result.maxx()-x+midx); mx++) {
		  if (kernel(mx,my,mz)>0.5) {
		    if ((filtertype!="dilateM" && filtertype!="dilateD") || source(x+mx-midx,y+my-midy,z+mz-midz)) vals(count++) = source(x+mx-midx,y+my-midy,z+mz-midz);
		  }
		}
	      }
	    }
	    if (count>1) {
	      NEWMAT::ColumnVector littlevals;
	      littlevals = vals.SubMatrix(1,count-1,1,1);
	      if (filtertype=="median") {
		SortAscending(littlevals);                         //count/2 works for odd kernel (even count) count+1 gives edge compatibility
		result(x,y,z) = (T)littlevals(MISCMATHS::Max(1,(count+1)/2)); //with steves IP, otherwise gives the IP median-1 element
	      } else if ((filtertype=="max") || (filtertype=="dilate"))
            result(x,y,z) = (T)littlevals.Maximum();
	        else if ((filtertype=="min") || (filtertype=="erode"))
            result(x,y,z) = (T)littlevals.Minimum();
          else if (filtertype=="erodeS") {
            if (source(x,y,z)!=0 && littlevals.Minimum()==0) result(x,y,z) = 0; else result(x,y,z) = source(x,y,z);}
          else if (filtertype=="dilateM") {
            if (source(x,y,z)==0) result(x,y,z) = (T)(littlevals.Sum()/--count); else result(x,y,z) = source(x,y,z);}
          else if (filtertype=="dilateD") {
            if (source(x,y,z)==0) {
		          SortDescending(littlevals);
              double mode=littlevals(1);
              int maxn=1;
              int currentn=1;
              for(int i=2;i<count;i++) {
                currentn = (littlevals(i) == littlevals(i-1)) ? ++currentn : 1;
                if (currentn>maxn) {
                  mode=littlevals(i);
                  maxn=currentn;
                }
              }
		          result(x,y,z) = (T)mode;
            } else result(x,y,z) = source(x,y,z);
		      }
	        else imthrow("morphfilter:: Filter type " + filtertype + "unsupported",7);
	    }   else result(x,y,z) = source(x,y,z);  // THE DEFAULT CASE (leave alone)
	  }
	}
      }
      source.setextrapolationmethod(oldex);
      return result;
    }

  template <class T>
  double dilateval(const volume<T>& im, const volume<T>& mask, int x, int y, int z)
  {
    double sum=0;
    int n=0;
    for (int zz=MISCMATHS::Max(0,z-1); zz<=MISCMATHS::Min(z+1,im.maxz()); zz++) {
      for (int yy=MISCMATHS::Max(0,y-1); yy<=MISCMATHS::Min(y+1,im.maxy()); yy++) {
	for (int xx=MISCMATHS::Max(0,x-1); xx<=MISCMATHS::Min(x+1,im.maxx()); xx++) {
	  if ((mask(xx,yy,zz)>(T) 0.5) && !((xx==x) && (yy==y) && (zz==z))) {
	    sum += im(xx,yy,zz);
	    n++;
	  }
	}
      }
    }
    return sum/(MISCMATHS::Max(n,1));
  }


  template <class T>
  bool ext_edge(const volume<T>& mask, int x, int y, int z)
  {
    if (mask(x,y,z)>(T) 0.5) { return false; }
    for (int zz=MISCMATHS::Max(0,z-1); zz<=MISCMATHS::Min(z+1,mask.maxz()); zz++) {
      for (int yy=MISCMATHS::Max(0,y-1); yy<=MISCMATHS::Min(y+1,mask.maxy()); yy++) {
	for (int xx=MISCMATHS::Max(0,x-1); xx<=MISCMATHS::Min(x+1,mask.maxx()); xx++) {
	  if ((mask(xx,yy,zz)>(T) 0.5) && !((xx==x) && (yy==y) && (zz==z))) {
	    return true;
	  }
	}
      }
    }
    return false;
  }

  struct vec3_temp { int x; int y; int z; };


  template <class T>
  volume<T> dilall(const volume<T>& input, volume<T>& mask)
  {
    volume<T> im(input);
    if (!samesize(input,mask)) { std::cerr << "ERROR::dilall::image are not the same size" << std::endl; }
    std::deque<vec3_temp> ptlist;
    vec3_temp v, newv;
    // initial pass
    for (int z=0; z<=im.maxz(); z++) {
      for (int y=0; y<=im.maxy(); y++) {
	for (int x=0; x<=im.maxx(); x++) {
	  if (ext_edge(mask,x,y,z)) {
	    v.x=x; v.y=y; v.z=z;
	    ptlist.push_front(v);
	  }
	}
      }
    }
    while (!ptlist.empty()) {
      v = ptlist.back();
      ptlist.pop_back();
      if (mask(v.x,v.y,v.z)<=(T) 0.5) {  // only do it if the voxel is still unset
	im(v.x,v.y,v.z)=dilateval(im,mask,v.x,v.y,v.z);
	mask(v.x,v.y,v.z)=(T)1;
	// check neighbours and add them to the list if necessary
	for (int zz=MISCMATHS::Max(0,v.z-1); zz<=MISCMATHS::Min(v.z+1,im.maxz()); zz++) {
	  for (int yy=MISCMATHS::Max(0,v.y-1); yy<=MISCMATHS::Min(v.y+1,im.maxy()); yy++) {
	    for (int xx=MISCMATHS::Max(0,v.x-1); xx<=MISCMATHS::Min(v.x+1,im.maxx()); xx++) {
	      if (ext_edge(mask,xx,yy,zz)) { newv.x=xx; newv.y=yy; newv.z=zz; ptlist.push_front(newv); }
	    }
	  }
	}
      }
    }
    return im;
  }


  template <class T>
  volume<T> fill_holes(const volume<T>& im, const int connectivity)
  {
    volume<int> mask;
    volume<T> workingImage(im);
    std::set<int> edgelabs{0};  //For later we include 0 as an edge-label

    workingImage.binarise(0,0,inclusive,false);

    mask = connected_components(workingImage,connectivity);
    // Store all labels at image edges ( including repeated 0s )
    // Edge defined as at least one of x,y,z at limit
    // The modulo operator will only be 0 at either edge
    for (int z=0; z<=im.maxz(); z++)
      for (int y=0; y<=im.maxy(); y++)
        for (int x=0; x<=im.maxx(); x++)
	        if ( !(x%im.maxx() && y%im.maxy() && z%im.maxz()))
            edgelabs.insert(mask(x,y,z));

    //Fill any masked voxel which isn't connected to the edge ( i.e. a "hole" )
    for (int z=0; z<=im.maxz(); z++)
	    for (int y=0; y<=im.maxy(); y++)
	      for (int x=0; x<=im.maxx(); x++)
          workingImage(x,y,z) = ( edgelabs.find(mask(x,y,z)) == edgelabs.end() ) ? 1 : im(x,y,z);

    return workingImage;
  }

 ///////////////////////////////////////////////////////////////////////////
  // RESAMPLE

  template <class T>
  volume<T> upsample_by_2(const volume<T>& lowresvol, bool centred)
  {
    volume<T> res;
    upsample_by_2(res,lowresvol,centred);
    return res;
  }



  ///////////////////////////////////////////////////////////////////////////
  // BLURRING
  template <class T>
  volume<T> blur(const volume<T>& source, const NEWMAT::ColumnVector& resel_size)
    {
      NEWMAT::ColumnVector bmaskx, bmasky, bmaskz;
      make_blur_mask(bmaskx,resel_size(1),source.xdim());
      make_blur_mask(bmasky,resel_size(2),source.ydim());
      make_blur_mask(bmaskz,resel_size(3),source.zdim());
      return convolve_separable(source,bmaskx,bmasky,bmaskz);
    }


  template <class T>
  volume<T> blur(const volume<T>& source, float iso_resel_size)
    {
      NEWMAT::ColumnVector resel_size(3);
      resel_size = iso_resel_size;
      return blur(source,resel_size);
    }

  template <class T>
  volume<T> smooth(const volume<T>& source, float sigma_mm)
    {
      float sigmax, sigmay, sigmaz;
      sigmax = sigma_mm/source.xdim();
      sigmay = sigma_mm/source.ydim();
      sigmaz = sigma_mm/source.zdim();
      int nx=((int) (sigmax-0.001))*2 + 3;
      int ny=((int) (sigmay-0.001))*2 + 3;
      int nz=((int) (sigmaz-0.001))*2 + 3;
      NEWMAT::ColumnVector kernelx, kernely, kernelz;
      kernelx = gaussian_kernel1D(sigmax,nx);
      kernely = gaussian_kernel1D(sigmay,ny);
      kernelz = gaussian_kernel1D(sigmaz,nz);
      return convolve_separable(source,kernelx,kernely,kernelz);
    }




  template <class T>
  volume<T> smooth2D(const volume<T>& source, float sigma_mm, int nulldir=3)
    {
      float sigmax, sigmay, sigmaz;
      sigmax = sigma_mm/source.xdim();
      sigmay = sigma_mm/source.ydim();
      sigmaz = sigma_mm/source.zdim();
      int nx=((int) (sigmax-0.001))*2 + 3;
      int ny=((int) (sigmay-0.001))*2 + 3;
      int nz=((int) (sigmaz-0.001))*2 + 3;
      NEWMAT::ColumnVector kernelx, kernely, kernelz, nullker(1);
      kernelx = gaussian_kernel1D(sigmax,nx);
      kernely = gaussian_kernel1D(sigmay,ny);
      kernelz = gaussian_kernel1D(sigmaz,nz);
      nullker = 1;
      if (nulldir==1) {
	return convolve_separable(source,nullker,kernely,kernelz);
      } else if (nulldir==2) {
	return convolve_separable(source,kernelx,nullker,kernelz);
      } else {
	// Smoothing in the x-y plane is the default!
	return convolve_separable(source,kernelx,kernely,nullker);
      }
    }

template <class T, class S>
  int insertpart(volume<T>& v1, const volume<S>& v2)  //N.B. This has superficial similarities
{                                                     //to copyconvert, but is in fact quite
  for (int z=v2.minz(); z<=v2.maxz(); z++) {          //different...
    for (int y=v2.miny(); y<=v2.maxy(); y++) {
      for (int x=v2.minx(); x<=v2.maxx(); x++) {
	v1(x,y,z)=(T) v2(x,y,z);
      }
    }
  }
  return 0;
}


 template <class T, class S,class U>
volume<S> extractpart(const volume<T>& v1, const volume<S>& v2, const volume<U>& kernel)
{
  volume<S> vout=v2;
  vout = (S) 0.0;
  int kxoff = (kernel.xsize()-1)/2;
  int kyoff = (kernel.ysize()-1)/2;
  int kzoff = (kernel.zsize()-1)/2;
  for (int z=v2.minz(); z<=v2.maxz(); z++) {
    for (int y=v2.miny(); y<=v2.maxy(); y++) {
      for (int x=v2.minx(); x<=v2.maxx(); x++) {
	vout(x,y,z)=(S) v1(x+kxoff,y+kyoff,z+kzoff);
      }
    }
  }
  return vout;
}

////////////////////////////////////////////////////////////////////////////
// Efficient FFT-based convolve
template <class T, class S>
volume<T> efficient_convolve(const volume<T>& vin, const volume<S>& vker)
{
  bool usefft=true;
  // estimate calculation time for the two methods and pick the best
  float offt = 2 * vin.nvoxels() * fsllog2(2 * vin.nvoxels());
  float osum = (float)vin.nvoxels() * (float)vker.nvoxels();
  // float cast to avoud overflow for large int multiplication
  float scalefactor = 44;  // relative unit operation cost for fft vs sum
  usefft = (osum > offt * scalefactor);
  //cout << usefft << endl;
  if (usefft) {
    int sx = MISCMATHS::Max(vin.xsize(),vker.xsize())*2;
    int sy = MISCMATHS::Max(vin.ysize(),vker.ysize())*2;
    int sz = MISCMATHS::Max(vin.zsize(),vker.zsize())*2;
    complexvolume vif, vkf;
    vif.re().reinitialize(sx,sy,sz);
    //vif.re().copyproperties(vin);     Is this needed check with MJ...
    vif.re() = 0.0;
    vif.im() = vif.re();
    vkf = vif;
    insertpart(vif.re(),vin);
    insertpart(vkf.re(),vker);
    fft3(vif);
    fft3(vkf);
    vif *= vkf;
    ifft3(vif);
    return extractpart(vif.re(),vin,vker);
  } else return convolve(vin,vker);
}

  template <class T, class S>
    volume4D<T> generic_convolve(const volume4D<T>& source, const volume<S>& kernel, bool seperable, bool renormalise)
  {
      volume4D<T> result(source);
      if (seperable)
      {
        volume<double> kerx(kernel.xsize(),1,1);
        volume<double> kery(1,kernel.ysize(),1);
        volume<double> kerz(1,1,kernel.zsize());
        for (int n=0; n<kernel.xsize(); n++)  kerx.value(n,0,0) = kernel(n,kernel.ysize()/2,kernel.zsize()/2);
        for (int n=0; n<kernel.ysize(); n++)  kery.value(0,n,0) = kernel(kernel.xsize()/2,n,kernel.zsize()/2);
        for (int n=0; n<kernel.zsize(); n++)  kerz.value(0,0,n) = kernel(kernel.xsize()/2,kernel.ysize()/2,n);
        volume<T> mask(1,1,1);
	for (int t=source.mint(); t<=source.maxt(); t++)  result.replaceSubVolume(t,convolve(result.constSubVolume(t),kerx,mask,true,renormalise));
        for (int t=source.mint(); t<=source.maxt(); t++)  result.replaceSubVolume(t,convolve(result.constSubVolume(t),kery,mask,true,renormalise));
        for (int t=source.mint(); t<=source.maxt(); t++)  result.replaceSubVolume(t,convolve(result.constSubVolume(t),kerz,mask,true,renormalise));
      }
      else
      {
        volume<S> norm_kernel(kernel);
        if (kernel.sum()) norm_kernel/=kernel.sum();
	for (int t=source.mint(); t<=source.maxt(); t++) result.replaceSubVolume(t,efficient_convolve(source.constSubVolume(t),norm_kernel));
        result.copyproperties(source);
        if(renormalise)
	{
	  volume4D<T> unitary_mask(source);
          unitary_mask=1;
          for (int t=source.mint(); t<=source.maxt(); t++) unitary_mask.replaceSubVolume(t,efficient_convolve(unitary_mask.constSubVolume(t),norm_kernel));
          result/=unitary_mask;
        }
      }
    return result;
  }


  ///////////////////////////////////////////////////////////////////////////
  // GRADIENT


  template <class T>
  volume<float> gradient(const volume<T>& source)
    {
      // in voxel coordinates (not mm)
      volume<float> maskx,masky,maskz;
      make_grad_masks(maskx,masky,maskz);
      volume<float> grad(source);
      float valx, valy, valz;
      int midx, midy, midz;
      midz=maskx.xsize()/2;
      midy=maskx.ysize()/2;
      midx=maskx.zsize()/2;
      for (int z=0; z<grad.zsize(); z++) {
	for (int y=0; y<grad.ysize(); y++) {
	  for (int x=0; x<grad.xsize(); x++) {
	    valx=0.0; valy=0.0; valz=0.0;
	    for (int mz=-midz; mz<=midz; mz++) {
	      for (int my=-midy; my<=midy; my++) {
		for (int mx=-midx; mx<=midx; mx++) {
		  valx+=source(x+mx,y+my,z+mz) * maskx(mx+midx,my+midy,mz+midz);
		  valy+=source(x+mx,y+my,z+mz) * masky(mx+midx,my+midy,mz+midz);
		  valz+=source(x+mx,y+my,z+mz) * maskz(mx+midx,my+midy,mz+midz);
		}
	      }
	    }
	    grad(x,y,z)=sqrt(MISCMATHS::Sqr(valx) + MISCMATHS::Sqr(valy) + MISCMATHS::Sqr(valz));
	  }
	}
      }
      return grad;
    }



   template <class T>
   void gradient(const volume<T>& source,volume<float>& grad)
    {
      volume<float> maskx,masky,maskz;
      make_grad_masks(maskx,masky,maskz);
      grad.reinitialize(source.xsize(),source.ysize(),source.zsize(),3);
      copybasicproperties(source,grad);
      float valx, valy, valz;
      int midx, midy, midz;
      midz=maskx.xsize()/2;
      midy=maskx.ysize()/2;
      midx=maskx.zsize()/2;
      for (int z=0; z<grad.zsize(); z++) {
	for (int y=0; y<grad.ysize(); y++) {
	  for (int x=0; x<grad.xsize(); x++) {
	    valx=0.0; valy=0.0; valz=0.0;
	    for (int mz=-midz; mz<=midz; mz++) {
	      for (int my=-midy; my<=midy; my++) {
		for (int mx=-midx; mx<=midx; mx++) {
		  valx+=source(x+mx,y+my,z+mz) * maskx(mx+midx,my+midy,mz+midz);
		  valy+=source(x+mx,y+my,z+mz) * masky(mx+midx,my+midy,mz+midz);
		  valz+=source(x+mx,y+my,z+mz) * maskz(mx+midx,my+midy,mz+midz);
		}
	      }
	    }
	    grad(x,y,z,0)=valx;
	    grad(x,y,z,1)=valy;
	    grad(x,y,z,2)=valz;
	  }
	}
      }

    }



  template <class T>
  volume<float> lrxgrad(const volume<float>& im, const volume<T>& mask)
  {
    // calculates separate left and right gradients (with copying when at
    //  borders of the mask) or zero for points outside of the mask
    // returns the two as part of a volume4D (vol[0] = left grad, vol[1] = right grad)
    volume4D<float> grad(im);
    grad.addvolume(im);
    if (!samesize(im,mask)) imthrow("Mask and image not the same size",20);

    for (int z=0; z<im.zsize(); z++) {
      for (int y=0; y<im.ysize(); y++) {
	for (int x=0; x<im.xsize(); x++) {
	  // default 0, only overridden if inside mask and can calculate the gradients
	  if (mask(x,y,z)>0.5) {
	    // left gradient
	    if (x>0) {
	      if (mask(x-1,y,z)>0.5) {
		grad(x,y,z,0) = im(x,y,z) - im(x-1,y,z);
	      }
	    }
	    // right gradient
	    if (x<im.xsize()-1) {
	      if (mask(x+1,y,z)>0.5) {
		grad(x,y,z,1) = im(x+1,y,z) - im(x,y,z);
	      } else {
		// if couldn't calculate right grad, then copy left
		grad(x,y,z,1) = grad(x,y,z,0);
	      }
	    }
	    // if couldn't calculate left grad, then copy right
	    if ( (x>0) && (mask(x-1,y,z)<=0.5) ) {
	      grad(x,y,z,0) = grad(x,y,z,1);
	    }
	    // NB: if couldn't calculate either (but mask>0.5) then it is still 0
	  }
	}
      }
    }
    return grad;
  }


  template <class T>
  volume<float> lrygrad(const volume<float>& im, const volume<T>& mask)
  {
    // calculates separate left and right gradients (with copying when at
    //  borders of the mask) or zero for points outside of the mask
    // returns the two as part of a volume4D (vol[0] = left grad, vol[1] = right grad)
    volume4D<float> grad(im);
    grad.addvolume(im);
    if (!samesize(im,mask)) imthrow("Mask and image not the same size",20);

    for (int z=0; z<im.zsize(); z++) {
      for (int y=0; y<im.ysize(); y++) {
	for (int x=0; x<im.xsize(); x++) {
	  // default 0, only overridden if inside mask and can calculate the gradients
	  if (mask(x,y,z)>0.5) {
	    // left gradient
	    if (y>0) {
	      if (mask(x,y-1,z)>0.5) {
		grad(x,y,z,0) = im(x,y,z) - im(x,y-1,z);
	      }
	    }
	    // right gradient
	    if (y<im.ysize()-1) {
	      if (mask(x,y+1,z)>0.5) {
		grad(x,y,z,1) = im(x,y+1,z) - im(x,y,z);
	      } else {
		// if couldn't calculate right grad, then copy left
		grad(x,y,z,1) = grad(x,y,z,0);
	      }
	    }
	    // if couldn't calculate left grad, then copy right
	    if ( (y>0) && (mask(x,y-1,z)<=0.5) ) {
	      grad(x,y,z,0) = grad(x,y,z,1);
	    }
	    // NB: if couldn't calculate either (but mask>0.5) then it is still 0
	  }
	}
      }
    }
    return grad;
  }


  template <class T>
  volume<float> lrzgrad(const volume<float>& im, const volume<T>& mask)
  {
    // calculates separate left and right gradients (with copying when at
    //  borders of the mask) or zero for points outside of the mask
    // returns the two as part of a volume4D (vol[0] = left grad, vol[1] = right grad)
    volume4D<float> grad(im);
    grad.addvolume(im);
    if (!samesize(im,mask)) imthrow("Mask and image not the same size",20);
    for (int z=0; z<im.zsize(); z++) {
      for (int y=0; y<im.ysize(); y++) {
	for (int x=0; x<im.xsize(); x++) {
	  // default 0, only overridden if inside mask and can calculate the gradients
	  if (mask(x,y,z)>0.5) {
	    // left gradient
	    if (z>0) {
	      if (mask(x,y,z-1)>0.5) {
		grad(x,y,z,0) = im(x,y,z) - im(x,y,z-1);
	      }
	    }
	    // right gradient
	    if (z<im.zsize()-1) {
	      if (mask(x,y,z+1)>0.5) {
		grad(x,y,z,1) = im(x,y,z+1) - im(x,y,z);
	      } else {
		// if couldn't calculate right grad, then copy left
		grad(x,y,z,1) = grad(x,y,z,0);
	      }
	    }
	    // if couldn't calculate left grad, then copy right
	    if ( (z>0) && (mask(x,y,z-1)<=0.5) ) {
	      grad(x,y,z,0) = grad(x,y,z,1);
	    }
	    // NB: if couldn't calculate either (but mask>0.5) then it is still 0
	  }
	}
      }
    }
    return grad;
  }




   template <class T>
   volume<T> bandpass_temporal_filter(volume<T>& source,double hp_sigma, double lp_sigma)
      {
         int hp_mask_size_PLUS, lp_mask_size_PLUS, hp_mask_size_MINUS, lp_mask_size_MINUS;
         double *hp_exp=NULL, *lp_exp=NULL, *array, *array2;
         volume<T> result(source);

         if (hp_sigma<=0) hp_mask_size_MINUS=0;
         else hp_mask_size_MINUS=(int)(hp_sigma*3);   /* this isn't a linear filter, so small hard cutoffs at ends don't matter */
         hp_mask_size_PLUS=hp_mask_size_MINUS;
         if (lp_sigma<=0) lp_mask_size_MINUS=0;
         else lp_mask_size_MINUS=(int)(lp_sigma*20)+2; /* this will be small, so we might as well be careful */
         lp_mask_size_PLUS=lp_mask_size_MINUS;

         array=new double[source.tsize()];
         array2=new double[source.tsize()];

         if (hp_sigma>0)
         {
            hp_exp=new double[hp_mask_size_MINUS+hp_mask_size_PLUS+1];
            hp_exp+=hp_mask_size_MINUS;
            for(int t=-hp_mask_size_MINUS; t<=hp_mask_size_PLUS; t++)
            hp_exp[t] = exp( -0.5 * ((double)(t*t)) / (hp_sigma*hp_sigma) );
         }

         if (lp_sigma>0)
         {
            double total=0;
            lp_exp=new double[lp_mask_size_MINUS+lp_mask_size_PLUS+1];
            lp_exp+=lp_mask_size_MINUS;
            for(int t=-lp_mask_size_MINUS; t<=lp_mask_size_PLUS; t++)
            {
              lp_exp[t] = exp( -0.5 * ((double)(t*t)) / (lp_sigma*lp_sigma) );
              total += lp_exp[t];
            }

            for(int t=-lp_mask_size_MINUS; t<=lp_mask_size_PLUS; t++)
              lp_exp[t] /= total;
         }
         for(int z=0;z<source.zsize();z++)
           for(int y=0;y<source.ysize();y++)
	     for(int x=0;x<source.xsize();x++)
	     {
               for(int t=0; t<source.tsize(); t++)
		 array[t] = (double)source.value(x,y,z,t);
               if (hp_sigma>0)
               {
                 int done_c0=0;
                 double c0=0;
		 double mean(0);
                 for(int t=0; t<source.tsize(); t++)
                 {
                    int tt;
                    double c, w, A=0, B=0, C=0, D=0, N=0, tmpdenom;
                    for(tt=MAX(t-hp_mask_size_MINUS,0); tt<=MIN(t+hp_mask_size_PLUS,source.tsize()-1); tt++)
                    {
                      int dt=tt-t;
                      w = hp_exp[dt];
                      A += w * dt;
                      B += w * array[tt];
                      C += w * dt * dt;
                      D += w * dt * array[tt];
                      N += w;
                    }
                    tmpdenom=C*N-A*A;
                    if (tmpdenom!=0)
	            {
	               c = (B*C-A*D) / tmpdenom;
	               if (!done_c0)
	               {
	                 c0=c;
	                 done_c0=1;
	               }
	               array2[t] = c0 + array[t] - c;
	             }
	             else  array2[t] = array[t];
		    mean+=array2[t];
		 }
		 //Demean timeseries
		 mean/=source.tsize();
		 for(int t=0; t<source.tsize(); t++)
		   array2[t]-=mean;
                  memcpy(array,array2,sizeof(double)*source.tsize());
	       }
	       /* {{{ apply lowpass filter to 1D array */
               if (lp_sigma>0)
               {
                 for(int t=0; t<source.tsize(); t++)
                 {
                    double total=0;
		    double sum(0);
                    int tt;
                    for(tt=MAX(t-lp_mask_size_MINUS,0); tt<=MIN(t+lp_mask_size_PLUS,source.tsize()-1); tt++) {
		      total += array[tt] * lp_exp[tt-t];
		      sum+=lp_exp[tt-t];
		    }
		    if (sum>0)
		      array2[t] = total/sum;
		    else
		      array2[t] = total;
                 }
                 memcpy(array,array2,sizeof(double)*source.tsize());
	       }
	  /* {{{ write 1D array back to input 4D data */
               for(int t=0; t<source.tsize(); t++)
		 result.value(x,y,z,t)= (T)array[t] ;

	     }
	 return result;
      }


  ///////////////////////////////////////////////////////////////////////////
  // EDGE DETECTION

  /* detects closed contour edges in a volume as the zero crossings of the
     Laplacian of a Gaussian (LoG). This is implemented by convolving the
     data with a kernel formed by subtracting a Gaussian of radius sigma1
     from a second Gaussian kernel of radius sigma2 (sigma1 < sigma2) */

 template <class T>
  volume<T> log_edge_detect(const volume<T>& source,
			    float sigma1, float sigma2, int mode)
    {

      /*
       * mode:
       *   0 : 3D edge detection
       *   1 : 2D edge detection in X plane
       *   2 : 2D edge detection in Y plane
       *   3 : 2D edge detection in Z plane
       */

      int radius1;
      volume<float> log_kern, temp_kern;
      volume<T> log_result(source);
      volume<T> zero_crossing_result(source);
      zero_crossing_result = 0;

      radius1 = (int)(4*sigma2);
      if (mode==0) {
	log_kern = gaussian_kernel2D(sigma2, radius1);
	temp_kern = gaussian_kernel2D(sigma1, radius1);
      } else {
	log_kern = gaussian_kernel3D(sigma2, radius1);
	temp_kern = gaussian_kernel3D(sigma1, radius1);
      }

      log_kern -= temp_kern;

      log_result = convolve(source, log_kern);
      for(int t=0;t<log_result.tsize();t++)
        for(int z=1;z<log_result.zsize()-1;z++)
          for(int y=1;y<log_result.ysize()-1;y++)
	          for(int x=1;x<log_result.xsize()-1;x++) {
              float val = log_result.value(x,y,z,t);
              float ival = -val;

              if (
                (val > 0.0) && (
                  // 6 connectivity
                  (log_result.value(x-1,y,z,t) <= ival) ||  // 1
                  (log_result.value(x+1,y,z,t) <= ival) ||  // 2
                  (log_result.value(x,y-1,z,t) <= ival) ||  // 3
                  (log_result.value(x,y+1,z,t) <= ival) ||  // 4
                  (log_result.value(x,y,z-1,t) <= ival) ||  // 5
                  (log_result.value(x,y,z+1,t) <= ival)     // 6
        	  )
		){

                  zero_crossing_result(x,y,z,t) = 1.0;
	      	}

              if (
                (val < 0.0) && (
                  // 6 connectivity
                  (log_result.value(x-1,y,z,t) > ival) ||  // 1
                  (log_result.value(x+1,y,z,t) > ival) ||  // 2
                  (log_result.value(x,y-1,z,t) > ival) ||  // 3
                  (log_result.value(x,y+1,z,t) > ival) ||  // 4
                  (log_result.value(x,y,z-1,t) > ival) ||  // 5
                  (log_result.value(x,y,z+1,t) > ival)     // 6
                  )
		){

                  zero_crossing_result(x,y,z,t) = 1.0;
	      	}
            }

      return zero_crossing_result;
    }

   template <class T>
   volume<T> fixed_edge_detect(const volume<T>& source, float threshold,
			       bool twodimensional)
   {
     volume<T> result = source;
     int zsize = 3;
     if (twodimensional) zsize=1;

     volume<float> log_kern(3,3,zsize);
     log_kern = -1;
     log_kern(1,1,(zsize-1)/2) = 8;

     extrapolation oldex = source.getextrapolationmethod();
     source.setextrapolationmethod(mirror);
     result = convolve(source, log_kern);
     source.setextrapolationmethod(oldex);
     result.binarise(threshold);

     return result;
   }



  ///////////////////////////////////////////////////////////////////////////
  // EDGE STRENGTHEN (from avwmaths)
  template <class T>
  volume<T> edge_strengthen(const volume<T>& source)
  {
    float tmpf=2*sqrt(1/(pow((double)source.xdim(),2.0)) + 1/(pow((double)source.ydim(),2.0)) + 1/(pow((double)source.zdim(),2.0)));
       volume<T> result;
       result=source;
       if (source.zsize()>2)
       {
          for(int t=0;t<source.tsize();t++)
           for(int z=1;z<source.zsize()-1;z++)
             for(int y=1;y<source.ysize()-1;y++)
	       for(int x=1;x<source.xsize()-1;x++)
	       {
                 double temp1 = pow(double(source.value(x,y,z+1,t)-source.value(x,y,z-1,t)),2.0)/ pow((double)source.zdim(),2.0);
                 double temp2 = pow(double(source.value(x,y+1,z,t)-source.value(x,y-1,z,t)),2.0)/ pow((double)source.ydim(),2.0);
		 double temp3 = pow(double(source.value(x+1,y,z,t)-source.value(x-1,y,z,t)),2.0)/ pow((double)source.xdim(),2.0);
		 result(x,y,z,t)=(T)(sqrt(temp1+temp2+temp3)/tmpf);
	      }
       }
       else
       {
         for(int t=0;t<source.tsize();t++)
           for(int z=0;z<source.zsize();z++)
             for(int y=1;y<source.ysize()-1;y++)
	       for(int x=1;x<source.xsize()-1;x++)
	       {
                 double temp1=pow(double(source.value(x,y+1,z,t)-source.value(x,y-1,z,t)),2.0)/ pow((double)source.ydim(),2.0);
		 double temp2=pow(double(source.value(x+1,y,z,t)-source.value(x-1,y,z,t)),2.0)/ pow((double)source.xdim(),2.0);
		 result(x,y,z,t) = (T)(sqrt (temp1+temp2) / tmpf);
               }
       }
       return result;
  }


  ///////////////////////////////////////////////////////////////////////////
  // CONNECTED COMPONENTS

  // support functions for connected components

  int find_first_nonzero(const NEWMAT::Matrix& mat);

  void addpair2set(int x, int y, std::vector<int>& sx, std::vector<int>& sy);

 void relabel_components_uniquely(volume<int>& labelvol,
				   const std::vector<int>& equivlista,
				   const std::vector<int>& equivlistb, NEWMAT::ColumnVector& clustersizes);

  void relabel_components_uniquely(volume<int>& labelvol,
				   const std::vector<int>& equivlista,
				   const std::vector<int>& equivlistb);

  ////////////////////////////////////////////////////////////////////////////
  struct offset {
    int64_t x,y,z;
  offset(const int64_t ix,const int64_t iy,const int64_t iz) : x(ix),y(iy),z(iz) {}
  };

  std::vector<offset> backConnectivity(int nDirections);

  template <class T>
  void nonunique_component_labels(const volume<T>& vol,
				  volume<int>& labelvol,
				  std::vector<int>& equivlista,
				  std::vector<int>& equivlistb,
				  int numconnected)
    {
      copyconvert(vol,labelvol);
      labelvol = 0;

      int labelnum(0);
      equivlista.erase(equivlista.begin(),equivlista.end());
      equivlistb.erase(equivlistb.begin(),equivlistb.end());
      std::vector<offset> neighbours(backConnectivity(numconnected));
      for (int z=vol.minz(); z<=vol.maxz(); z++)
	for (int y=vol.miny(); y<=vol.maxy(); y++)
	  for (int x=vol.minx(); x<=vol.maxx(); x++) {
	    T val(vol(x,y,z));
	    if (val>0.5) {  // The eligibility test
	      int lval(labelvol(x,y,z));
	      for (std::vector<offset>::iterator it=neighbours.begin();it!=neighbours.end();it++) {
		int xnew(x+it->x),ynew(y+it->y),znew(z+it->z);
		if ( (xnew>=vol.minx()) && (ynew>=vol.miny()) && (znew>=vol.minz())
 		     && (MISCMATHS::round(vol(xnew,ynew,znew)))==val) {
		  // Binary relation
		  int lval2 = labelvol(xnew,ynew,znew);
		  if (lval != lval2) {
		    if (lval!=0)
		      addpair2set(lval2,lval,equivlista,equivlistb);
		    labelvol(x,y,z) = lval2;
		    lval = lval2;
		  }
		}
	      }
	      if (lval==0)
		labelvol(x,y,z) = ++labelnum;
	    }
	  }
    }

  template <class T>
  void nonunique_component_labels(const volume<T>& vol,
				  volume<int>& labelvol,
				  std::vector<int>& equivlista,
				  std::vector<int>& equivlistb,
				  bool (*binaryrelation)(T , T),
				  int numconnected)
    {
      copyconvert(vol,labelvol);
      labelvol = 0;

      int labelnum(0);
      equivlista.erase(equivlista.begin(),equivlista.end());
      equivlistb.erase(equivlistb.begin(),equivlistb.end());
      std::vector<offset> neighbours(backConnectivity(numconnected));
      for (int z=vol.minz(); z<=vol.maxz(); z++)
	for (int y=vol.miny(); y<=vol.maxy(); y++)
	  for (int x=vol.minx(); x<=vol.maxx(); x++) {
	    T val(vol(x,y,z));
	    if (val>0.5) {  // The eligibility test
	      int lval(labelvol(x,y,z));
	      for (std::vector<offset>::iterator it=neighbours.begin();it!=neighbours.end();it++) {
		int xnew(x+it->x),ynew(y+it->y),znew(z+it->z);
		if ((xnew>=vol.minx()) && (ynew>=vol.miny()) && (znew>=vol.minz())
		     && ((*binaryrelation)(vol(xnew,ynew,znew),val))) {
		  // Binary relation
		  int lval2 = labelvol(xnew,ynew,znew);
		  if (lval != lval2) {
		    if (lval!=0)
		      addpair2set(lval2,lval,equivlista,equivlistb);
		    labelvol(x,y,z) = lval2;
		    lval = lval2;
		  }
		}
	      }
	      if (lval==0)
		labelvol(x,y,z) = ++labelnum;
	    }
	  }
    }

  template <class T>
  volume<int> connected_components(const volume<T>& vol, NEWMAT::ColumnVector& clustersize, int numconnected)
    {
      volume<int> labelvol;
      copyconvert(vol,labelvol);
      std::vector<int> equivlista, equivlistb;
      nonunique_component_labels(vol,labelvol,
				 equivlista,equivlistb,numconnected);
      relabel_components_uniquely(labelvol,equivlista,equivlistb,clustersize);
      return labelvol;
    }


  template <class T>
  volume<int> connected_components(const volume<T>& vol,
                                   const volume<T>& mask,
                                   bool (*binaryrelation)(T , T), NEWMAT::ColumnVector& clustersize)
    {
      volume<int> labelvol;
      copyconvert(vol,labelvol);
      std::vector<int> equivlista, equivlistb;
      nonunique_component_labels(vol,mask,labelvol,equivlista,equivlistb,
				 binaryrelation);
      relabel_components_uniquely(labelvol,equivlista,equivlistb,clustersize);
      return labelvol;
    }


  template <class T>
  volume<int> connected_components(const volume<T>& vol, int numconnected){
    NEWMAT::ColumnVector clustersize;
    return connected_components(vol,clustersize,numconnected);
  }


  template <class T>
  volume<int> connected_components(const volume<T>& vol,const volume<T>& mask,
                                   bool (*binaryrelation)(T , T)){
    NEWMAT::ColumnVector clustersize;
    return connected_components(vol,mask,binaryrelation,clustersize);
  }


//////////////////////////////////// distancemap-related functions //////////////////////


class rowentry { public: int x; int y; int z; float d; } ;

bool rowentry_lessthan(const rowentry& r1, const rowentry& r2);

template <class T>
class distancemapper {
private:
  const volume<T> &bvol;
  const volume<T> &mask;
  const volume<T> &bvolneg;
  std::vector<rowentry> schedule;
  std::vector<offset> octantsign;
  bool dualmasks;
public:
  // basic constructor takes binaryvol (mask of valid values)
  //   + maskvol (non-zero at desired calculated points only)
  //   Dual binaryvol version (maskpos and maskneg) for signed distance to nearest mask
  distancemapper(const volume<T>& binarypos, const volume<T>& binaryneg, const volume<T>& maskvol);
  distancemapper(const volume<T>& binaryvol, const volume<T>& maskvol);
  ~distancemapper();
  const volume<float> distancemap();
  volume4D<float> sparseinterpolate(const volume4D<float>& values,
				    const std::string& interpmethod="general");
private:
  int setup_globals();
  int find_nearest(int x, int y, int z, int& x1, int& y1, int& z1,
                   bool findav, NEWMAT::ColumnVector& localav, const volume4D<float>& vals);
  int find_nearest(int x, int y, int z, int& x1, int& y1, int& z1);
  int create_distancemap(volume4D<float>& vout, const volume4D<float>& valim,
			  const std::string& interpmethod="none");
  int basic_create_distancemap(volume4D<float>& vout,
			       const volume4D<float>& valim,
			       const std::string& interpmethod);
};

template <class T>
distancemapper<T>::distancemapper(const volume<T>& binarypos, const volume<T>& binaryneg, const volume<T>& maskvol) :
  bvol(binarypos), mask(maskvol), bvolneg(binaryneg)
{
  if (!samesize(bvol,mask)) imthrow("Mask and image not the same size",20);
  if (!samesize(bvol,bvolneg)) imthrow("Two binary images are not the same size",20);
  dualmasks=true;
  setup_globals();
}


template <class T>
distancemapper<T>::distancemapper(const volume<T>& binaryvol, const volume<T>& maskvol) :
  bvol(binaryvol), mask(maskvol), bvolneg(binaryvol)
{
  if (!samesize(bvol,mask)) imthrow("Mask and image not the same size",20);
  dualmasks=false;
  setup_globals();
}

template <class T>
distancemapper<T>::~distancemapper()
{
  // destructors for the private members should do the job
}

template <class T>
int distancemapper<T>::setup_globals()
{
  // octantsign gives the 8 different octant sign combinations for coord
  //  offsets
  for (int p=-1; p<=1; p+=2)
    for (int q=-1; q<=1; q+=2)
      for (int r=-1; r<=1; r+=2)
	octantsign.push_back(offset(p,q,r));

  // construct list of displacements (in one octant) in ascending
  // order of distance
  for (int z=bvol.minz(); z<=bvol.maxz(); z++) {
    for (int y=bvol.miny(); y<=bvol.maxy(); y++) {
      for (int x=bvol.minx(); x<=bvol.maxx(); x++) {
	   rowentry newrow;
	   newrow.x=x;
	   newrow.y=y;
	   newrow.z=z;
	   float d2 = MISCMATHS::norm2sq(x*bvol.xdim(),y*bvol.ydim(),z*bvol.zdim());
	   newrow.d=d2;
	   schedule.push_back(newrow);
      }
    }
  }

  // sort schedule to get ascending d2
  sort(schedule.begin(),schedule.end(),NEWIMAGE::rowentry_lessthan);
  return 0;
}

// findav determines whether to do interpolation calculations or just
//  return the location only
template <class T>
int distancemapper<T>::find_nearest(int x, int y, int z, int& x1, int& y1, int& z1,
				 bool findav, NEWMAT::ColumnVector& localav,
				 const volume4D<float>& vals)
{
  float sumw=0.0, mindist=0, maxdist=0.0, weight;
  float minVoxSize = std::min(bvol.xdim(),std::min(bvol.ydim(),bvol.zdim()));
  NEWMAT::ColumnVector sumvw;
  if (findav) {
    localav.ReSize(vals.tsize());
    localav=0.0;
    sumvw.ReSize(vals.tsize());
    sumvw=0.0;
  }
  for (std::vector<rowentry>::iterator it=schedule.begin(); it!=schedule.end();it++) {
    for (std::vector<offset>::iterator oit=octantsign.begin(); oit!=octantsign.end();oit++) {
	x1=x+it->x*oit->x;
	y1=y+it->y*oit->y;
	z1=z+it->z*oit->z;
	if (bvol.in_bounds(x1,y1,z1)) {
	    bool primaryValid( bvol.value(x1,y1,z1) > 0.5);
	    if ( primaryValid || ( dualmasks && (bvolneg.value(x1,y1,z1)>0.5))) {
	      if ((!findav) || (dualmasks)) {
		return primaryValid ? 1 : -1; //Either primary or secondary is valid, so if primary is false...
	      } else if (mindist==0.0) {  // first time a point is encountered
		  mindist=it->d;
		  // select distance band to average over (farther -> more)
		  maxdist=MISCMATHS::Max(mindist+minVoxSize,mindist*1.5);
	      } else if (it->d>maxdist) {  // stop after maxdist reached
		localav=sumvw/sumw;
		return 1;
	      }
	      weight = it->d ? mindist/it->d : 1;
	      sumw += weight;
	      for (int t=0; t<vals.tsize(); t++)
		  sumvw(t+1) += weight * vals.value(x1,y1,z1,t);
	    }
	}
    }
  }
  // return most distant point (as last resort)
  if (findav) { if (sumw>0) { localav=sumvw/sumw; return 1; }}//should be (findav) and not not?
  return 0;  // not found any binary voxel (true for a zero image)
}


template <class T>
int distancemapper<T>::find_nearest(int x, int y, int z, int& x1, int& y1, int& z1)
{
  NEWMAT::ColumnVector dummy;
  volume4D<float> dummyvol;
  return this->find_nearest(x,y,z,x1,y1,z1,false,dummy,dummyvol);
}


// create the distance map as vout
// if return_distance is false then interpolate the value of the input volume
//  at the output location rather than store the distance to it
template <class T>
int distancemapper<T>::create_distancemap(volume4D<float>& vout,
					  const volume4D<float>& valim,
					  const std::string& interpmethod)
{
  if (interpmethod!="general")
    return this->basic_create_distancemap(vout,valim,interpmethod);
  // only get to here for sparseinterpolation
  float meanvoxsize = pow(valim.xdim()*valim.ydim()*valim.zdim(),1.0/3.0);
  int nsubsamp=0;
  if (meanvoxsize<4.0) { nsubsamp=1; }
  if (meanvoxsize<2.0) { nsubsamp=2; }
  if (meanvoxsize<1.0) { nsubsamp=3; }
  // for the straightforward case (>4mm resolution)
  if (nsubsamp==0) {
    return basic_create_distancemap(vout,valim,interpmethod);
  }
  // otherwise run subsamplings, sparseinterp and upsamplings
  //   this is to fill in the large gaps in the image, since this
  //   takes a *huge* amount of time otherwise
  // NB: mask in distancemapper is 1 where the new values are to go
  //     but in this function we use mask=1 where valid samples are
  std::vector<volume4D<float> > sampledData;
  std::vector<volume<float> > sampledMask;
  sampledData.resize(nsubsamp+1);
  sampledMask.resize(nsubsamp+1);
  sampledMask[0]=1.0f - mask;  // sampledMask is now 1 for all valid sample points
  sampledData[0]=valim*sampledMask[0];
  for (int n=1; n<=nsubsamp; n++) {
    sampledMask[n]=subsample_by_2(sampledMask[n-1],false);
    sampledData[n].reinitialize(sampledMask[n].xsize(),sampledMask[n].ysize(),sampledMask[n].zsize(),valim.tsize());
    sampledData[n].copyproperties(sampledMask[n]);
    for (int t=0; t<valim.tsize(); t++)
      (sampledData[n])[t]=divide(subsample_by_2((sampledData[n-1])[t],false),sampledMask[n],sampledMask[n]);

    sampledMask[n].binarise(1e-4);   // include any voxel with any partial overlap
  }
  // run the sparseinterpolate function (at 8mm-ish resolution)
  //   - not the method in this object, but a new one
  sampledData[nsubsamp]=NEWIMAGE::sparseinterpolate(sampledData[nsubsamp],sampledMask[nsubsamp]);
  // dilate and invert mask
  volume<float> kernel(3,3,3);
  kernel=1.0f;
  sampledMask[nsubsamp]=morphfilter(sampledMask[nsubsamp],kernel,"dilate");
  sampledMask[nsubsamp]= 1.0f - sampledMask[nsubsamp];  // now the mask is 1 for all the more distant "gaps"
  // upsample mask and spareinterpolated result
  volume<float> tmp;
  for (int n=nsubsamp; n>0; n--) {
    for (int t=0; t<valim.tsize(); t++) {
      ShadowVolume<float> shadow(sampledData[n-1][t]);
      upsample_by_2(shadow,sampledData[n][t],false);
    }
    ShadowVolume<float> shadow(sampledMask[n-1]);
    upsample_by_2(shadow,sampledMask[n],false);
  }
  sampledMask[0].binarise(0.5f);
  // only keep results in zero part of mask as well as the
  //  original points
  volume4D<float> vres(valim*(1.0f-mask) + sampledData[0]*sampledMask[0]);
  sampledMask[0] = mask - sampledMask[0];   // only calculate new points in the "gap"
  // now run basic_create_distancemap at full resolution, but with the
  //  new mask extra filled-in areas from above
  volume<float> invMask(1.0f - sampledMask[0]);
  distancemapper<float> newdmapper(invMask,sampledMask[0]);
  return newdmapper.basic_create_distancemap(vout,vres,interpmethod);
}

// The following function creates the distancemap or interpolated image
// from valim (previously the mask of valid voxels = binaryvol and the
// mask of voxels to calculate = maskvol, must have been specified)
template <class T>
int distancemapper<T>::basic_create_distancemap(volume4D<float>& vout,
						const volume4D<float>& valim,
						const std::string& interpmethod)
{
  int x1, y1, z1;
  NEWMAT::ColumnVector localav;
  int interp=0, distsign=1;
  if ((interpmethod=="nn") || (interpmethod=="nearestneighbour")) interp=1;
  if (interpmethod=="general") interp=2;
  if (interp>0) { vout = valim; } else { vout = bvol; vout *= 0.0f; }
  if ((interp>0) && (!samesize(bvol,valim.constSubVolume(0))))
    {  print_volume_info(bvol,std::string("bvol"),std::cerr); print_volume_info(mask,std::string("mask"),std::cerr); print_volume_info(valim.constSubVolume(0),std::string("valim"),std::cerr);
       imthrow("Binary image and interpolant not the same size",21); }
  for (int z=0; z<vout.zsize(); z++) {
    for (int y=0; y<vout.ysize(); y++) {
      for (int x=0; x<vout.xsize(); x++) {
	if (mask(x,y,z)>((T) 0.5)) {
	  distsign=find_nearest(x,y,z,x1,y1,z1,interp>=2,localav,valim);
	  switch (interp) {
	  case 2:
	    for (int t=0;t<valim.tsize();t++) { vout(x,y,z,t)=localav(t+1); }
	    break;
	  case 1:
	    for (int t=0;t<valim.tsize();t++) { vout(x,y,z,t)=valim(x1,y1,z1,t); }
	    break;
	  case 0:
	  default:
	    vout(x,y,z,0)=distsign*sqrt(MISCMATHS::norm2sq((x1-x)*bvol.xdim(),
				       (y1-y)*bvol.ydim(),(z1-z)*bvol.zdim()));
	  }
	}
      }
    }
  }
  return 0;
}


template <class T>
const volume<float> distancemapper<T>::distancemap()
{
  volume4D<float> dmap;
  create_distancemap(dmap,dmap,"none");
  return dmap[0];
}

template <class T>
volume4D<float> distancemapper<T>::sparseinterpolate(const volume4D<float>& values,
						     const std::string& interpmethod)
{
  volume4D<float> vout;
  create_distancemap(vout,values,interpmethod);
  return vout;
}


template <class T>
volume<float> distancemap(const volume<T>& binaryvol)
{
  volume<T> mask;
  mask = ((T) 1) - binarise(binaryvol,((T) 0.5));
  return distancemap(binaryvol,mask);
}

template <class T>
volume<float> distancemap(const volume<T>& binaryvol, const volume<T>& mask)
{
  distancemapper<T> dmapper(binaryvol,mask);
  return dmapper.distancemap();
}

template <class T>
volume<float> distancemap(const volume<T>& binarypos, const volume<T>& binaryneg, const volume<T>& maskvol)
{
  distancemapper<T> dmapper(binarypos,binaryneg,maskvol);
  return dmapper.distancemap();
}

template <class T>
volume4D<float> sparseinterpolate(const volume4D<T>& sparsesamps,
				  const volume<T>& mask,
				  const std::string& interpmethod)
{
  // can have "general" or "nearestneighbour" (or "nn") for interpmethod
  volume<T> invmask;
  invmask=((T) 1) - mask;
  distancemapper<T> dmapper(mask,invmask);
  return dmapper.sparseinterpolate(sparsesamps,interpmethod);
}


//////////////////////////////////// TFCE-related functions //////////////////////

template <class T>
void tfce_orig_slow(volume<T>& VolIntn, float H, float E, int NumConn, float minT, float deltaT)
{
  float maxval=VolIntn.max();
  volume<float> clusterenhance;
  copyconvert(VolIntn,clusterenhance);
  clusterenhance=0;

  if (deltaT==0)
    deltaT = (maxval - minT)/100.0;   // this needs fixing!!

  for (float thresh=minT+deltaT; thresh<=maxval; thresh+=deltaT)
    {
      volume<float> clusters;
      copyconvert(VolIntn,clusters);
      clusters.binarise(thresh);

      NEWMAT::ColumnVector clustersizes;
      volume<int>tmpvol=connected_components(clusters,clustersizes,NumConn);
      clustersizes = pow(clustersizes,E) * pow(thresh,H);
      for(int z=0;z<VolIntn.zsize();z++)
	for(int y=0;y<VolIntn.ysize();y++)
	  for(int x=0;x<VolIntn.xsize();x++)
	    if (tmpvol.value(x,y,z)>0)
	      clusterenhance.value(x,y,z) += clustersizes(tmpvol.value(x,y,z));
    }
  copyconvert(clusterenhance,VolIntn);
  return;
}

class VecSort{
 public:
  int Sx, Sy, Sz, Sl;
  double Sv;
  bool operator<(const VecSort& other) const{
    return Sv < other.Sv;
  }
};

template <class T>
void tfce(volume<T>& data, float H, float E, int NumConn, float minT, float deltaT)
{
  volume<int> VolLabl; copyconvert(data, VolLabl);
  volume<float> VolEnhn; copyconvert(data, VolEnhn); VolEnhn=0;
  bool doIT=false;
  const int INIT=-1, MASK=-2;
  int curlab=0;
  int pX, pY, pZ, qX, qY, qZ, rX, rY, rZ;
  int FldCntr=0, FldCntri=0, tfceCntr=0, xFC = 0;
  int minX=1, minY=1, minZ=1, maxX=data.maxx()-1, maxY=data.maxy()-1, maxZ=data.maxz()-1;
  int sizeC=maxX*maxY*maxZ;
  int counter=0, edsta[27];
  float maxT=data.max();
  if(deltaT==0)
    deltaT=maxT/100.0;
  if(deltaT<=0)
    throw std::runtime_error("Error: tfce requires a positive deltaT input.");
  if ( data.xsize() < 3 || data.ysize() < 3 || data.zsize() < 3 )
    throw std::runtime_error("Error: tfce currently requires an input with at least 3 voxels extent into each dimension.");
  if( data.max()/deltaT > 10000 )
    std::cout << "Warning: tfce has detected a large number of integral steps. This operation may require a great deal of time to complete." << std::endl;
  std::vector<VecSort> VecSortI(sizeC);
  std::queue<int> Qx, Qy, Qz;
  for(int z0=-1; z0<=1; z0++)
    for(int y0=-1; y0<=1; y0++)
      for(int x0=-1; x0<=1; x0++){
	edsta[counter++] = (x0*x0+y0*y0+z0*z0);
	if (edsta[counter-1]<2 && edsta[counter-1]!=0) edsta[counter-1]=6;
	if (edsta[counter-1]<3 && edsta[counter-1]!=0) edsta[counter-1]=18;
	if (edsta[counter-1]<4 && edsta[counter-1]!=0) edsta[counter-1]=26;
	if (edsta[counter-1]==0) edsta[counter-1]=100;
      }
  FldCntr=0;
  for(int z=minZ; z<=maxZ; z++)
    for(int y=minY; y<=maxY; y++)
      for(int x=minX; x<=maxX; x++) {
	float iVal=data.value(x,y,z);
	if( iVal > minT ) {
	  VecSortI[FldCntr].Sx=x; VecSortI[FldCntr].Sy=y; VecSortI[FldCntr].Sz=z;
	  VecSortI[FldCntr++].Sv=iVal;
	}
      }
  sizeC=FldCntr;
  VecSortI.resize(sizeC);
  sort(VecSortI.begin(), VecSortI.end());
  for(float curThr=minT; curThr<(maxT+deltaT);curThr+=deltaT){
    VolLabl = INIT; FldCntr = xFC; curlab = 0;
    while( (VecSortI[FldCntr].Sv<=curThr) && (FldCntr<sizeC) ){
      VecSortI[FldCntr].Sl=0; VecSortI[FldCntr++].Sv=0;
    }
    xFC=FldCntr;
    while( VecSortI[FldCntr].Sv>curThr  && (FldCntr<sizeC) ){
      pX=VecSortI[FldCntr].Sx; pY=VecSortI[FldCntr].Sy; pZ=VecSortI[FldCntr++].Sz;
      VolLabl.value(pX,pY,pZ)=MASK;
    }
    for(FldCntri=xFC; FldCntri<FldCntr; FldCntri++){
      pX=VecSortI[FldCntri].Sx; pY=VecSortI[FldCntri].Sy; pZ=VecSortI[FldCntri].Sz;
      if(VolLabl.value(pX, pY, pZ)==MASK){//sI
	curlab+=1;
	Qx.push(pX); Qy.push(pY); Qz.push(pZ);
	VolLabl.value(pX, pY, pZ)=curlab;
	while(!Qx.empty()){//sW
	  qX=Qx.front(); qY=Qy.front(); qZ=Qz.front();
	  Qx.pop(); Qy.pop(); Qz.pop();
	  counter=0;
	  for(int z0=-1; z0<=1; z0++)
	    for(int y0=-1; y0<=1; y0++)
	      for(int x0=-1; x0<=1; x0++){
		rX=qX+x0; rY=qY+y0; rZ=qZ+z0;
		doIT =(NumConn>=edsta[counter++]);
		if(doIT && (VolLabl.value(rX, rY, rZ)==MASK)){
		  Qx.push(rX); Qy.push(rY);Qz.push(rZ);
		  VolLabl.value(rX, rY, rZ)=curlab;
		}
	      }
	}//eW
      }//eI
    }
    NEWMAT::ColumnVector ClusterSizes(curlab), ClusterSizesI(curlab);
    ClusterSizes=0;
    for(tfceCntr=xFC; tfceCntr<sizeC; tfceCntr++){
      VecSortI[tfceCntr].Sl=VolLabl.value(VecSortI[tfceCntr].Sx, VecSortI[tfceCntr].Sy, VecSortI[tfceCntr].Sz);
      if ( VecSortI[tfceCntr].Sl>0 )
	ClusterSizes(VecSortI[tfceCntr].Sl)+=1;
    }
    float HH=pow(curThr, H);
    ClusterSizesI=pow(ClusterSizes, E)*HH;
    for(tfceCntr=xFC; tfceCntr<sizeC; tfceCntr++){
      if ( VecSortI[tfceCntr].Sl>0 ){
	VolEnhn.value(VecSortI[tfceCntr].Sx, VecSortI[tfceCntr].Sy, VecSortI[tfceCntr].Sz) += ClusterSizesI(VecSortI[tfceCntr].Sl);
      }
    }
  }//end curThr
  copyconvert(VolEnhn,data);
  return;
}

template <class T>
void tfce_support(volume<T>& VolIntn, float H, float E, int NumConn, float minT, float deltaT, int Xoi, int Yoi, int Zoi, float threshTFCE)
{
  volume<float> VolEnhn; copyconvert(VolIntn, VolEnhn);
  volume<float> VolTemp; copyconvert(VolIntn, VolTemp);
  int minX=0, minY=0, minZ=0, maxX=VolIntn.maxx(), maxY=VolIntn.maxy(), maxZ=VolIntn.maxz();
  float maxT = VolIntn.value(Xoi,Yoi,Zoi);
  int thrCntr = 0; int thrNum = 0;
  if(deltaT==0){
    deltaT = maxT/100;
    thrNum = 101;
  }
  else
    thrNum = int(maxT/deltaT)+1;
  NEWMAT::ColumnVector Clusters(thrNum), Thresholds(thrNum), ClusterSizes;
  for(float thresh=minT; thresh<maxT; thresh+=deltaT){
    copyconvert(VolEnhn, VolTemp);
    VolTemp.binarise(thresh);
    volume<int> VolLabl = connected_components(VolTemp, ClusterSizes, NumConn);
    Clusters(thrCntr+1) = ClusterSizes(VolLabl.value(Xoi, Yoi, Zoi));
    Thresholds(thrCntr+1) = thresh;
    for(int z=minZ; z<=maxZ; z++)
      for(int y=minY; y<=maxY; y++)
	for(int x=minX; x<=maxX; x++)
	  if( VolLabl.value(x,y,z) != VolLabl.value(Xoi,Yoi,Zoi) )
	    VolEnhn.value(x,y,z) = 0;
    thrCntr++;
  }
  float  summ=0, thre=0;
  for (int i=0; i<thrCntr; i++){
    summ+=std::pow((float)Clusters(thrCntr-i), E)*std::pow((float)Thresholds(thrCntr-i), H);
    thre = Thresholds(thrCntr-i);
    if(summ>=threshTFCE)
      break;
  }
  if(summ<threshTFCE)
  std::cout<<"it doesn't reach to specified threshold"<<std::endl;
  copyconvert(VolIntn, VolEnhn);
  copyconvert(VolIntn, VolTemp);     VolTemp.binarise(thre);
  volume<int> VolLabl = connected_components(VolTemp, ClusterSizes, NumConn);
  for(int z=minZ; z<=maxZ; z++)
    for(int y=minY; y<=maxY; y++)
      for(int x=minX; x<=maxX; x++)
	if( VolLabl.value(x,y,z)!=VolLabl(Xoi, Yoi, Zoi) )
	  VolEnhn.value(x,y,z) = 0;
  copyconvert(VolEnhn,VolIntn);
  return;
}




////////////////////////////////////////////////////////////////////////////

}


#endif
