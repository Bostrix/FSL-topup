// Definitions for class splinefield
//
// splinefield.cpp
//
// Jesper Andersson and Matthew Webster, FMRIB Image Analysis Group
//
// Copyright (C) 2007-2014 University of Oxford
//
//     CCOPYRIGHT
//

#include <ctime>
#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <thread>
#include <memory>
#include "armawrap/newmat.h"
#include "utils/threading.h"
#include "miscmaths/bfmatrix.h"
#include "fsl_splines.h"
#include "splinefield.h"

using namespace std;
using namespace NEWMAT;
using namespace MISCMATHS;

namespace BASISFIELD {

// Constructor, assignement and destructor

// Plain vanilla constructors
splinefield::splinefield(const std::vector<unsigned int>& psz,
			 const std::vector<double>&       pvxs,
			 const std::vector<unsigned int>& ksp,
			 int                              order,
			 Utilities::NoOfThreads           nthr)
: basisfield(psz,pvxs), _sp(order,ksp), _dsp(3), _nthr(nthr._n)
{
  if (order < 2 || order > 3) {throw BasisfieldException("splinefield::splinefield: Only quadratic and cubic splines implemented yet");}
  if (ksp.size() != NDim()) {throw BasisfieldException("splinefield::splinefield: Dimensionality mismatch");}
  /*
  if (pksp[0]<0 || pksp[0]>FieldSz_x() || (NDim()>1 && (pksp[1]<0 || pksp[1]>FieldSz_y())) || (NDim()>2 && (pksp[2]<0 || pksp[2]>FieldSz_z()))) {
    throw BasisfieldException("splinefield::splinefield: Invalid knot-spacing");
  }
  */
  for (unsigned int i=0; i<3; i++) {
    std::vector<unsigned int>  deriv(3,0);
    deriv[i] = 1;
    _dsp[i] = std::shared_ptr<Spline3D<double> >(new Spline3D<double>(order,ksp,deriv));
  }

  std::shared_ptr<NEWMAT::ColumnVector>  lcoef(new NEWMAT::ColumnVector(CoefSz()));
  *lcoef = 0.0;
  set_coef_ptr(lcoef);
}

// Copy constructor
splinefield::splinefield(const splinefield& inf) : basisfield(inf), _sp(inf._sp), _nthr(inf._nthr), _dsp(inf._dsp)
{
  // basisfield::assign(inf);
  // splinefield::assign(inf);
}

void splinefield::assign_splinefield(const splinefield& inf)
{
  _sp = inf._sp;
  for (unsigned int i=0; i<3; i++) {
    _dsp[i] = std::shared_ptr<Spline3D<double> >(new Spline3D<double>(*(inf._dsp[i])));
  }
  _nthr = inf._nthr;
}

splinefield& splinefield::operator=(const splinefield& inf)
{
  if (&inf == this) {return(*this);} // Detect self

  basisfield::assign_basisfield(inf);   // Assign common part
  splinefield::assign_splinefield(inf);  // Assign splinefield specific bits

  return(*this);
}

// General utility functions

// Functions that actually do some work

double splinefield::Peek(double x, double y, double z, FieldIndex fi) const
{
  const Spline3D<double>   *sp_ptr = 0;
  if (fi == FIELD) sp_ptr = &_sp;
  else {
    std::vector<unsigned int>   deriv(3,0);
    switch (fi) {
    case DFDX:
      deriv[0] = 1;
      break;
    case DFDY:
      deriv[1] = 1;
      break;
    case DFDZ:
      deriv[2] = 1;
      break;
    default:
      throw BasisfieldException("Peek: Invalid FieldIndex value");
    }
    sp_ptr = new Spline3D<double>(_sp.Order(),_sp.KnotSpacing(),deriv);
  }

  std::vector<double>         vox(3);
  vox[0]=x; vox[1]=y; vox[2]=z;
  std::vector<unsigned int>   coefsz(3);
  coefsz[0] = CoefSz_x(); coefsz[1] = CoefSz_y(); coefsz[2] = CoefSz_z();
  std::vector<unsigned int>   first_cindx(3);
  std::vector<unsigned int>   last_cindx(3);
  sp_ptr->RangeOfSplines(vox,coefsz,first_cindx,last_cindx);

  double rval = 0.0;
  std::vector<unsigned int>  cindx(3,0);
  for (cindx[2]=first_cindx[2]; cindx[2]<last_cindx[2]; cindx[2]++) {
    for (cindx[1]=first_cindx[1]; cindx[1]<last_cindx[1]; cindx[1]++) {
      for (cindx[0]=first_cindx[0]; cindx[0]<last_cindx[0]; cindx[0]++) {
        rval += GetCoef(cindx[0],cindx[1],cindx[2]) * sp_ptr->SplineValueAtVoxel(vox,cindx);
      }
    }
  }

  if (fi != FIELD) delete sp_ptr;

  return(rval);
}

void splinefield::Update(FieldIndex fi)
{
  if (fi>int(NDim())) {throw BasisfieldException("splinefield::Update: Cannot take derivative in singleton direction");}

  if (UpToDate(fi)) {return;} // Field already fine.

  const std::shared_ptr<NEWMAT::ColumnVector> lcoef = GetCoef();
  if (!lcoef) {throw BasisfieldException("splinefield::Update: No coefficients set yet");}

  // Get spline kernel
  std::vector<unsigned int> deriv(3,0);
  if (fi) deriv[fi-1] = 1;
  Spline3D<double> spline(_sp.Order(),_sp.KnotSpacing(),deriv);

  std::vector<unsigned int> coefsz(3,0);
  coefsz[0] = CoefSz_x(); coefsz[1] = CoefSz_y(); coefsz[2] = CoefSz_z();
  std::vector<unsigned int> fieldsz(3,0);
  fieldsz[0] = FieldSz_x(); fieldsz[1] = FieldSz_y(); fieldsz[2] = FieldSz_z();
  std::shared_ptr<NEWMAT::ColumnVector>  fptr = get_ptr(fi);

  get_field(spline,static_cast<double *>(lcoef->Store()),coefsz,fieldsz,static_cast<double *>(fptr->Store()));
  set_update_flag(true,fi);
}

void splinefield::Update(const std::vector<FieldIndex>& fiv)
{
  if (std::any_of(fiv.cbegin(),fiv.cend(),[this](const FieldIndex& fi){ return(fi>int(this->NDim())); })) { throw BasisfieldException("splinefield::Update: Cannot take derivative in singleton direction"); }

  if (std::all_of(fiv.cbegin(),fiv.cend(),[this](const FieldIndex& fi){ return(this->UpToDate(fi)); })) { return; }

  const std::shared_ptr<NEWMAT::ColumnVector> lcoef = GetCoef();
  if (!lcoef) {throw BasisfieldException("splinefield::Update: No coefficients set yet");}

  // Get spline kernels and pointers for the indices that are not up to date
  std::vector<Spline3D<double> > splines;
  std::vector<double *> fptrs;
  for (unsigned int i=0; i<fiv.size(); i++) {
    if (!UpToDate(fiv[i])) {
      std::vector<unsigned int> deriv(3,0);
      if (fiv[i] != 0) deriv[fiv[i]-1] = 1;
      splines.push_back(Spline3D<double>(_sp.Order(),_sp.KnotSpacing(),deriv));
      std::shared_ptr<NEWMAT::ColumnVector>  fptr = get_ptr(fiv[i]);
      fptrs.push_back(static_cast<double *>(fptr->Store()));
    }
  }

  // Some house keeping.  
  std::vector<unsigned int> coefsz(3,0);
  coefsz[0] = CoefSz_x(); coefsz[1] = CoefSz_y(); coefsz[2] = CoefSz_z();
  std::vector<unsigned int> fieldsz(3,0);
  fieldsz[0] = FieldSz_x(); fieldsz[1] = FieldSz_y(); fieldsz[2] = FieldSz_z();

  // Get the fields
  get_fields(splines,static_cast<double *>(lcoef->Store()),coefsz,fieldsz,fptrs);
  std::for_each(fiv.cbegin(),fiv.cend(),[this](const FieldIndex& fi){ this->set_update_flag(true,fi); });
}


NEWMAT::ReturnMatrix splinefield::Jte(const NEWIMAGE::volume<float>&  ima1,
                                      const NEWIMAGE::volume<float>&  ima2,
                                      const NEWIMAGE::volume<char>    *mask) const
{
  std::vector<unsigned int> deriv(NDim(),0);
  NEWMAT::ColumnVector      tmp = Jte(deriv,ima1,ima2,mask);
  tmp.Release();
  return(tmp);
}

NEWMAT::ReturnMatrix splinefield::Jte(const std::vector<unsigned int>&  deriv,
                                      const NEWIMAGE::volume<float>&    ima1,
                                      const NEWIMAGE::volume<float>&    ima2,
                                      const NEWIMAGE::volume<char>      *mask) const
{
  if (deriv.size() != 3) throw BasisfieldException("splinefield::Jte: Wrong size deriv vector");
  if (!samesize(ima1,ima2,3,true) || (mask && !samesize(ima1,*mask,3))) {
    throw BasisfieldException("splinefield::Jte: Image dimensionality mismatch");
  }
  if (static_cast<unsigned int>(ima1.xsize()) != FieldSz_x() ||
      static_cast<unsigned int>(ima1.ysize()) != FieldSz_y() ||
      static_cast<unsigned int>(ima1.zsize()) != FieldSz_z()) {
    throw BasisfieldException("splinefield::Jte: Image-Field dimensionality mismatch");
  }
  float *prodima = new float[FieldSz()];
  hadamard(ima1,ima2,mask,prodima);

  Spline3D<double>           spline(_sp.Order(),_sp.KnotSpacing(),deriv);
  std::vector<unsigned int>  coefsz(3,0);
  coefsz[0] = CoefSz_x(); coefsz[1] = CoefSz_y(); coefsz[2] = CoefSz_z();
  std::vector<unsigned int>  imasz(3,0);
  imasz[0] = FieldSz_x(); imasz[1] = FieldSz_y(); imasz[2] = FieldSz_z();
  NEWMAT::ColumnVector       ovec(CoefSz());

  if (_nthr==1) get_jte(spline,coefsz,prodima,imasz,static_cast<double *>(ovec.Store()));
  else get_jte_para(spline,coefsz,prodima,imasz,static_cast<double *>(ovec.Store()),_nthr);

  delete[] prodima;
  ovec.Release();
  return(ovec);
}

NEWMAT::ReturnMatrix splinefield::Jte(const NEWIMAGE::volume<float>&    ima,
                                      const NEWIMAGE::volume<char>      *mask) const
{
  std::vector<unsigned int> deriv(NDim(),0);
  NEWMAT::ColumnVector tmp = Jte(deriv,ima,mask);
  return(tmp);
}

NEWMAT::ReturnMatrix splinefield::Jte(const std::vector<unsigned int>&  deriv,
                                      const NEWIMAGE::volume<float>&    ima,
                                      const NEWIMAGE::volume<char>      *mask) const
{
  if (deriv.size() != 3) throw BasisfieldException("splinefield::Jte: Wrong size if deriv vector");
  if (mask && !samesize(ima,*mask,3)) {
    throw BasisfieldException("splinefield::Jte: Image-Mask dimensionality mismatch");
  }
  if (static_cast<unsigned int>(ima.xsize()) != FieldSz_x() ||
      static_cast<unsigned int>(ima.ysize()) != FieldSz_y() ||
      static_cast<unsigned int>(ima.zsize()) != FieldSz_z()) {
    throw BasisfieldException("splinefield::Jte: Image-Field dimensionality mismatch");
  }

  float *fima = new float[FieldSz()];
  float *fiptr = fima;
  if (mask) {
    NEWIMAGE::volume<char>::fast_const_iterator itm = mask->fbegin();
    for (NEWIMAGE::volume<float>::fast_const_iterator it=ima.fbegin(), it_end=ima.fend(); it!=it_end; ++it, ++itm, ++fiptr) {
      if (*itm) *fiptr = *it;
      else *fiptr = 0.0;
    }
  }
  else {
    for (NEWIMAGE::volume<float>::fast_const_iterator it=ima.fbegin(), it_end=ima.fend(); it!=it_end; ++it, ++fiptr) *fiptr = *it;
  }

  Spline3D<double>           spline(_sp.Order(),_sp.KnotSpacing(),deriv);
  std::vector<unsigned int>  coefsz(3,0);
  coefsz[0] = CoefSz_x(); coefsz[1] = CoefSz_y(); coefsz[2] = CoefSz_z();
  std::vector<unsigned int>  imasz(3,0);
  imasz[0] = FieldSz_x(); imasz[1] = FieldSz_y(); imasz[2] = FieldSz_z();
  NEWMAT::ColumnVector       ovec(CoefSz());

  if (_nthr==1) get_jte(spline,coefsz,fima,imasz,static_cast<double *>(ovec.Store()));
  else get_jte_para(spline,coefsz,fima,imasz,static_cast<double *>(ovec.Store()),_nthr);

  delete[] fima;
  ovec.Release();
  return(ovec);
}

std::shared_ptr<SpMat<float> > splinefield::J() const
{
  std::vector<unsigned int> deriv(3,0);
  NEWIMAGE::volume<float>   ones(static_cast<int>(FieldSz_x()),static_cast<int>(FieldSz_y()),static_cast<int>(FieldSz_z())); ones = 1.0;
  std::shared_ptr<SpMat<float> > tmp = J(deriv,ones,NULL);
  return(tmp);
}

std::shared_ptr<SpMat<float> > splinefield::J(const std::vector<unsigned int>&  deriv) const
{
  NEWIMAGE::volume<float>   ones(static_cast<int>(FieldSz_x()),static_cast<int>(FieldSz_y()),static_cast<int>(FieldSz_z())); ones = 1.0;
  std::shared_ptr<SpMat<float> > tmp = J(deriv,ones,NULL);
  return(tmp);
}

std::shared_ptr<SpMat<float> > splinefield::J(const NEWIMAGE::volume<float>&    ima,
					      const NEWIMAGE::volume<char>      *mask) const
{
  std::vector<unsigned int>  deriv(3,0);
  std::shared_ptr<SpMat<float> >  tmp = J(deriv,ima,mask);
  return(tmp);
}

std::shared_ptr<SpMat<float> > splinefield::J(const std::vector<unsigned int>&  deriv,
					      const NEWIMAGE::volume<float>&    ima,
					      const NEWIMAGE::volume<char>      *mask) const
{
  if (deriv.size() != 3) throw BasisfieldException("splinefield::J: Wrong size derivative vector");
  if (mask != NULL && !samesize(ima,*mask,3)) throw BasisfieldException("splinefield::J: Mismatch between image and mask");
  if (FieldSz_x() != static_cast<unsigned int>(ima.xsize()) ||
      FieldSz_y() != static_cast<unsigned int>(ima.ysize()) ||
      FieldSz_z() != static_cast<unsigned int>(ima.zsize())) throw BasisfieldException("splinefield::J: Mismatch between image and field");

  float *pima = new float[FieldSz()];
  copy_and_mask_ima(ima,mask,pima);
  std::vector<unsigned int>  isz(3,0);
  isz[0] = FieldSz_x(); isz[1] = FieldSz_y(); isz[2] = FieldSz_z();
  std::vector<unsigned int>  csz(3,0);
  csz[0] = CoefSz_x(); csz[1] = CoefSz_y(); csz[2] = CoefSz_z();

  std::shared_ptr<SpMat<float> > tmp;
  if (deriv[0]==0 && deriv[1]==0 && deriv[2]==0) {
    tmp = get_j(_sp,csz,pima,isz);
  }
  else {
    Spline3D<double>   sp(_sp.Order(),_sp.KnotSpacing(),deriv);
    tmp = get_j(sp,csz,pima,isz);
  }

  delete[] pima;
  return(tmp);

}

std::shared_ptr<BFMatrix> splinefield::JtJ(const NEWIMAGE::volume<float>&     ima,
                                           const NEWIMAGE::volume<char>       *mask,
                                           MISCMATHS::BFMatrixPrecisionType   prec)
const
{
  std::vector<unsigned int>  deriv(3,0);
  std::shared_ptr<BFMatrix>  tmp = JtJ(deriv,ima,ima,mask,prec);
  return(tmp);
}

std::shared_ptr<BFMatrix> splinefield::JtJ(const NEWIMAGE::volume<float>&     ima1,
                                           const NEWIMAGE::volume<float>&     ima2,
                                           const NEWIMAGE::volume<char>       *mask,
                                           MISCMATHS::BFMatrixPrecisionType   prec)
const
{
  std::vector<unsigned int>  deriv(3,0);
  std::shared_ptr<BFMatrix>  tmp = JtJ(deriv,ima1,ima2,mask,prec);
  return(tmp);
}

std::shared_ptr<BFMatrix> splinefield::JtJ(const std::vector<unsigned int>&   deriv,
                                           const NEWIMAGE::volume<float>&     ima,
                                           const NEWIMAGE::volume<char>       *mask,
                                           MISCMATHS::BFMatrixPrecisionType   prec)
const
{
  std::shared_ptr<BFMatrix>  tmp = JtJ(deriv,ima,ima,mask,prec);
  return(tmp);
}

std::shared_ptr<BFMatrix> splinefield::JtJ(const std::vector<unsigned int>&   deriv,
                                           const NEWIMAGE::volume<float>&     ima1,
                                           const NEWIMAGE::volume<float>&     ima2,
                                           const NEWIMAGE::volume<char>       *mask,
                                           MISCMATHS::BFMatrixPrecisionType   prec)
const
{
  if (deriv.size() != 3) throw BasisfieldException("splinefield::JtJ: Wrong size derivative vector");
  if (!samesize(ima1,ima2,3,true)) throw BasisfieldException("splinefield::JtJ: Image dimension mismatch");
  if (mask && !samesize(ima1,*mask,3)) throw BasisfieldException("splinefield::JtJ: Mismatch between image and mask");
  if (FieldSz_x() != static_cast<unsigned int>(ima1.xsize()) ||
      FieldSz_y() != static_cast<unsigned int>(ima1.ysize()) ||
      FieldSz_z() != static_cast<unsigned int>(ima1.zsize())) throw BasisfieldException("splinefield::JtJ: Mismatch between image and field");

  float *prodima = new float[FieldSz()];
  hadamard(ima1,ima2,mask,prodima);
  std::vector<unsigned int>  isz(3,0);
  isz[0] = FieldSz_x(); isz[1] = FieldSz_y(); isz[2] = FieldSz_z();
  std::vector<unsigned int>  csz(3,0);
  csz[0] = CoefSz_x(); csz[1] = CoefSz_y(); csz[2] = CoefSz_z();

  std::shared_ptr<BFMatrix> tmp;
  if (deriv[0]==0 && deriv[1]==0 && deriv[2]==0) {
    if (_nthr==1) tmp = make_fully_symmetric_jtj(_sp,csz,prodima,isz,prec);
    else tmp = make_fully_symmetric_jtj_para(_sp,csz,prodima,isz,prec,_nthr);
  }
  else {
    Spline3D<double>   sp(_sp.Order(),_sp.KnotSpacing(),deriv);
    if (_nthr==1) tmp = make_fully_symmetric_jtj(sp,csz,prodima,isz,prec);
    else tmp = make_fully_symmetric_jtj_para(sp,csz,prodima,isz,prec,_nthr);
  }

  delete[] prodima;
  return(tmp);
}

std::shared_ptr<BFMatrix> splinefield::JtJ(const std::vector<unsigned int>&   deriv1,
                                             const NEWIMAGE::volume<float>&     ima1,
                                             const std::vector<unsigned int>&   deriv2,
                                             const NEWIMAGE::volume<float>&     ima2,
                                             const NEWIMAGE::volume<char>       *mask,
                                             MISCMATHS::BFMatrixPrecisionType   prec)
const
{
  if (deriv1.size() != 3 || deriv2.size() != 3) throw BasisfieldException("splinefield::JtJ: Wrong size derivative vector");

  std::shared_ptr<BFMatrix>  tmp;
  if (deriv1 == deriv2) tmp = JtJ(deriv1,ima1,ima2,mask,prec);
  else {
    if (!samesize(ima1,ima2,3,true)) throw BasisfieldException("splinefield::JtJ: Image dimension mismatch");
    if (mask && !samesize(ima1,*mask,3)) throw BasisfieldException("splinefield::JtJ: Mismatch between image and mask");
    if (FieldSz_x() != static_cast<unsigned int>(ima1.xsize()) ||
        FieldSz_y() != static_cast<unsigned int>(ima1.ysize()) ||
        FieldSz_z() != static_cast<unsigned int>(ima1.zsize())) throw BasisfieldException("splinefield::JtJ: Mismatch between image and field");
    float *prodima = new float[FieldSz()];
    hadamard(ima1,ima2,mask,prodima);
    std::vector<unsigned int>  isz(3,0);
    isz[0] = FieldSz_x(); isz[1] = FieldSz_y(); isz[2] = FieldSz_z();
    std::vector<unsigned int>  csz(3,0);
    csz[0] = CoefSz_x(); csz[1] = CoefSz_y(); csz[2] = CoefSz_z();
    Spline3D<double>           sp1(_sp.Order(),_sp.KnotSpacing(),deriv1);
    Spline3D<double>           sp2(_sp.Order(),_sp.KnotSpacing(),deriv2);
    if (_nthr==1) tmp = make_asymmetric_jtj(sp1,csz,sp2,csz,prodima,isz,prec);
    else tmp = make_asymmetric_jtj_para(sp1,csz,sp2,csz,prodima,isz,prec,_nthr);
    delete[] prodima;
  }
  return(tmp);
}

std::shared_ptr<BFMatrix> splinefield::JtJ(const NEWIMAGE::volume<float>&        ima1,
                                             const basisfield&                     bf2,      // Spline that determines column in JtJ
                                             const NEWIMAGE::volume<float>&        ima2,
                                             const NEWIMAGE::volume<char>          *mask,
                                             MISCMATHS::BFMatrixPrecisionType      prec)
const
{
  if (!samesize(ima1,ima2,3,true)) throw BasisfieldException("splinefield::JtJ: Image dimension mismatch");
  if (mask && !samesize(ima1,*mask,3)) throw BasisfieldException("splinefield::JtJ: Mismatch between image and mask");
  if (FieldSz_x() != static_cast<unsigned int>(ima1.xsize()) ||
      FieldSz_y() != static_cast<unsigned int>(ima1.ysize()) ||
      FieldSz_z() != static_cast<unsigned int>(ima1.zsize())) throw BasisfieldException("splinefield::JtJ: Mismatch between image and field");
  if (FieldSz_x() != bf2.FieldSz_x() || FieldSz_y() != bf2.FieldSz_y() || FieldSz_z() != FieldSz_z()) {
    throw BasisfieldException("splinefield::JtJ: Mismatch between fields");
  }

  std::shared_ptr<BFMatrix>   tmp;
  try {
    const splinefield&  tbf2 = dynamic_cast<const splinefield &>(bf2);

    float *prodima = new float[FieldSz()];
    hadamard(ima1,ima2,mask,prodima);

    std::vector<unsigned int>  isz(3,0);
    isz[0] = FieldSz_x(); isz[1] = FieldSz_y(); isz[2] = FieldSz_z();
    std::vector<unsigned int>  csz1(3,0);
    csz1[0] = CoefSz_x(); csz1[1] = CoefSz_y(); csz1[2] = CoefSz_z();
    std::vector<unsigned int>  csz2(3,0);
    csz2[0] = tbf2.CoefSz_x(); csz2[1] = tbf2.CoefSz_y(); csz2[2] = tbf2.CoefSz_z();

    if (_nthr==1) tmp = make_asymmetric_jtj(_sp,csz1,tbf2._sp,csz2,prodima,isz,prec);
    else tmp = make_asymmetric_jtj_para(_sp,csz1,tbf2._sp,csz2,prodima,isz,prec,_nthr);
    delete[] prodima;
  }
  catch (bad_cast) {
    throw BasisfieldException("splinefield::JtJ: Must pass like to like field");
  }

  return(tmp);

}


double splinefield::MemEnergy() const // Membrane energy of field
{
  const std::shared_ptr<NEWMAT::ColumnVector> lcoef = GetCoef();
  if (!lcoef) {throw BasisfieldException("splinefield::MemEnergy: No coefficients set yet");}

  NEWMAT::ColumnVector AtAb = 0.5 * MemEnergyGrad();
  double me = DotProduct(*lcoef,AtAb);
  return(me);
}

double splinefield::BendEnergy() const // Bending energy of field
{
  const std::shared_ptr<NEWMAT::ColumnVector> lcoef = GetCoef();
  if (!lcoef) {throw BasisfieldException("splinefield::MemEnergy: No coefficients set yet");}

  NEWMAT::ColumnVector AtAb = 0.5 * BendEnergyGrad();
  double be = DotProduct(*lcoef,AtAb);
  return(be);
}


NEWMAT::ReturnMatrix splinefield::MemEnergyGrad() const // Gradient of bending energy of field
{
  const std::shared_ptr<NEWMAT::ColumnVector> lcoef = GetCoef();
  if (!lcoef) {throw BasisfieldException("splinefield::MemEnergyGrad: No coefficients set yet");}

  std::vector<unsigned int>  csz(3,0);
  csz[0] = CoefSz_x(); csz[1] = CoefSz_y(); csz[2] = CoefSz_z();
  std::vector<unsigned int>  isz(3,0);
  isz[0] = FieldSz_x(); isz[1] = FieldSz_y(); isz[2] = FieldSz_z();
  std::vector<unsigned int>  ksp(3,0);
  ksp[0] = Ksp_x(); ksp[1] = Ksp_y(); ksp[2] = Ksp_z();
  std::vector<double> vxs(3,0);
  vxs[0] = Vxs_x(); vxs[1] = Vxs_y(); vxs[2] = Vxs_z();

  NEWMAT::ColumnVector  grad(CoefSz());

  calculate_memen_AtAb(*lcoef,ksp,isz,vxs,csz,Order(),grad);
  grad *= 2.0;

  grad.Release();
  return(grad);
}

NEWMAT::ReturnMatrix splinefield::BendEnergyGrad() const // Gradient of bending energy of field
{
  const std::shared_ptr<NEWMAT::ColumnVector> lcoef = GetCoef();
  if (!lcoef) {throw BasisfieldException("splinefield::BendEnergyGrad: No coefficients set yet");}

  std::vector<unsigned int>  csz(3,0);
  csz[0] = CoefSz_x(); csz[1] = CoefSz_y(); csz[2] = CoefSz_z();
  std::vector<unsigned int>  isz(3,0);
  isz[0] = FieldSz_x(); isz[1] = FieldSz_y(); isz[2] = FieldSz_z();
  std::vector<unsigned int>  ksp(3,0);
  ksp[0] = Ksp_x(); ksp[1] = Ksp_y(); ksp[2] = Ksp_z();
  std::vector<double> vxs(3,0);
  vxs[0] = Vxs_x(); vxs[1] = Vxs_y(); vxs[2] = Vxs_z();

  NEWMAT::ColumnVector  grad(CoefSz());

  calculate_bender_AtAb(*lcoef,ksp,isz,vxs,csz,Order(),grad);
  grad *= 2.0;

  grad.Release();
  return(grad);
}

std::shared_ptr<BFMatrix> splinefield::MemEnergyHess(MISCMATHS::BFMatrixPrecisionType   prec) const  // Hessian of membrane energy
{
  std::vector<unsigned int>    lksp(3,0);
  lksp[0] = Ksp_x(); lksp[1] = Ksp_y(); lksp[2] = Ksp_z();
  std::vector<unsigned int>    csz(3,0);
  csz[0] = CoefSz_x(); csz[1] = CoefSz_y(); csz[2] = CoefSz_z();
  unsigned int ncoef = csz[2]*csz[1]*csz[0];
  std::vector<unsigned int>    isz(3,0);
  isz[0] = FieldSz_x(); isz[1] = FieldSz_y(); isz[2] = FieldSz_z();

  unsigned int *irp; unsigned int *jcp; double *valp;
  if (_nthr==1) calculate_memen_bender_H(lksp,csz,isz,MemEn,irp,jcp,valp);
  else calculate_memen_bender_H_para(lksp,csz,isz,MemEn,_nthr,irp,jcp,valp);
  std::shared_ptr<BFMatrix>  H;
  if (prec==BFMatrixFloatPrecision) {
    H = std::shared_ptr<BFMatrix>(new SparseBFMatrix<float>(ncoef,ncoef,irp,jcp,valp,Utilities::NoOfThreads(_nthr)));
  }
  else H = std::shared_ptr<BFMatrix>(new SparseBFMatrix<double>(ncoef,ncoef,irp,jcp,valp,Utilities::NoOfThreads(_nthr)));
  delete [] irp; delete [] jcp; delete [] valp;

  return(H);
}

std::shared_ptr<BFMatrix> splinefield::BendEnergyHess(MISCMATHS::BFMatrixPrecisionType   prec) const // Hessian of bending energy
{
  std::vector<unsigned int>    lksp(3,0);
  lksp[0] = Ksp_x(); lksp[1] = Ksp_y(); lksp[2] = Ksp_z();
  std::vector<unsigned int>    csz(3,0);
  csz[0] = CoefSz_x(); csz[1] = CoefSz_y(); csz[2] = CoefSz_z();
  unsigned int ncoef = csz[2]*csz[1]*csz[0];
  std::vector<unsigned int>    isz(3,0);
  isz[0] = FieldSz_x(); isz[1] = FieldSz_y(); isz[2] = FieldSz_z();

  unsigned int *irp; unsigned int *jcp; double *valp;
  if (_nthr==1) calculate_memen_bender_H(lksp,csz,isz,BendEn,irp,jcp,valp);
  else calculate_memen_bender_H_para(lksp,csz,isz,BendEn,_nthr,irp,jcp,valp);

  std::shared_ptr<BFMatrix>  H;
  if (prec==BFMatrixFloatPrecision) {
    H = std::shared_ptr<BFMatrix>(new SparseBFMatrix<float>(ncoef,ncoef,irp,jcp,valp,Utilities::NoOfThreads(_nthr)));
  }
  else H = std::shared_ptr<BFMatrix>(new SparseBFMatrix<double>(ncoef,ncoef,irp,jcp,valp,Utilities::NoOfThreads(_nthr)));
  delete [] irp; delete [] jcp; delete [] valp;

  return(H);
}

std::shared_ptr<MISCMATHS::SpMat<float> > splinefield::BendEnergyHessAsSpMat() const // Hessian of bending energy
{
  std::vector<unsigned int>    lksp(3,0);
  lksp[0] = Ksp_x(); lksp[1] = Ksp_y(); lksp[2] = Ksp_z();
  std::vector<unsigned int>    csz(3,0);
  csz[0] = CoefSz_x(); csz[1] = CoefSz_y(); csz[2] = CoefSz_z();
  unsigned int ncoef = csz[2]*csz[1]*csz[0];
  std::vector<unsigned int>    isz(3,0);
  isz[0] = FieldSz_x(); isz[1] = FieldSz_y(); isz[2] = FieldSz_z();

  unsigned int *irp; unsigned int *jcp; double *valp;
  if (_nthr==1) calculate_memen_bender_H(lksp,csz,isz,BendEn,irp,jcp,valp); // N.B. irp, jcp and valp passed by reference
  else calculate_memen_bender_H_para(lksp,csz,isz,BendEn,_nthr,irp,jcp,valp);
  std::shared_ptr<MISCMATHS::SpMat<float> > H = std::shared_ptr<MISCMATHS::SpMat<float> >(new SpMat<float>(ncoef,ncoef,irp,jcp,valp,Utilities::NoOfThreads(_nthr)));
  delete [] irp; delete [] jcp; delete [] valp;

  return(H);
}

/////////////////////////////////////////////////////////////////////
//
// Calculates field from cofficients. It loops over coefficients
// instead of over voxels in the resulting field. This makes it
// more difficult to parallelise. For that reason I have instead
// added a function get_fields (below) that allows one to calculate
// more than one field at a time, allowing for "re-use" of much of
// the calculations.
//
/////////////////////////////////////////////////////////////////////

void splinefield::get_field(const Spline3D<double>&           sp,
                            const double                      *coef,
                            const std::vector<unsigned int>&  csz,
                            const std::vector<unsigned int>&  fsz,
                            double                            *fld) const
{
  std::vector<unsigned int>    ff(3,0);    // First index into field
  std::vector<unsigned int>    lf(3,0);    // Last index into field
  std::vector<unsigned int>    os(3,0);    // Offset into spline
  std::vector<unsigned int>    ci(3,0);    // Coefficient/spline index
  std::vector<unsigned int>    ks(3,0);    // Kernel size

  for (int i=0; i<3; i++) ks[i]=sp.KernelSize(i);
  memset(fld,0,fsz[0]*fsz[1]*fsz[2]*sizeof(double));

  for (ci[2]=0; ci[2]<csz[2]; ci[2]++) {
    for (ci[1]=0; ci[1]<csz[1]; ci[1]++) {
      for (ci[0]=0; ci[0]<csz[0]; ci[0]++) {
        sp.RangeInField(ci,fsz,ff,lf);
        sp.OffsetIntoKernel(ci,fsz,os);
        double c = coef[ci[2]*csz[1]*csz[0] + ci[1]*csz[0] + ci[0]];
        for (unsigned int fk=ff[2], sk=os[2]; fk<lf[2]; fk++, sk++) {
          for (unsigned int fj=ff[1], sj=os[1]; fj<lf[1]; fj++, sj++) {
            unsigned int fbi = fk*fsz[1]*fsz[0] + fj*fsz[0];
            unsigned int sbi = sk*ks[1]*ks[0] + sj*ks[0];
            for (unsigned int fi=ff[0], si=os[0]; fi<lf[0]; fi++, si++) {
              fld[fbi+fi] += c * sp[sbi+si];
	    }
	  }
	}
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////
//
// Calculates fields from cofficients. It performs the same calculations
// as the function above. But it takes a vector of splines and a vector
// of fields, thereby allowing for reuse of much of the calculations.
//
/////////////////////////////////////////////////////////////////////

void splinefield::get_fields(const std::vector<Spline3D<double> >& spv,
                            const double                           *coef,
                            const std::vector<unsigned int>&       csz,
                            const std::vector<unsigned int>&       fsz,
			    std::vector<double *>                  fldv) const
{
  std::vector<unsigned int>    ff(3,0);    // First index into field
  std::vector<unsigned int>    lf(3,0);    // Last index into field
  std::vector<unsigned int>    os(3,0);    // Offset into spline
  std::vector<unsigned int>    ci(3,0);    // Coefficient/spline index
  std::vector<unsigned int>    ks(3,0);    // Kernel size

  for (int i=0; i<3; i++) ks[i]=spv[0].KernelSize(i);
  for (unsigned int i=0; i<fldv.size(); i++) memset(fldv[i],0,fsz[0]*fsz[1]*fsz[2]*sizeof(double));

  for (ci[2]=0; ci[2]<csz[2]; ci[2]++) {
    for (ci[1]=0; ci[1]<csz[1]; ci[1]++) {
      for (ci[0]=0; ci[0]<csz[0]; ci[0]++) {
        spv[0].RangeInField(ci,fsz,ff,lf);
        spv[0].OffsetIntoKernel(ci,fsz,os);
        double c = coef[ci[2]*csz[1]*csz[0] + ci[1]*csz[0] + ci[0]];
        for (unsigned int fk=ff[2], sk=os[2]; fk<lf[2]; fk++, sk++) {
          for (unsigned int fj=ff[1], sj=os[1]; fj<lf[1]; fj++, sj++) {
            unsigned int fbi = fk*fsz[1]*fsz[0] + fj*fsz[0];
            unsigned int sbi = sk*ks[1]*ks[0] + sj*ks[0];
            for (unsigned int fi=ff[0], si=os[0]; fi<lf[0]; fi++, si++) {
	      for (unsigned int vi=0; vi<spv.size(); vi++) {
		fldv[vi][fbi+fi] += c * spv[vi][sbi+si];
	      }
	    }
	  }
	}
      }
    }
  }
}

template<class T>
void splinefield::get_jte(const Spline3D<double>&            sp,
                          const std::vector<unsigned int>&   csz,
                          const T                            *ima,
                          const std::vector<unsigned int>&   isz,
                          double                             *jte) const
{
  std::vector<unsigned int>    fi(3,0);    // First index into images
  std::vector<unsigned int>    li(3,0);    // Last index into images
  std::vector<unsigned int>    os(3,0);    // Offset into spline
  std::vector<unsigned int>    ci(3,0);    // Coefficient/spline index
  std::vector<unsigned int>    ks(3,0);    // Kernel size

  for (int i=0; i<3; i++) ks[i]=sp.KernelSize(i);

  memset(jte,0,csz[0]*csz[1]*csz[2]*sizeof(double));

  for (ci[2]=0; ci[2]<csz[2]; ci[2]++) {
    for (ci[1]=0; ci[1]<csz[1]; ci[1]++) {
      for (ci[0]=0; ci[0]<csz[0]; ci[0]++) {
        sp.RangeInField(ci,isz,fi,li);
        sp.OffsetIntoKernel(ci,isz,os);
        double *jtep = &jte[ci[2]*csz[1]*csz[0] + ci[1]*csz[0] + ci[0]];
        for (unsigned int ik=fi[2], sk=os[2]; ik<li[2]; ik++, sk++) {
          for (unsigned int ij=fi[1], sj=os[1]; ij<li[1]; ij++, sj++) {
            unsigned int ibi = ik*isz[1]*isz[0] + ij*isz[0];
            unsigned int sbi = sk*ks[1]*ks[0] + sj*ks[0];
            for (unsigned int ii=fi[0], si=os[0]; ii<li[0]; ii++, si++) {
              *jtep += sp[sbi+si] * static_cast<double>(ima[ibi+ii]);
	    }
	  }
	}
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////
//
// This routine, and the next one, performs the exact same calculations
// as get_jte above, but is parallelised using the C++11 thread library.
// The main routine, make_jte_para, is responsible for allocating all the
// objects  that are to be used, and then to divide up the work between
// the threads that it will spawn.
//
/////////////////////////////////////////////////////////////////////

template<class T>
void splinefield::get_jte_para(const Spline3D<double>&            sp,   // Spline
			       const std::vector<unsigned int>&   csz,  // Size of Coefficient matrix for csp
			       const T                            *ima, // Image
			       const std::vector<unsigned int>&   isz,  // Image size
			       double                             *jte, // Jte
			       unsigned int                       nthr) // Number of threads
const
{
  // Vector of threads
  std::vector<std::thread> threads(nthr-1); // + main thread makes nthr

  // Decide no of columns per thread
  double cpt = static_cast<double>(csz[0]*csz[1]*csz[2]) / static_cast<double>(nthr);
  std::vector<unsigned int> ncols(nthr+1,0); //
  for (unsigned int i=1; i<nthr; i++) ncols[i] = static_cast<unsigned int>(std::floor(i*cpt));
  ncols[nthr] = csz[0]*csz[1]*csz[2];

  // Next spawn threads to do calculations
  for (unsigned int i=0; i<nthr-1; i++) {
    threads[i] = std::thread(&splinefield::calculate_jte_elements<T>,this,ncols[i],ncols[i+1],
			     std::ref(sp),std::ref(csz),ima,std::ref(isz),jte);
  }
  this->calculate_jte_elements(ncols[nthr-1],ncols[nthr],sp,csz,ima,isz,jte);
  std::for_each(threads.begin(),threads.end(),std::mem_fn(&std::thread::join));
}

template <class T>
void splinefield::calculate_jte_elements(// Input
					 unsigned int                      first_ci, // First column index
					 unsigned int                      last_ci,  // One past last column index
					 const Spline3D<double>&           sp,       // Spline
					 const std::vector<unsigned int>&  csz,      // Size of Coefficient matrix for csp
					 const T                           *ima,     // Image
					 const std::vector<unsigned int>&  isz,      // Image size
					 // Output
					 double                            *jte)     // Well, Jte.
const
{
  std::vector<unsigned int>    fi(3,0);    // First index into images
  std::vector<unsigned int>    li(3,0);    // Last index into images
  std::vector<unsigned int>    os(3,0);    // Offset into spline
  std::vector<unsigned int>    ci(3,0);    // Coefficient/spline index
  std::vector<unsigned int>    ks(3,0);    // Kernel size

  for (int i=0; i<3; i++) ks[i]=sp.KernelSize(i);

  for (unsigned int i=first_ci; i<last_ci; i++) {
    ci[2] = i / (csz[0]*csz[1]);                   // Truncated integer division
    ci[1] = (i-ci[2]*csz[0]*csz[1]) / csz[0];   // Truncated integer division
    ci[0] = i - ci[2]*csz[0]*csz[1] - ci[1]*csz[0];
    sp.RangeInField(ci,isz,fi,li);
    sp.OffsetIntoKernel(ci,isz,os);
    double *jtep = &jte[ci[2]*csz[1]*csz[0] + ci[1]*csz[0] + ci[0]];
    *jtep = 0.0;
    for (unsigned int ik=fi[2], sk=os[2]; ik<li[2]; ik++, sk++) {
      for (unsigned int ij=fi[1], sj=os[1]; ij<li[1]; ij++, sj++) {
	unsigned int ibi = ik*isz[1]*isz[0] + ij*isz[0];
	unsigned int sbi = sk*ks[1]*ks[0] + sj*ks[0];
	for (unsigned int ii=fi[0], si=os[0]; ii<li[0]; ii++, si++) {
	  *jtep += sp[sbi+si] * static_cast<double>(ima[ibi+ii]);
	}
      }
    }
  }
}

std::shared_ptr<SpMat<float> > splinefield::get_j(const Spline3D<double>&            sp,
                                                  const std::vector<unsigned int>&   csz,
                                                  const float                        *ima,
                                                  const std::vector<unsigned int>&   isz) const
{
  unsigned int                 m = isz[2]*isz[1]*isz[0];     // # of rows in J
  unsigned int                 n = csz[2]*csz[1]*csz[0];     // # of columns in J
  unsigned int                 nnz = sp.NzMaxInA(isz);       // Max # of non-zero elements
  unsigned int                 *irp = new unsigned int[nnz]; // Row indices
  unsigned int                 *jcp = new unsigned int[n+1]; // Indicies into irp indicating start/stop of columns
  double                       *valp = new double[nnz];      // The values of the matrix
  std::vector<unsigned int>    ks(3,0);                      // Kernel size
  std::vector<unsigned int>    ci(3,0);                      // Coefficient/spline index
  std::vector<unsigned int>    fi(3,0);                      // First index into images
  std::vector<unsigned int>    li(3,0);                      // Last index into images
  std::vector<unsigned int>    os(3,0);                      // Offset into spline

  for (int i=0; i<3; i++) ks[i]=sp.KernelSize(i);

  unsigned int ivindx = 0;  // Index into irp and valp
  unsigned int jindx = 0;   // Index into jcp
  for (ci[2]=0; ci[2]<csz[2]; ci[2]++) {
    for (ci[1]=0; ci[1]<csz[1]; ci[1]++) {
      for (ci[0]=0; ci[0]<csz[0]; ci[0]++) {
	jcp[jindx++] = ivindx;
	sp.RangeInField(ci,isz,fi,li);
	sp.OffsetIntoKernel(ci,isz,os);
        // ik denotes index into ima for third dimension, sk denotes index into spline
	// kernel for third dimension. ij ans sj similarly for second dimension etc.
	for (unsigned int ik=fi[2], sk=os[2]; ik<li[2]; ik++, sk++) {
	  for (unsigned int ij=fi[1], sj=os[1]; ij<li[1]; ij++, sj++) {
	    unsigned int ibi = ik*isz[1]*isz[0] + ij*isz[0];  // ImageBaseIndex
	    unsigned int sbi = sk*ks[1]*ks[0] + sj*ks[0];     // SplineBaseIndex
	    for (unsigned int ii=fi[0], si=os[0]; ii<li[0]; ii++, si++) {
	      irp[ivindx] = ibi+ii;
	      valp[ivindx++] = sp[sbi+si] * static_cast<double>(ima[ibi+ii]);
	    }
	  }
	}
      }
    }
  }
  jcp[jindx] = ivindx;

  std::shared_ptr<SpMat<float> > j = std::shared_ptr<SpMat<float> >(new SpMat<float>(m,n,irp,jcp,valp,Utilities::NoOfThreads(_nthr)));

  delete [] irp; delete [] jcp; delete [] valp;

  return(j);
}

/////////////////////////////////////////////////////////////////////
//
// This routines calculates A'*B where A is an nxm matrix where n is the
// number of voxels in ima and m is the number of splines of one kind
// and B is an nxl matrix where l is the number of splines of a different
// kind.
//
// The splines may e.g. be a spline of some ksp in A and the same
// ksp spline differentiated in one direction in B. This can then be
// used for modelling effects of Jacobian modulation in distortion
// correction. In this first case jtj is still square, though not
// symmetrical.
//
// The other possibility is that A has splines of ksp1 modelling e.g.
// displacements and B has splines of ksp2 modelling e.g. an intensity
// bias field. In this case jtj is no longer square.
//
// This routine does not utilise symmetries/repeated values at any
// level. For the second case above there is nothing to utilise (as
// far as I can tell). For the former case there are complicated
// patterns of repeated values that it should be possible to utilise
// in order to speed things up. Future improvments.
//
/////////////////////////////////////////////////////////////////////

template<class T>
std::shared_ptr<BFMatrix> splinefield::make_asymmetric_jtj(const Spline3D<double>&           rsp,   // Spline that determines row in JtJ
                                                             const std::vector<unsigned int>&  cszr,  // Coefficient matrix size for rsp
                                                             const Spline3D<double>&           sp2,   // Spline that determines column in JtJ
                                                             const std::vector<unsigned int>&  cszc,  // Coefficient matrix size for csp/sp2
                                                             const T                           *ima,  // Image (typically product of two images)
                                                             const std::vector<unsigned int>&  isz,   // Matrix size of image
                                                             MISCMATHS::BFMatrixPrecisionType  prec)  // Precision (float/double) of resulting matrix
const
{
  Spline3D<double>             csp(sp2);    // Read-write copy of spline that determines column
  std::vector<unsigned int>    cindx(3,0);  // Index of spline that determines column in JtJ
  std::vector<unsigned int>    rindx(3,0);  // Index of spline that determines row in JtJ
  std::vector<unsigned int>    fo(3,0);     // First index of overlapping spline in x-, y- and z-direction
  std::vector<unsigned int>    lo(3,0);     // Last index of overlapping spline in x-, y- and z-direction

  unsigned int                 m = cszr[2]*cszr[1]*cszr[0];  // # of rows in JtJ
  unsigned int                 n = cszc[2]*cszc[1]*cszc[0];  // # of columns in JtJ
  unsigned int                 nnz = rsp.NzMax(isz,csp);     // Max # of non-zero elements
  unsigned int                 *irp = new unsigned int[nnz]; // Row indices
  unsigned int                 *jcp = new unsigned int[n+1]; // Indicies into irp indicating start/stop of columns
  double                       *valp = new double[nnz];      // The values of the matrix

  unsigned int vali = 0;     // Index of present non-zero value (linear indexing)
  unsigned int ci = 0;       // Column index

  ZeroSplineMap   r_zeromap(rsp,cszr,ima,isz); // Keeps track of all splines for which ima is zero over whole support
  ZeroSplineMap   c_zeromap(csp,cszc,ima,isz); // Keeps track of all splines for which ima is zero over whole support

  // Same Kernel Size?
  bool sks = (rsp.KernelSize(0) == csp.KernelSize(0) && rsp.KernelSize(1) == csp.KernelSize(1) && rsp.KernelSize(2) == csp.KernelSize(2));

  for (cindx[2]=0; cindx[2]<cszc[2]; cindx[2]++) {
    for (cindx[1]=0; cindx[1]<cszc[1]; cindx[1]++) {
      for (cindx[0]=0; cindx[0]<cszc[0]; cindx[0]++) {
        ci = cindx[2]*cszc[1]*cszc[0] + cindx[1]*cszc[0] + cindx[0];
        jcp[ci] = vali;
        bool c_is_zero = c_zeromap(cindx);
        if (!c_is_zero) csp.Premul(cindx,isz,ima);
        if (sks) csp.RangeOfOverlappingSplines(cindx,isz,fo,lo);
        else csp.RangeOfOverlappingSplines(cindx,isz,rsp,fo,lo);
        for (unsigned int k=fo[2]; k<lo[2]; k++) {
          for (unsigned int j=fo[1]; j<lo[1]; j++) {
            for (unsigned int i=fo[0]; i<lo[0]; i++) {
              unsigned int ri = k*cszr[1]*cszr[0] + j*cszr[0] + i;
              rindx[0]=i; rindx[1]=j; rindx[2]=k;
              irp[vali] = ri;
              if (c_is_zero || r_zeromap(rindx)) valp[vali++] = 0;
              else valp[vali++] = csp.MulByOther(cindx,rindx,isz,rsp);
	    }
	  }
	}
      }
    }
  }
  jcp[ci+1] = vali;

  std::shared_ptr<BFMatrix>   jtj;
  if (prec==BFMatrixFloatPrecision) jtj = std::shared_ptr<BFMatrix>(new SparseBFMatrix<float>(m,n,irp,jcp,valp));
  else jtj = std::shared_ptr<BFMatrix>(new SparseBFMatrix<double>(m,n,irp,jcp,valp));

  delete [] irp; delete [] jcp; delete [] valp;

  return(jtj);
}

/////////////////////////////////////////////////////////////////////
//
// This routine, and the next two, performs the exact same calculations
// as make_asymmetric_jtj above, but is parallelised using the C++11
// thread library. The main routine, make_asymmetric_jtj_para, is
// responsible for allocating all the structures that are to be used,
// and then to divide up the work between the threads that it will
// spawn.
// It will use two functions to achieve this, both of which are spawned
// by the main routine. The first is used to calculate the total number
// of non-zero elements that each thread will be responsible for
// calculating. These numbers will be passed into the second function
// that will use that as an offset for where to start inserting the
// non-zero values that it is also responsible for calculating.
//
/////////////////////////////////////////////////////////////////////

template<class T>
std::shared_ptr<BFMatrix> splinefield::make_asymmetric_jtj_para(const Spline3D<double>&           rsp,   // Spline that determines row in JtJ
                                                                const std::vector<unsigned int>&  cszr,  // Coefficient matrix size for rsp
                                                                const Spline3D<double>&           csp,   // Spline that determines column in JtJ
                                                                const std::vector<unsigned int>&  cszc,  // Coefficient matrix size for csp/sp2
                                                                const T                           *ima,  // Image (typically product of two images)
                                                                const std::vector<unsigned int>&  isz,   // Matrix size of image
                                                                MISCMATHS::BFMatrixPrecisionType  prec,  // Precision (float/double) of resulting matrix
                                                                unsigned int                      nthr)  // Number of threads
const
{
  unsigned int                 m = cszr[2]*cszr[1]*cszr[0];  // # of rows in JtJ
  unsigned int                 n = cszc[2]*cszc[1]*cszc[0];  // # of columns in JtJ
  unsigned int                 nnz = rsp.NzMax(isz,csp);     // Max # of non-zero elements
  unsigned int                 *irp = new unsigned int[nnz]; // Row indices
  unsigned int                 *jcp = new unsigned int[n+1]; // Indicies into irp indicating start/stop of columns
  double                       *valp = new double[nnz];      // The values of the matrix

  ZeroSplineMap   r_zeromap(rsp,cszr,ima,isz); // Keeps track of all splines for which ima is zero over whole support
  ZeroSplineMap   c_zeromap(csp,cszc,ima,isz); // Keeps track of all splines for which ima is zero over whole support

  // Same Kernel Size?
  bool sks = (rsp.KernelSize(0) == csp.KernelSize(0) && rsp.KernelSize(1) == csp.KernelSize(1) && rsp.KernelSize(2) == csp.KernelSize(2));

  // Vector of threads
  std::vector<std::thread> threads(nthr-1); // + main thread makes nthr

  // Start by deciding no of columns per thread
  double cpt = static_cast<double>(n) / static_cast<double>(nthr);
  std::vector<unsigned int> ncols(nthr+1,0); //
  for (unsigned int i=1; i<nthr; i++) ncols[i] = static_cast<unsigned int>(std::floor(i*cpt));
  ncols[nthr] = n;

  // Next determine how many non-zero elements each thread will calculate
  // We could possibly improve performance by having equal amounts of non-zero
  // elements per thread, instead of as now equal number of columns. Future work.
  std::vector<unsigned int> nnz_pt(nthr,0); // No. of non-zero elements per thread
  for (unsigned int i=0; i<nthr-1; i++) {
    threads[i] = std::thread(&splinefield::no_of_non_zero_elements_asym,this,ncols[i],ncols[i+1],std::ref(cszc),
			     std::ref(isz),std::ref(csp),std::ref(rsp),std::ref(nnz_pt[i]));
  }
  this->no_of_non_zero_elements_asym(ncols[nthr-1],ncols[nthr],cszc,isz,csp,rsp,nnz_pt[nthr-1]);
  std::for_each(threads.begin(),threads.end(),std::mem_fn(&std::thread::join));

  // Convert to start-index per thread
  std::vector<unsigned int> vali(nthr,0);
  for (unsigned int i=1; i<nthr; i++) vali[i] = vali[i-1] + nnz_pt[i-1];

  // Then actually calculate those non-zero elements
  for (unsigned int i=0; i<nthr-1; i++) {
    threads[i] = std::thread(&splinefield::calculate_non_zero_elements_asym<T>,this,ncols[i],ncols[i+1],vali[i],
			     std::ref(cszc),std::ref(cszr),std::ref(isz),csp,std::ref(rsp),sks,
			     std::ref(r_zeromap),std::ref(c_zeromap),ima,irp,jcp,valp);
  }
  this->calculate_non_zero_elements_asym(ncols[nthr-1],ncols[nthr],vali[nthr-1],cszc,cszr,isz,
					 csp,rsp,sks,r_zeromap,c_zeromap,ima,irp,jcp,valp);
  std::for_each(threads.begin(),threads.end(),std::mem_fn(&std::thread::join));

  std::shared_ptr<BFMatrix>   jtj;
  if (prec==BFMatrixFloatPrecision) {
    jtj = std::shared_ptr<BFMatrix>(new SparseBFMatrix<float>(m,n,irp,jcp,valp,Utilities::NoOfThreads(nthr)));
  }
  else jtj = std::shared_ptr<BFMatrix>(new SparseBFMatrix<double>(m,n,irp,jcp,valp,Utilities::NoOfThreads(nthr)));

  delete [] irp; delete [] jcp; delete [] valp;

  return(jtj);
}

void splinefield::no_of_non_zero_elements_asym(// Input
					       unsigned int  first_ci,                // First column index
					       unsigned int  last_ci,                 // One past last column index
					       const std::vector<unsigned int>&  csz, // Size of Coefficient matrix for csp
					       const std::vector<unsigned int>&  isz, // Image size
					       const Spline3D<double>&           csp, // Spline that determines column in JtJ
					       const Spline3D<double>&           rsp, // Spline that determines row in JtJ
					       // Output
					       unsigned int& nnz)                     // Total no. of non-zero elements in columns first_ci--last_ci-1
const
{
  nnz = 0;
  std::vector<unsigned int> cindx(3,0);
  for (unsigned int i=first_ci; i<last_ci; i++) {
    cindx[2] = i / (csz[0]*csz[1]);                   // Truncated integer division
    cindx[1] = (i-cindx[2]*csz[0]*csz[1]) / csz[0];   // Truncated integer division
    cindx[0] = i - cindx[2]*csz[0]*csz[1] - cindx[1]*csz[0];
    nnz += csp.NumberOfOverlappingSplines(cindx,isz,rsp);
  }
}

template<class T>
void splinefield::calculate_non_zero_elements_asym(// Input
						   unsigned int                      first_ci, // First column index
						   unsigned int                      last_ci,  // One past last column index
						   unsigned int                      vali,     // First index into irp and valp
						   const std::vector<unsigned int>&  cszc,     // Size of coefficient matrix for "column splines"
						   const std::vector<unsigned int>&  cszr,     // Size of coefficient matrix for "row splines"
						   const std::vector<unsigned int>&  isz,      // Size of image matrix
						   Spline3D<double>                  csp,      // Spline defining columns
						   const Spline3D<double>&           rsp,      // Spline defining rows
						   bool                              sks,      // True if csp and rsp has same size
						   const ZeroSplineMap&              rzm,      // Keeps track of splines for which ima is zero over entire support
						   const ZeroSplineMap&              czm,      // Keeps track of splines for which ima is zero over entire support
						   const T                           *ima,     // Image
						   // Output
						   unsigned int                      *irp,     // nzmax length array of row indicies
						   unsigned int                      *jcp,     // n+1 length array of starts of column indicies
						   double                            *valp)    // nzmax length array of values
const
{
  std::vector<unsigned int>    cindx(3,0);  // Index into coefficient matrix of spline that determines column in JtJ
  std::vector<unsigned int>    rindx(3,0);  // Index into coefficient matrix of spline that determines row in JtJ
  std::vector<unsigned int>    fo(3,0);     // First index of overlapping spline in x-, y- and z-direction
  std::vector<unsigned int>    lo(3,0);     // Last index of overlapping spline in x-, y- and z-direction
  unsigned int                 ci = 0;      // Column index

  for (ci=first_ci; ci<last_ci; ci++) {
    cindx[2] = ci / (cszc[0]*cszc[1]);                   // Truncated integer division
    cindx[1] = (ci-cindx[2]*cszc[0]*cszc[1]) / cszc[0];   // Truncated integer division
    cindx[0] = ci - cindx[2]*cszc[0]*cszc[1] - cindx[1]*cszc[0];
    jcp[ci] = vali;
    bool c_is_zero = czm(cindx);
    if (!c_is_zero) csp.Premul(cindx,isz,ima);
    if (sks) csp.RangeOfOverlappingSplines(cindx,isz,fo,lo);
    else csp.RangeOfOverlappingSplines(cindx,isz,rsp,fo,lo);
    for (unsigned int k=fo[2]; k<lo[2]; k++) {
      for (unsigned int j=fo[1]; j<lo[1]; j++) {
	for (unsigned int i=fo[0]; i<lo[0]; i++) {
	  unsigned int ri = k*cszr[1]*cszr[0] + j*cszr[0] + i;
	  rindx[0]=i; rindx[1]=j; rindx[2]=k;
	  irp[vali] = ri;
	  if (c_is_zero || rzm(rindx)) valp[vali++] = 0;
	  else valp[vali++] = csp.MulByOther(cindx,rindx,isz,rsp);
	}
      }
    }
  }
  if (ci == cszc[2]*cszc[1]*cszc[0]) jcp[ci] = vali; // If this is the "last" thread
}

/////////////////////////////////////////////////////////////////////
//
// This routines calculates A'*A where A is an nxm matrix where n is the
// number of voxels in ima and m is the number of splines. Each column
// of A is a spline elementwise multiplied by by the image. A is
// typically too large to be represented, each column being the size
// of the image volume and there typicall being tens of thousands of
// columns. Therefore only A'*A is explicitly represented, and even that
// as a sparse matrix.
// The routine looks a little complicated. If not using symmetry it is
// actually quite straightforward. However, only ~1/8 of the elements
// are unique due to there being three levels of symmetry. In order to
// maximise efficiency I have utilised this symmetry, which sadly leads
// to lots of book keeping.
//
/////////////////////////////////////////////////////////////////////

template<class T>
std::shared_ptr<BFMatrix> splinefield::make_fully_symmetric_jtj(const Spline3D<double>&            sp2,
                                                                  const std::vector<unsigned int>&   csz,
                                                                  const T                            *ima,
                                                                  const std::vector<unsigned int>&   isz,
                                                                  MISCMATHS::BFMatrixPrecisionType   prec)
const
{
  Spline3D<double>           sp1(sp2);      // Another copy of spline
  std::vector<unsigned int>  indx1(3,0);    // Index of spline that determines column in JtJ
  std::vector<unsigned int>  indx2(3,0);    // Index of spline that determines row in JtJ
  std::vector<unsigned int>  fo(3,0);       // First index of overlapping spline in x-, y- and z-direction
  std::vector<unsigned int>  lo(3,0);       // Last index of overlapping spline in x-, y- and z-direction

  unsigned int               ncoef = csz[2]*csz[1]*csz[0];     // Size of JtJ
  unsigned int               nnz = sp1.NzMax(isz);             // Max # of non-zero elements
  unsigned int               *irp = new unsigned int[nnz];     // Row indicies
  unsigned int               *jcp = new unsigned int[ncoef+1]; // Indicies into irp indicating start/stop of columns
  double                     *valp = new double[nnz];          // The values of the matrix

  unsigned int vali = 0;     // Index of present non-zero value (linear indexing)
  unsigned int ci = 0;       // Column index

  ZeroSplineMap  zeromap(sp1,csz,ima,isz);

  for (indx1[2]=0; indx1[2]<csz[2]; indx1[2]++) {
    for (indx1[1]=0; indx1[1]<csz[1]; indx1[1]++) {
      for (indx1[0]=0; indx1[0]<csz[0]; indx1[0]++) {
        ci = indx1[2]*csz[1]*csz[0] + indx1[1]*csz[0] + indx1[0];
        jcp[ci] = vali;
        bool indx1_is_zero = zeromap(indx1);
        if (!indx1_is_zero) sp1.Premul(indx1,isz,ima);
        sp1.RangeOfOverlappingSplines(indx1,isz,fo,lo);
        //
        // Fill in values above the main diagonal,
        // utilising the top level of symmetry.
        //
        for (unsigned int k=fo[2]; k<indx1[2]; k++) {
          for (unsigned int j=fo[1]; j<lo[1]; j++) {
            for (unsigned int i=fo[0]; i<lo[0]; i++) {
              unsigned int ri = k*csz[1]*csz[0] + j*csz[0] + i;
              irp[vali] = ri;
              if (indx1_is_zero) valp[vali++] = 0;
              else valp[vali++] = get_val(ci,ri,irp,jcp,valp);
	    }
	  }
	}
        for (unsigned int k=indx1[2]; k<lo[2]; k++) {
          //
          // Fill in values above the "main" diagonals at the
          // 2nd level of symmetry. N.B. that the values should
          // NOT be mirrored around the main diagonal.
          //
          for (unsigned int j=fo[1]; j<indx1[1]; j++) {
            for (unsigned int i=fo[0]; i<lo[0]; i++) {
              unsigned int ri = k*csz[1]*csz[0] + j*csz[0] + i;
              unsigned int cri = indx1[2]*csz[1]*csz[0] + j*csz[0] + i;
              unsigned int cci = k*csz[1]*csz[0] + indx1[1]*csz[0] + indx1[0];
              irp[vali] = ri;
              if (indx1_is_zero) valp[vali++] = 0;
              else valp[vali++] = get_val(cci,cri,irp,jcp,valp);
	    }
	  }

	  for (unsigned int j=indx1[1]; j<lo[1]; j++) {
            //
            // Fill in values above the "main" diagonals at the third
            // and final level of symmetry. Same N.B. as above applies.
            //
	    for (unsigned int i=fo[0]; i<indx1[0]; i++) {
              unsigned int ri = k*csz[1]*csz[0] + j*csz[0] + i;
              unsigned int cri = indx1[2]*csz[1]*csz[0]+indx1[1]*csz[0]+i;
              unsigned int cci = k*csz[1]*csz[0]+j*csz[0]+indx1[0];
              irp[vali] = ri;
              if (indx1_is_zero) valp[vali++] = 0;
              else valp[vali++] = get_val(cci,cri,irp,jcp,valp);
	    }
            //
            // And these are the positions for which we actually need to
            // calculate new values. Roughly ~1/8 of the total.
            //
	    for (unsigned int i=indx1[0]; i<lo[0]; i++) {
              unsigned int ri = k*csz[1]*csz[0] + j*csz[0] + i;
              indx2[0]=i; indx2[1]=j; indx2[2]=k;
              irp[vali] = ri;
              if (indx1_is_zero || zeromap(indx2)) valp[vali++] = 0;
              else valp[vali++] = sp1.MulByOther(indx1,indx2,isz,sp2);
	    }
	  }
	}
      }
    }
  }
  jcp[ci+1] = vali;

  std::shared_ptr<BFMatrix>   jtj;
  if (prec==BFMatrixFloatPrecision) jtj = std::shared_ptr<BFMatrix>(new SparseBFMatrix<float>(ncoef,ncoef,irp,jcp,valp));
  else jtj = std::shared_ptr<BFMatrix>(new SparseBFMatrix<double>(ncoef,ncoef,irp,jcp,valp));

  delete [] irp; delete [] jcp; delete [] valp;

  return(jtj);
}

/////////////////////////////////////////////////////////////////////
//
// Helper routine to find the value for a given row-index in a
// (possibly unfinished) compressed column storage format. Uses
// bisection.
//
/////////////////////////////////////////////////////////////////////

double splinefield::get_val(unsigned int           row,    // The row we want to find the value of
                            unsigned int           col,    // The column we want to find the value of
                            const unsigned int     *irp,   // Array of row-indicies
                            const unsigned int     *jcp,   // Array of indicies into irp
                            const double           *val)   // Array of values sorted as irp
const
{
   const unsigned int  *a = &(irp[jcp[col]]);
   const double        *v = &(val[jcp[col]]);
   int                 n = jcp[col+1]-jcp[col];
   int                 j = 0;
   int                 jlo = -1;
   int                 jup = n;

   if (row < a[0] || row > a[n-1]) {return(0.0);}

   while (jup-jlo > 1)
   {
      j = (jlo+jup) >> 1;
      if (row >= a[j]) {jlo = j;}
      else {jup = j;}
   }

   if (a[jlo] == row) {return(v[jlo]);}
   else return(0.0);
}

/////////////////////////////////////////////////////////////////////
//
// This routine, and the next three, performs the exact same calculations
// as make_fully_symmetric_jtj above, but is parallelised using the C++11
// thread library. The main routine, make_fully_symmetric_jtj_para, is
// responsible for allocating all the structures that are to be used,
// and then to divide up the work between the threads that it will
// spawn.
// It will use three functions to achieve this, all of which are spawned
// by the main routine.
// The first is used to calculate the total number of non-zero elements
// that each thread will be responsible for calculating.
// These numbers will be passed into the second function that will use
// that as an offset for where to start inserting the non-zero values that
// it is responsible for calculating. The second function will _only_
// calculate the unique non-zero elements.
// The third routine finally is responsible for finding and filling in
// the numbers that have already been calculated. It takes advantage of
// several levels of symmetry to do that.
//
/////////////////////////////////////////////////////////////////////

template<class T>
std::shared_ptr<BFMatrix> splinefield::make_fully_symmetric_jtj_para(const Spline3D<double>&            sp,
                                                                     const std::vector<unsigned int>&   csz,
                                                                     const T                            *ima,
                                                                     const std::vector<unsigned int>&   isz,
                                                                     MISCMATHS::BFMatrixPrecisionType   prec,
                                                                     unsigned int                       nthr)  // Number of threads)
const
{
  unsigned int               ncoef = csz[2]*csz[1]*csz[0];     // Size of JtJ
  unsigned int               nnz = sp.NzMax(isz);             // Max # of non-zero elements
  unsigned int               *irp = new unsigned int[nnz];     // Row indicies
  unsigned int               *jcp = new unsigned int[ncoef+1]; // Indicies into irp indicating start/stop of columns
  double                     *valp = new double[nnz];          // The values of the matrix

  ZeroSplineMap  zeromap(sp,csz,ima,isz);

  // Vector of threads
  std::vector<std::thread> threads(nthr-1); // + main thread makes nthr

  // Start by deciding no of columns per thread
  double cpt = static_cast<double>(ncoef) / static_cast<double>(nthr);
  std::vector<unsigned int> ncols(nthr+1,0); //
  for (unsigned int i=1; i<nthr; i++) ncols[i] = static_cast<unsigned int>(std::floor(i*cpt));
  ncols[nthr] = ncoef;

  // Next determine how many non-zero elements each thread will be responsible
  // for (calculate or fill in).
  // We could possibly improve performance by having equal amounts of non-zero
  // elements per thread, instead of as now equal number of columns. Future work.
  std::vector<unsigned int> nnz_pt(nthr,0); // No. of non-zero elements per thread
  for (unsigned int i=0; i<nthr-1; i++) {
    threads[i] = std::thread(&splinefield::no_of_non_zero_elements_sym,this,ncols[i],ncols[i+1],
			     std::ref(csz),std::ref(isz),std::ref(sp),std::ref(nnz_pt[i]));
  }
  this->no_of_non_zero_elements_sym(ncols[nthr-1],ncols[nthr],csz,isz,sp,nnz_pt[nthr-1]);
  std::for_each(threads.begin(),threads.end(),std::mem_fn(&std::thread::join));

  // Convert to start-index per thread
  std::vector<unsigned int> vali(nthr,0);
  for (unsigned int i=1; i<nthr; i++) vali[i] = vali[i-1] + nnz_pt[i-1];

  // Calculate the non-zero elements that need calculating
  for (unsigned int i=0; i<nthr-1; i++) {
    threads[i] = std::thread(&splinefield::calculate_non_zero_elements_sym<T>,this,ncols[i],ncols[i+1],vali[i],
			     std::ref(csz),std::ref(isz),std::ref(sp),std::ref(zeromap),ima,irp,jcp,valp);
  }
  this->calculate_non_zero_elements_sym(ncols[nthr-1],ncols[nthr],vali[nthr-1],
					csz,isz,sp,zeromap,ima,irp,jcp,valp);
  std::for_each(threads.begin(),threads.end(),std::mem_fn(&std::thread::join));

  // Fill in rest of elements using symmetry
  for (unsigned int i=0; i<nthr-1; i++) {
    threads[i] = std::thread(&splinefield::fill_in_non_zero_elements_sym,this,ncols[i],
			     ncols[i+1],std::ref(csz),std::ref(isz),sp,irp,jcp,valp);
  }
  this->fill_in_non_zero_elements_sym(ncols[nthr-1],ncols[nthr],
				      csz,isz,sp,irp,jcp,valp);
  std::for_each(threads.begin(),threads.end(),std::mem_fn(&std::thread::join));

  std::shared_ptr<BFMatrix>   jtj;
  if (prec==BFMatrixFloatPrecision) {
    jtj = std::shared_ptr<BFMatrix>(new SparseBFMatrix<float>(ncoef,ncoef,irp,jcp,valp,Utilities::NoOfThreads(nthr)));
  }
  else jtj = std::shared_ptr<BFMatrix>(new SparseBFMatrix<double>(ncoef,ncoef,irp,jcp,valp,Utilities::NoOfThreads(nthr)));

  delete [] irp; delete [] jcp; delete [] valp;

  return(jtj);
}

void splinefield::no_of_non_zero_elements_sym(// Input
					      unsigned int  first_ci,                // First column index
					      unsigned int  last_ci,                 // One past last column index
					      const std::vector<unsigned int>&  csz, // Size of Coefficient matrix for csp
					      const std::vector<unsigned int>&  isz, // Image size
					      const Spline3D<double>&           sp,  // Spline that determines column/row in JtJ
					      // Output
					      unsigned int& nnz)                     // Total no. of non-zero elements in columns first_ci--last_ci-1
const
{
  nnz = 0;
  std::vector<unsigned int> cindx(3,0);
  for (unsigned int i=first_ci; i<last_ci; i++) {
    cindx[2] = i / (csz[0]*csz[1]);                   // Truncated integer division
    cindx[1] = (i-cindx[2]*csz[0]*csz[1]) / csz[0];   // Truncated integer division
    cindx[0] = i - cindx[2]*csz[0]*csz[1] - cindx[1]*csz[0];
    nnz += sp.NumberOfOverlappingSplines(cindx,isz);
  }
}

template<class T>
void splinefield::calculate_non_zero_elements_sym(// Input
						  unsigned int                      first_ci, // First column index
						  unsigned int                      last_ci,  // One past last column index
						  unsigned int                      vali,     // First index into irp and valp
						  const std::vector<unsigned int>&  csz,      // Size of coefficient matrix
						  const std::vector<unsigned int>&  isz,      // Size of image matrix
						  const Spline3D<double>&           rsp,      // Spline determining row in JtJ
						  const ZeroSplineMap&              zeromap,  // Keeps track of splines for which ima is zero over entire support
						  const T                           *ima,     // Image
						  // Output
						  unsigned int                      *irp,     // nzmax length array of row indicies
						  unsigned int                      *jcp,     // n+1 length array of starts of column indicies
						  double                            *valp)    // nzmax length array of values
const
{
  Spline3D<double>             csp = rsp;   // Read/write copy of read only input spline. Determines column in JtJ
  std::vector<unsigned int>    cindx(3,0);  // Index into coefficient matrix of spline that determines column in JtJ
  std::vector<unsigned int>    rindx(3,0);  // Index into coefficient matrix of spline that determines row in JtJ
  std::vector<unsigned int>    fo(3,0);     // First index of overlapping spline in x-, y- and z-direction
  std::vector<unsigned int>    lo(3,0);     // Last index of overlapping spline in x-, y- and z-direction
  unsigned int                 ci = 0;      // Column index

  for (ci=first_ci; ci<last_ci; ci++) {
    cindx[2] = ci / (csz[0]*csz[1]);                   // Truncated integer division
    cindx[1] = (ci-cindx[2]*csz[0]*csz[1]) / csz[0];   // Truncated integer division
    cindx[0] = ci - cindx[2]*csz[0]*csz[1] - cindx[1]*csz[0];
    jcp[ci] = vali;
    bool cindx_is_zero = zeromap(cindx);
    if (!cindx_is_zero) csp.Premul(cindx,isz,ima);
    csp.RangeOfOverlappingSplines(cindx,isz,fo,lo);
    //
    // Deal with the values above the diagonal.
    // irp is populated and vali is incremented,
    // but valp is set to zero for now.
    //
    for (unsigned int k=fo[2]; k<cindx[2]; k++) {
      for (unsigned int j=fo[1]; j<lo[1]; j++) {
	for (unsigned int i=fo[0]; i<lo[0]; i++) {
	  unsigned int ri = k*csz[1]*csz[0] + j*csz[0] + i;
	  irp[vali] = ri;
	  valp[vali++] = 0;
	}
      }
    }
    for (unsigned int k=cindx[2]; k<lo[2]; k++) {
      //
      // Deal with the values above the "main" diagonals at the
      // 2nd level of symmetry. Again populating irp,
      // incrementing vali and setting valp to zero.
      //
      for (unsigned int j=fo[1]; j<cindx[1]; j++) {
	for (unsigned int i=fo[0]; i<lo[0]; i++) {
	  unsigned int ri = k*csz[1]*csz[0] + j*csz[0] + i;
	  irp[vali] = ri;
	  valp[vali++] = 0;
	}
      }
      for (unsigned int j=cindx[1]; j<lo[1]; j++) {
	//
	// Fill in values above the "main" diagonals at the third
	// and final level of symmetry. Again populating irp,
	// incrementing vali and setting valp to zero.
	//
	for (unsigned int i=fo[0]; i<cindx[0]; i++) {
	  unsigned int ri = k*csz[1]*csz[0] + j*csz[0] + i;
	  irp[vali] = ri;
	  valp[vali++] = 0;
	}
	//
	// And these are the positions for which we actually need to
	// calculate new values. Roughly ~1/8 of the total.
	//
	for (unsigned int i=cindx[0]; i<lo[0]; i++) {
	  unsigned int ri = k*csz[1]*csz[0] + j*csz[0] + i;
	  rindx[0]=i; rindx[1]=j; rindx[2]=k;
	  irp[vali] = ri;
	  if (cindx_is_zero || zeromap(rindx)) valp[vali++] = 0;
	  else valp[vali++] = csp.MulByOther(cindx,rindx,isz,rsp);
	}
      }
    }
  }
  if (ci == csz[2]*csz[1]*csz[0]) jcp[ci] = vali; // If this is the "last" thread
}


void splinefield::fill_in_non_zero_elements_sym(// Input
						unsigned int                      first_ci, // First column index
						unsigned int                      last_ci,  // One past last column index
						const std::vector<unsigned int>&  csz,      // Size of coefficient matrix
						const std::vector<unsigned int>&  isz,      // Size of image matrix
						const Spline3D<double>&           sp,       // Spline
						// Output
						const unsigned int                *irp,     // nzmax length array of row indicies
						const unsigned int                *jcp,     // n+1 length array of starts of column indicies
						double                            *valp)    // nzmax length array of values
const
{
  std::vector<unsigned int>    cindx(3,0);  // Index into coefficient matrix of spline that determines column in JtJ
  std::vector<unsigned int>    fo(3,0);     // First index of overlapping spline in x-, y- and z-direction
  std::vector<unsigned int>    lo(3,0);     // Last index of overlapping spline in x-, y- and z-direction

  for (unsigned int ci=first_ci; ci<last_ci; ci++) { // ci: column index
    cindx[2] = ci / (csz[0]*csz[1]);                   // Truncated integer division
    cindx[1] = (ci-cindx[2]*csz[0]*csz[1]) / csz[0];   // Truncated integer division
    cindx[0] = ci - cindx[2]*csz[0]*csz[1] - cindx[1]*csz[0];
    sp.RangeOfOverlappingSplines(cindx,isz,fo,lo);
    // Loop over all the values we have already calculated
    for (unsigned int k=cindx[2]; k<lo[2]; k++) {
      for (unsigned int j=cindx[1]; j<lo[1]; j++) {
	for (unsigned int i=cindx[0]; i<lo[0]; i++) {
	  int level = static_cast<int>(i!=cindx[0]) + static_cast<int>(j!=cindx[1]) + static_cast<int>(k!=cindx[2]);
	  if (level > 0) { // level==0 would mean we are on main diagonal and no copying should be done
	    unsigned int ri = k*csz[1]*csz[0] + j*csz[0] + i; // row index
	    double val = get_val(ri,ci,irp,jcp,valp);
	    if (level == 1) { // Only one copy across main diagonal
	      set_val(ci,ri,irp,jcp,valp,val);
	    }
	    else if (level == 2) { // One copy across sub-diagonal followed by two copies across main diagonal
	      unsigned int sd_offset = 0; // Sub-diagonal offset
	      if (k != cindx[2]) sd_offset = (k-cindx[2])*csz[0]*csz[1];
	      else if (j != cindx[1]) sd_offset = (j-cindx[1])*csz[0];
	      unsigned int ci2 = ri - sd_offset;
	      unsigned int ri2 = ci + sd_offset;
	      set_val(ri2,ci2,irp,jcp,valp,val);
	      set_val(ci,ri,irp,jcp,valp,val);
	      set_val(ci2,ri2,irp,jcp,valp,val);
	    }
	    else { // One copy across first sub-diagonal, two copies across the second and four across the third.
	      unsigned int sd_offset = (k-cindx[2])*csz[1]*csz[0] + (j-cindx[1])*csz[0]; // First subdiagonal
	      unsigned int ci2 = ri - sd_offset;
	      unsigned int ri2 = ci + sd_offset;
	      set_val(ri2,ci2,irp,jcp,valp,val);
	      sd_offset = (k-cindx[2])*csz[1]*csz[0]; // Second subdiagonal
	      unsigned int ci3 = ri - sd_offset;
	      unsigned int ri3 = ci + sd_offset;
	      unsigned int ci4 = ri2 - sd_offset;
	      unsigned int ri4 = ci2 + sd_offset;
	      set_val(ri3,ci3,irp,jcp,valp,val);
	      set_val(ri4,ci4,irp,jcp,valp,val);
	      // And finally across the main diagonal
	      set_val(ci,ri,irp,jcp,valp,val);
	      set_val(ci2,ri2,irp,jcp,valp,val);
	      set_val(ci3,ri3,irp,jcp,valp,val);
	      set_val(ci4,ri4,irp,jcp,valp,val);
	    }
	  }
	}
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////
//
// Helper routine to set the value for a given row-col in a
// compressed column storage format. Uses bisection. Throws if
// the rwo-col combination doesn't exist.
//
/////////////////////////////////////////////////////////////////////

void splinefield::set_val(unsigned int           row,    // The row we want to find the value of
			  unsigned int           col,    // The column we want to find the value of
			  const unsigned int     *irp,   // Array of row-indicies
			  const unsigned int     *jcp,   // Array of indicies into irp
			  double                 *valp,  // Array of values sorted as irp
			  double                 val)    // The value to set
const
{
   const unsigned int  *a = &(irp[jcp[col]]);
   double              *v = &(valp[jcp[col]]);
   int                 n = jcp[col+1]-jcp[col];
   int                 j = 0;
   int                 jlo = -1;
   int                 jup = n;

   if (row < a[0] || row > a[n-1]) {
     std::ostringstream my_os;
     my_os << "splinefield::set_val: attempt to reference non-existent element at row = " << row << ", col = " << col << std::endl;
     throw BasisfieldException(my_os.str());
   }

   while (jup-jlo > 1)
   {
      j = (jlo+jup) >> 1;
      if (row >= a[j]) {jlo = j;}
      else {jup = j;}
   }

   if (a[jlo] != row) {
     std::ostringstream my_os;
     my_os << "splinefield::set_val: attempt to reference non-existent element at row = " << row << ", col = " << col << std::endl;
     throw BasisfieldException(my_os.str());
   }
   else v[jlo] = val;
}

//
// Calculates 0.5 times the gradient of membrane energy
//
void splinefield::calculate_memen_AtAb(const NEWMAT::ColumnVector&       b,
                                       const std::vector<unsigned int>&  lksp,
                                       const std::vector<unsigned int>&  isz,
                                       const std::vector<double>&        vxs,
                                       const std::vector<unsigned int>&  csz,
                                       unsigned int                      sp_ord,
                                       NEWMAT::ColumnVector&             grad)
const
{
  // Get helper that is the sum of all the component
  // helpers (dd/dx + dd/dy + dd/dy)
  std::vector<unsigned int>  deriv(3,0);
  deriv[0] = 1;
  Spline3D<double>     sp1(sp_ord,lksp,deriv);
  Memen_H_Helper       sum_hlpr(sp1);
  for (unsigned int d=1; d<3; d++) {
    deriv[d-1]=0; deriv[d]=1;
    Spline3D<double>   sp2(sp_ord,lksp,deriv);
    Memen_H_Helper     hlpr(sp2);
    sum_hlpr += hlpr;
  }
  calculate_AtAb(b,isz,csz,sp1,sum_hlpr,grad);
}

//
// Calculates 0.5 times the gradient of bending energy
//
void splinefield::calculate_bender_AtAb(const NEWMAT::ColumnVector&       b,
                                        const std::vector<unsigned int>&  lksp,
                                        const std::vector<unsigned int>&  isz,
                                        const std::vector<double>&        vxs,
                                        const std::vector<unsigned int>&  csz,
                                        unsigned int                      sp_ord,
                                        NEWMAT::ColumnVector&             grad)
const
{
  // Get helper that is the sum of all the component
  // helpers (d2d/dx2 + d2d/dy2 + d2d/dz2 + 2*d2d/dxdy + ...
  Spline3D<double>     sp(sp_ord,lksp);
  Memen_H_Helper       sum_hlpr(sp);
  sum_hlpr *= 0.0;
  for (unsigned int d1=0; d1<3; d1++) {
    for (unsigned int d2=d1; d2<3; d2++) {
      std::vector<unsigned int>   deriv(3,0);
      deriv[d1]++; deriv[d2]++;
      Spline3D<double>   spd(sp_ord,lksp,deriv);   // Spline differentiated in 2 directions
      spd /= (vxs[d1]*vxs[d2]);                         // Derivative in mm^{-1}
      Memen_H_Helper   hlpr(spd);
      if (d1 != d2) hlpr *= 2.0;
      sum_hlpr += hlpr;
    }
  }
  // Use the helper to calculate (A_1^T*A_1 + A_2^T*A_2 ...)*b
  // where A_1 is a matrix with translated d2d/dx2 kernels etc
  calculate_AtAb(b,isz,csz,sp,sum_hlpr,grad);
}

//
// Heart of the calculation of the
// gradient of membrane/bending-energy
//
void splinefield::calculate_AtAb(const NEWMAT::ColumnVector&       b,
                                 const std::vector<unsigned int>&  isz,
                                 const std::vector<unsigned int>&  csz,
                                 const Spline3D<double>&           sp,
                                 const Memen_H_Helper&             hlpr,
                                 NEWMAT::ColumnVector&             grad)
const
{
  grad = 0;
  double *db=static_cast<double *>(b.Store());    // double[] representation of coefficients
  double *dg=static_cast<double *>(grad.Store()); // double[] representation of gradient
  std::vector<unsigned int>   ci(3,0);            // Index of coefficient
  std::vector<unsigned int>   first(3,0);         // index of first overlapping coefficient
  std::vector<unsigned int>   last(3,0);          // Index of last overlapping coefficient
  unsigned int                li=0;               // Linear index of coefficient
  for (ci[2]=0; ci[2]<csz[2]; ci[2]++) {          // First three for loop over all coefficients
    for (ci[1]=0; ci[1]<csz[1]; ci[1]++) {
      for (ci[0]=0; ci[0]<csz[0]; ci[0]++, li++) {
        double *optr = &(dg[li]);
        sp.RangeOfOverlappingSplines(ci,isz,first,last);
        int lastk=last[2]-ci[2]; int lastj=last[1]-ci[1]; int lasti=last[0]-ci[0];
        for (int k=first[2]-ci[2]; k<lastk; k++) {
          int offset = k*csz[1]*csz[0];
          for (int j=first[1]-ci[1]; j<lastj; j++) {
            int offset2 = li + offset + j*csz[0];
            for (int i=first[0]-ci[0]; i<lasti; i++) {
	      *optr += hlpr.Peek(i,j,k) * db[offset2+i];
	      // grad.element(li) += hlpr(i,j,k) * b.element(offset2+i);
	    }
	  }
	}
      }
    }
  }
}

//
// Calculates the contribution to the Hessian from the membrane-energy
// or the bending-energy depending on the parameter et (energy type).
//
void splinefield::calculate_memen_bender_H(// Input
					   const std::vector<unsigned int>&  lksp, // Knot-spacing
					   const std::vector<unsigned int>&  csz,  // Size of coefficent grid
					   const std::vector<unsigned int>&  isz,  // Image size
					   EnergyType                        et,   // Membrane or bending energy
					   // Output. References because they are allocated here.
					   unsigned int*&                    irp,  // Row index vector
					   unsigned int*&                    jcp,  // Vector with start indicies of columns
					   double*&                          valp) // Values
const
{
  // Get Helpers with values for all possible overlaps.
  // For Membrane energy we need 3 helpers, and for
  // bending energy we need 6.
  double vxs[] = {Vxs_x(), Vxs_y(), Vxs_z()};
  std::shared_ptr<Memen_H_Helper>   helpers[6];  // Always room for 6 helpers
  unsigned int nh = 0;                             // Number of helpers
  if (et == MemEn) {
    for (unsigned int d=0; d<3; d++) {
      std::vector<unsigned int>   deriv(3,0);
      deriv[d] = 1;
      Spline3D<double>            spd(_sp.Order(),lksp,deriv);  // Spline differentiated in one direction
      spd /= vxs[d];                                  // Derivative in mm^{-1}
      helpers[nh] = std::shared_ptr<Memen_H_Helper>(new Memen_H_Helper(spd));
      *(helpers[nh++]) *= 2.0;                        // To get factor of 2 of entire matrix.
    }
  }
  else if (et == BendEn) {
    for (unsigned int d1=0; d1<3; d1++) {
      for (unsigned int d2=d1; d2<3; d2++) {
        std::vector<unsigned int>   deriv(3,0);
        deriv[d1]++;
        deriv[d2]++;
        Spline3D<double>            spd(_sp.Order(),lksp,deriv);  // Spline twice differentiated
        spd /= (vxs[d1]*vxs[d2]);                       // Derivative in mm^{-1}
        helpers[nh] = std::shared_ptr<Memen_H_Helper>(new Memen_H_Helper(spd));
        if (d1 == d2) *(helpers[nh++]) *= 2.0;       // To get factor of 2 of entire matrix
        else *(helpers[nh++]) *= 4.0;                // Count cross-terms twice
      }
    }
  }

  // Build compressed column storage representation of H

  std::vector<unsigned int>  cindx(3,0);    // Index of spline
  std::vector<unsigned int>  fo(3,0);       // First index of overlapping spline in x-, y- and z-direction
  std::vector<unsigned int>  lo(3,0);       // Last index of overlapping spline in x-, y- and z-direction

  Spline3D<double>  sp(_sp.Order(),lksp);
  unsigned int      nnz = sp.NzMax(isz);                // Max # of non-zero elements
  irp = new unsigned int[nnz];                          // Row indicies
  jcp = new unsigned int[csz[0]*csz[1]*csz[2]+1];                      // Indicies into irp indicating start/stop of columns
  valp = new double[nnz];                               // The values of the matrix

  unsigned int      vali = 0;                           // Index of present non-zero value (linear indexing)
  unsigned int      ci = 0;                             // Column index

  for (cindx[2]=0; cindx[2]<csz[2]; cindx[2]++) {
    for (cindx[1]=0; cindx[1]<csz[1]; cindx[1]++) {
      for (cindx[0]=0; cindx[0]<csz[0]; cindx[0]++) {
        ci = cindx[2]*csz[1]*csz[0] + cindx[1]*csz[0] + cindx[0];
        jcp[ci] = vali;
        sp.RangeOfOverlappingSplines(cindx,isz,fo,lo);
        for (unsigned int k=fo[2]; k<lo[2]; k++) {
          for (unsigned int j=fo[1]; j<lo[1]; j++) {
	    unsigned int bi = k*csz[1]*csz[0] + j*csz[0];
            for (unsigned int i=fo[0]; i<lo[0]; i++) {
              irp[vali] = bi+i;
              valp[vali] = 0.0;
              for (unsigned int d=0; d<nh; d++) {
                valp[vali] += helpers[d]->Peek(i-cindx[0],j-cindx[1],k-cindx[2]);
	      }
              vali++;
	    }
	  }
	}
      }
    }
  }
  jcp[ci+1] = vali;

  return;
}

//
// Calculates the contribution to the Hessian from the membrane-energy
// or the bending-energy depending on the parameter et (energy type).
// Does the same thing as the function above, but parallelised.
//
void splinefield::calculate_memen_bender_H_para(// Input
						const std::vector<unsigned int>&  lksp, // Knot-spacing
						const std::vector<unsigned int>&  csz,  // Size of coefficent grid
						const std::vector<unsigned int>&  isz,  // Image size
						EnergyType                        et,   // Membrane or bending energy
						unsigned int                      nthr, // Number of threads
						// Output. References because they are allocated here.
						unsigned int*&                    irp,  // Row index vector
						unsigned int*&                    jcp,  // Vector with start indicies of columns
						double*&                          valp) // Values
const
{
  // Get Helpers with values for all possible overlaps.
  // For Membrane energy we need 3 helpers, and for
  // bending energy we need 6.
  double vxs[] = {Vxs_x(), Vxs_y(), Vxs_z()};
  std::shared_ptr<Memen_H_Helper>   helpers[6];  // Always room for 6 helpers
  unsigned int nh = 0;                             // Number of helpers
  if (et == MemEn) {
    for (unsigned int d=0; d<3; d++) {
      std::vector<unsigned int>   deriv(3,0);
      deriv[d] = 1;
      Spline3D<double>            spd(_sp.Order(),lksp,deriv);  // Spline differentiated in one direction
      spd /= vxs[d];                                  // Derivative in mm^{-1}
      helpers[nh] = std::shared_ptr<Memen_H_Helper>(new Memen_H_Helper(spd));
      *(helpers[nh++]) *= 2.0;                        // To get factor of 2 of entire matrix.
    }
  }
  else if (et == BendEn) {
    for (unsigned int d1=0; d1<3; d1++) {
      for (unsigned int d2=d1; d2<3; d2++) {
        std::vector<unsigned int>   deriv(3,0);
        deriv[d1]++;
        deriv[d2]++;
        Spline3D<double>            spd(_sp.Order(),lksp,deriv);  // Spline twice differentiated
        spd /= (vxs[d1]*vxs[d2]);                       // Derivative in mm^{-1}
        helpers[nh] = std::shared_ptr<Memen_H_Helper>(new Memen_H_Helper(spd));
        if (d1 == d2) *(helpers[nh++]) *= 2.0;       // To get factor of 2 of entire matrix
        else *(helpers[nh++]) *= 4.0;                // Count cross-terms twice
      }
    }
  }

  Spline3D<double>  sp(_sp.Order(),lksp);
  unsigned int      nnz = sp.NzMax(isz);                // Max # of non-zero elements
  irp = new unsigned int[nnz];                          // Row indicies
  jcp = new unsigned int[csz[0]*csz[1]*csz[2]+1];       // Indicies into irp indicating start/stop of columns
  valp = new double[nnz];                               // The values of the matrix

  std::vector<std::thread> threads(nthr-1);
  // Determine range of columns per thread
  double cpt = static_cast<double>(csz[2]*csz[1]*csz[0]) / static_cast<double>(nthr);
  std::vector<unsigned int> ncols(nthr+1,0); //
  for (unsigned int i=1; i<nthr; i++) ncols[i] = static_cast<unsigned int>(std::floor(i*cpt));
  ncols[nthr] = csz[2]*csz[1]*csz[0];

  // Find number of elements calculated by each thread
  std::vector<unsigned int> nnz_pt(nthr,0); // No. of non-zero elements per thread
  for (unsigned int i=0; i<nthr-1; i++) {
    threads[i] = std::thread(&splinefield::calculate_memen_bender_hess_no_of_elements,this,ncols[i],
			     ncols[i+1],std::ref(csz),std::ref(isz),std::ref(sp),std::ref(nnz_pt[i]));
  }
  this->calculate_memen_bender_hess_no_of_elements(ncols[nthr-1],ncols[nthr],csz,isz,sp,nnz_pt[nthr-1]);
  std::for_each(threads.begin(),threads.end(),std::mem_fn(&std::thread::join));

  // Convert to start-index per thread
  std::vector<unsigned int> vali(nthr,0);
  for (unsigned int i=1; i<nthr; i++) vali[i] = vali[i-1] + nnz_pt[i-1];

  // Finally, do the calculations

  for (unsigned int i=0; i<nthr-1; i++) {
    threads[i] = std::thread(&splinefield::calculate_memen_bender_hess_columns,this,ncols[i],ncols[i+1],
			     vali[i],std::ref(csz),std::ref(isz),std::ref(sp),nh,helpers,irp,jcp,valp);
  }
  this->calculate_memen_bender_hess_columns(ncols[nthr-1],ncols[nthr],vali[nthr-1],csz,isz,sp,nh,helpers,irp,jcp,valp);
  std::for_each(threads.begin(),threads.end(),std::mem_fn(&std::thread::join));

  return;
}

void splinefield::calculate_memen_bender_hess_no_of_elements(//Input
                                                             unsigned int                      first_ci,  // First column index
							     unsigned int                      last_ci,   // One past last column index
							     const std::vector<unsigned int>&  csz,       // Size of coefficient grid
							     const std::vector<unsigned int>&  isz,       // Size of image
							     const Spline3D<double>&           sp,        // Spline
							     // Output
							     unsigned int&                     nnz)       // Number of non-zero elements
const
{
  nnz = 0;
  std::vector<unsigned int> cindx(3,0);
  for (unsigned int ci=first_ci; ci<last_ci; ci++) {
    cindx[2] = ci / (csz[0]*csz[1]);                   // Truncated integer division
    cindx[1] = (ci-cindx[2]*csz[0]*csz[1]) / csz[0];   // Truncated integer division
    cindx[0] = ci - cindx[2]*csz[0]*csz[1] - cindx[1]*csz[0];
    nnz += sp.NumberOfOverlappingSplines(cindx,isz);
  }
}

void splinefield::calculate_memen_bender_hess_columns(// Input
						      unsigned int                      first_ci,  // First column index
						      unsigned int                      last_ci,   // One past last column index
						      unsigned int                      vali,      // Starting index into irp and valp
						      const std::vector<unsigned int>&  csz,       // Size of coefficient grid
						      const std::vector<unsigned int>&  isz,       // Size of image
						      const Spline3D<double>&           sp,        // Spline
						      unsigned int                      nh,        // Number of helpers
						      std::shared_ptr<Memen_H_Helper>  *helpers,  // Array of helpers
						      // Output
						      unsigned int                     *irp,      // Row index vector
						      unsigned int                     *jcp,      // Vector with start indicies of columns
						      double                           *valp)     // Values
const
{
  std::vector<unsigned int>  fo(3,0);       // First index of overlapping spline in x-, y- and z-direction
  std::vector<unsigned int>  lo(3,0);       // Last index of overlapping spline in x-, y- and z-direction
  unsigned int ci = 0;                      // Column index
  for (ci=first_ci; ci<last_ci; ci++) {
    std::vector<unsigned int> cindx(3,0);
    cindx[2] = ci / (csz[0]*csz[1]);                   // Truncated integer division
    cindx[1] = (ci-cindx[2]*csz[0]*csz[1]) / csz[0];   // Truncated integer division
    cindx[0] = ci - cindx[2]*csz[0]*csz[1] - cindx[1]*csz[0];
    jcp[ci] = vali;
    sp.RangeOfOverlappingSplines(cindx,isz,fo,lo);
    for (unsigned int k=fo[2]; k<lo[2]; k++) {
      for (unsigned int j=fo[1]; j<lo[1]; j++) {
	unsigned int bi = k*csz[1]*csz[0] + j*csz[0];
	for (unsigned int i=fo[0]; i<lo[0]; i++) {
	  irp[vali] = bi+i;
	  valp[vali] = 0.0;
	  for (unsigned int d=0; d<nh; d++) {
	    valp[vali] += helpers[d]->Peek(i-cindx[0],j-cindx[1],k-cindx[2]);
	  }
	  vali++;
	}
      }
    }
  }
  if (ci == csz[2]*csz[1]*csz[0]) jcp[ci] = vali; // If this is the "last" thread
}

void splinefield::hadamard(const NEWIMAGE::volume<float>& ima1,
                           const NEWIMAGE::volume<float>& ima2,
                           float                          *prod) const
{
  if (!samesize(ima1,ima2,3,true)) throw BasisfieldException("hadamard: Image dimension mismatch");

  for (NEWIMAGE::volume<float>::fast_const_iterator it1=ima1.fbegin(), it2=ima2.fbegin(), it1_end=ima1.fend(); it1 != it1_end; ++it1, ++it2, ++prod) {
    *prod = (*it1) * (*it2);
  }
}

void splinefield::hadamard(const NEWIMAGE::volume<float>& ima1,
                           const NEWIMAGE::volume<float>& ima2,
                           const NEWIMAGE::volume<char>&  mask,
                           float                          *prod) const
{
  if (!samesize(ima1,ima2,3,true) || !samesize(ima1,mask,3)) throw BasisfieldException("hadamard: Image dimension mismatch");

  NEWIMAGE::volume<char>::fast_const_iterator itm = mask.fbegin();
  for (NEWIMAGE::volume<float>::fast_const_iterator it1=ima1.fbegin(), it2=ima2.fbegin(), it1_end=ima1.fend(); it1 != it1_end; ++it1, ++it2, ++itm, ++prod) {
    *prod = static_cast<float>(*itm) * (*it1) * (*it2);
  }
}

void splinefield::hadamard(const NEWIMAGE::volume<float>& ima1,
                           const NEWIMAGE::volume<float>& ima2,
                           const NEWIMAGE::volume<char>   *mask,
                           float                          *prod) const
{
  if (mask) hadamard(ima1,ima2,*mask,prod);
  else hadamard(ima1,ima2,prod);
}

void splinefield::copy_and_mask_ima(const NEWIMAGE::volume<float>& ima,
				    const NEWIMAGE::volume<char>   *mask,
				    float                          *mima) const
{
  if (mask != NULL) {
    NEWIMAGE::volume<char>::fast_const_iterator itm = mask->fbegin();
    for (NEWIMAGE::volume<float>::fast_const_iterator iti=ima.fbegin(), iti_end=ima.fend(); iti != iti_end; ++iti, ++itm, ++mima) *mima = static_cast<float>(*itm) * (*iti);
  }
  else {
    for (NEWIMAGE::volume<float>::fast_const_iterator iti=ima.fbegin(), iti_end=ima.fend(); iti != iti_end; ++iti, ++mima) *mima = *iti;
  }
}

void splinefield::Set(const NEWIMAGE::volume<float>& pfield)
{
  if (int(FieldSz_x()) != pfield.xsize() || int(FieldSz_y()) != pfield.ysize() || int(FieldSz_z()) != pfield.zsize()) {
    throw BasisfieldException("basisfield::Set:: Matrix size mismatch beween basisfield class and supplied field");
  }
  if (Vxs_x() != pfield.xdim() || Vxs_y() != pfield.ydim() || Vxs_z() != pfield.zdim()) {
    throw BasisfieldException("basisfield::Set:: Voxel size mismatch beween basisfield class and supplied field");
  }

  helper_vol orig(pfield);
  helper_vol conv = get_coefs_one_dim(orig,0,CoefSz_x(),Order(),Ksp_x());
  if (NDim() > 1) conv = get_coefs_one_dim(conv,1,CoefSz_y(),Order(),Ksp_y());
  if (NDim() > 2) conv = get_coefs_one_dim(conv,2,CoefSz_z(),Order(),Ksp_z());

  NEWMAT::ColumnVector lcoef = conv.AsNewmatVector();
  SetCoef(lcoef);
  return;
}


//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// The following is a set of routines that are used for zooming
// fields, and for deciding which zooms are valid and which are
// not. These will be almost excessively commented. The reason
// for that is that I have struggled so badly to get things clear
// in my own head, and I don't want to return in 6 months not
// understanding the code.
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

/////////////////////////////////////////////////////////////////////
//
// The following is the "main" zooming routine. It will return a
// completely new field, with new matrix size and/or voxel size
// and/or knot-spacing. The new field will have coefficients set
// so that the field is identical to the initial field (in the
// case of upsampling) or the "best field in a least squares sense"
// (in the case of downsampling).
//
/////////////////////////////////////////////////////////////////////

std::shared_ptr<BASISFIELD::basisfield> splinefield::ZoomField(const std::vector<unsigned int>&   nms,
                                                                 const std::vector<double>&         nvxs,
                                                                 std::vector<unsigned int>          nksp) const
{
  std::vector<unsigned int> oms(3), oksp(3);
  std::vector<double>       ovxs(3);
  oms[0] = FieldSz_x(); oms[1] = FieldSz_y(); oms[2] = FieldSz_z();
  oksp[0] = Ksp_x(); oksp[1] = Ksp_y(); oksp[2] = Ksp_z();
  ovxs[0] = Vxs_x(); ovxs[1] = Vxs_y(); ovxs[2] = Vxs_z();

  if (!nksp.size()) nksp = oksp;

  // cout << "Old matrix size: " << oms[0] << "  " << oms[1] << "  " << oms[2] << endl;
  // cout << "new matrix size: " << nms[0] << "  " << nms[1] << "  " << nms[2] << endl;
  // cout << "Old voxel size: " << ovxs[0] << "  " << ovxs[1] << "  " << ovxs[2] << endl;
  // cout << "New voxel size: " << nvxs[0] << "  " << nvxs[1] << "  " << nvxs[2] << endl;
  // cout << "Old knot-spacing: " << oksp[0] << "  " << oksp[1] << "  " << oksp[2] << endl;
  // cout << "New knot-spacing: " << nksp[0] << "  " << nksp[1] << "  " << nksp[2] << endl;

  // Check that new field is kosher
  // Make sure that we are not asked to change both voxel-size and knot-spacing
  if (nksp != oksp && nvxs != ovxs) throw BasisfieldException("ZoomField: Cannot change both voxel-size and knot-spacing");
  // If voxel-size changed, make sure that the new voxel-size allows us to estimate the new coefficients
  if (nvxs != ovxs && !new_vxs_is_ok(nvxs)) throw BasisfieldException("ZoomField: The requested voxel-size is invalid");

  // If voxel size changed, fake an old knot-spacing for the purpose of estimating new coefficients
  // This will always work when going from low->higher resolution (in which case the fake knot-spacing
  // will be some integer multiple of the old knot-spacing). It will not always work when going from
  // high->lower resolution since we may then get a fake knot-spacing that should really be a
  // non-integer #, which or present implementation cannot handle. This is really indicative of a
  // flawed implementation, but for the time being I will just work around it.

  std::vector<unsigned int> fksp = oksp;
  if (nvxs != ovxs) {
    if (faking_works(nvxs,nksp,ovxs)) {
      fksp = fake_old_ksp(nvxs,nksp,ovxs);
    }
    else {
      return(zoom_field_in_stupid_way(nms,nvxs,nksp));
    }
  }
  // cout << "fksp[0] = " << fksp[0] << ", fksp[1] = " << fksp[1] << ", fksp[2] = " << fksp[2] << endl;

  // Create the new field
  std::shared_ptr<BASISFIELD::splinefield> tptr(new BASISFIELD::splinefield(nms,nvxs,nksp,this->Order()));

  // Get original set of coefficients
  const std::shared_ptr<NEWMAT::ColumnVector> ocoef = GetCoef();

  // If not all coefficients zero, create the coefficients for the
  // new field from those of the old field
  NEWMAT::ColumnVector zerovec(ocoef->Nrows());
  zerovec = 0.0;

  if (*ocoef != zerovec) {

    // Repack coefficient sizes of new and old field (for convenience)
    std::vector<unsigned int> ocs(3);
    ocs[0]=CoefSz_x(); ocs[1]=CoefSz_y(); ocs[2]=CoefSz_z();
    std::vector<unsigned int> ncs(3);
    ncs[0]=tptr->CoefSz_x(); ncs[1]=tptr->CoefSz_y(); ncs[2]=tptr->CoefSz_z();

    // Resample x-direction
    BASISFIELD::Spline1D<double> osp(3,fksp[0]);                      // Spline object for old spline
    BASISFIELD::Spline1D<double> nsp(3,nksp[0]);                      // Spline object for new spline
    NEWMAT::Matrix M = nsp.GetMMatrix(osp,nms[0],ncs[0],ocs[0]);      // Resampling matrix
    double *tmp_coef_x = new double[ncs[0]*ocs[1]*ocs[2]];            // Temporary coefficient matrix
    NEWMAT::ColumnVector iv(ocs[0]);                                  // Vector holding one "column" running in x-direction
    NEWMAT::ColumnVector ov(ncs[0]);                                  // Dito, after resampling
    for (unsigned int k=0; k<ocs[2]; k++) {
      for (unsigned int j=0; j<ocs[1]; j++) {
        for (unsigned int i=0; i<ocs[0]; i++) {
          iv.element(i) = ocoef->element(k*ocs[0]*ocs[1]+j*ocs[0]+i);  // Collect old column
        }
        ov = M*iv;                                                    // Calculate new column
        for (unsigned int i=0; i<ncs[0]; i++) {
	  tmp_coef_x[k*ncs[0]*ocs[1]+j*ncs[0]+i] = ov.element(i);     // Put it into temporary volume
        }
      }
    }

    // Resample y-direction
    osp = BASISFIELD::Spline1D<double>(3,fksp[1]);                    // Spline object for old spline
    nsp = BASISFIELD::Spline1D<double>(3,nksp[1]);                    // Spline object for new spline
    M = nsp.GetMMatrix(osp,nms[1],ncs[1],ocs[1]);                     // Resampling matrix
    double *tmp_coef_y = new double[ncs[0]*ncs[1]*ocs[2]];            // Temporary coefficient matrix
    iv.ReSize(ocs[1]);                                                // Vector holding one "column" running in y-direction
    ov.ReSize(ncs[1]);                                                // Dito, after resampling
    for (unsigned int k=0; k<ocs[2]; k++) {
      for (unsigned int i=0; i<ncs[0]; i++) {
        for (unsigned int j=0; j<ocs[1]; j++) {
          iv.element(j) = tmp_coef_x[k*ncs[0]*ocs[1]+j*ncs[0]+i];     // Collect old column
        }
        ov = M*iv;                                                    // Calculate new column
        for (unsigned int j=0; j<ncs[1]; j++) {
          tmp_coef_y[k*ncs[0]*ncs[1]+j*ncs[0]+i] = ov.element(j);     // Put it into temporary volume
        }
      }
    }
    delete[] tmp_coef_x;

    // Resample z-direction
    osp = BASISFIELD::Spline1D<double>(3,fksp[2]);                    // Spline object for old spline
    nsp = BASISFIELD::Spline1D<double>(3,nksp[2]);                    // Spline object for new spline
    M = nsp.GetMMatrix(osp,nms[2],ncs[2],ocs[2]);                     // Resampling matrix
    NEWMAT::ColumnVector ncoef(ncs[0]*ncs[1]*ncs[2]);                 // Temporary coefficient matrix, now as ColumnVector
    iv.ReSize(ocs[2]);                                                // Vector holding one "column" running in z-direction
    ov.ReSize(ncs[2]);                                                // Dito, after resampling
    for (unsigned int j=0; j<ncs[1]; j++) {
      for (unsigned int i=0; i<ncs[0]; i++) {
        for (unsigned int k=0; k<ocs[2]; k++) {
          iv.element(k) = tmp_coef_y[k*ncs[0]*ncs[1]+j*ncs[0]+i];     // Collect old column
        }
        ov = M*iv;                                                    // Calculate new column.
        for (unsigned int k=0; k<ncs[2]; k++) {
          ncoef.element(k*ncs[0]*ncs[1]+j*ncs[0]+i) = ov.element(k);  // Put into final set of coefficients
        }
      }
    }
    delete[] tmp_coef_y;

    tptr->SetCoef(ncoef);
  }

  return(tptr);
}

/*
void splinefield::SetToConstant(double fv)
{
  NEWMAT::ColumnVector  lcoef(CoefSz());
  lcoef = fv;
  SetCoef(lcoef);
}
*/

void splinefield::SetToConstant(double fv)
{
  vector<unsigned int>  sz(3,0);
  sz[0] = FieldSz_x(); sz[1] = FieldSz_y(); sz[2] = FieldSz_z();
  vector<double>        vxs(3,0.0);
  vxs[0] = Vxs_x(); vxs[1] = Vxs_y(); vxs[2] = Vxs_z();
  NEWIMAGE::volume<float>  vol(sz[0],sz[1],sz[2]);
  vol.setdims(vxs[0],vxs[1],vxs[2]);
  vol = fv;

  Set(vol);
}

/////////////////////////////////////////////////////////////////////
//
// Returns the new matrix size for a given level of subsampling
// provided that this is a power of two. It is done by calling
// "next_size_down" recursively so that at each subsequent upsampling
// step the previous field should have a slightly larger FOV than
// the new one. This guarantees that as we go to higher resolutions
// we will not be in a position of having to extrapolate from a
// previous resolution.
//
// The "power of 2 constraint" is not strictly necessary, and the
// routines for doing the actual zooming has less severe constraints.
// However, I have decided my life is a bit easier if I enforce this
// constraint at this level.
//
/////////////////////////////////////////////////////////////////////

std::vector<unsigned int> splinefield::SubsampledMatrixSize(const std::vector<unsigned int>&  ss,        // Subsampling factor
                                                            std::vector<unsigned int>         oms) const // Old Matrix Size
{
  std::vector<unsigned int>   nms;  // New matrix size
  if (!oms.size()) {nms.resize(NDim()); nms[0]=FieldSz_x(); if (NDim()>1) nms[1]=FieldSz_y(); if (NDim()>2) nms[2]=FieldSz_z();}
  else nms = oms;
  if (nms.size() != ss.size()) throw BasisfieldException("splinefield::SubsampledMatrixSize: Size mismatch between ss and oms");
  for (unsigned int i=0; i<ss.size(); i++) nms[i]=subsampled_matrix_size(ss[i],nms[i]);

  return(nms);
}

unsigned int splinefield::subsampled_matrix_size(unsigned int  ss,        // Subsampling factor
                                                 unsigned int  oms) const // Old matrix size
{
  if (!is_a_power_of_2(ss)) throw BasisfieldException("splinefield::SubsampledMatrixSize: Subsampling factor not a power of 2");
  while (ss > 1) {
    oms = next_size_down(oms);
    ss /= 2;
  }
  return(oms);
}

/////////////////////////////////////////////////////////////////////
//
// Returns the new voxel-size for a given level of subsampling.
// For splinefields this is simply the old voxel-size divided by
// the subsampling factor, guaranteeing that every voxel centre
// in the original (low res) field is represented by a voxel centre
// in the new field.
// For DCT-fields it is a little less straightforward, which is
// the reason we have declared the routine in the way it is.
//
/////////////////////////////////////////////////////////////////////

std::vector<double> splinefield::SubsampledVoxelSize(const std::vector<unsigned int>&  ss,       // Subsampling factor
                                                     std::vector<double>               ovxs,     // Old voxel size
                                                     std::vector<unsigned int>         ms) const // Matrix size
{
  std::vector<double> nvxs;
  if (!ovxs.size()) {nvxs.resize(NDim()); nvxs[0]=Vxs_x(); if (NDim()>1) nvxs[1]=Vxs_y(); if (NDim()>2) nvxs[2]=Vxs_z();}
  else nvxs = ovxs;
  if (ovxs.size() != ss.size()) throw BasisfieldException("splinefield::SubsampledVoxelSize: Size mismatch between ss and ovxs");
  for (unsigned int i=0; i<ss.size(); i++) nvxs[i]=subsampled_voxel_size(ss[i],nvxs[i]);

  return(nvxs);
}

double splinefield::subsampled_voxel_size(unsigned int   ss,         // Subsampling factor
                                          double         ovxs) const // Old voxel size
{
  if (!is_a_power_of_2(ss)) throw BasisfieldException("splinefield::subsampled_voxel_size: Subsampling factor not a power of 2");
  return(double(ss)*ovxs);
}
/////////////////////////////////////////////////////////////////////
//
// Returns the "next size down" in a sampling pyramide where at each
// stage subsampling is done with a factor of two. The guiding principle
// is that int the new (downsampled) field the centre of the first
// voxel should coincide with the centre of the first voxel in the
// previous field. This means that the edge of that first voxel will
// extend beyond the edge of the previous field. At the other end
// the last voxel will extend a similar amount, or by an additional
// voxel depending on if original size is odd or even.
//
/////////////////////////////////////////////////////////////////////

std::vector<unsigned int> splinefield::next_size_down(const std::vector<unsigned int>& isize) const
{
  std::vector<unsigned int>  osize(isize.size(),0);
  for (int i=0; i<int(isize.size()); i++) osize[i] = next_size_down(isize[i]);
  return(osize);
}

unsigned int splinefield::next_size_down(unsigned int isize) const
{
  if (isize%2) return((isize+1)/2);  // if odd
  else return(isize/2 + 1);          // if even
}

/////////////////////////////////////////////////////////////////////
//
// Routine to make sure that a given subsampling factor is a power
// of two.
//
/////////////////////////////////////////////////////////////////////

bool splinefield::is_a_power_of_2(double fac) const
{
  double candidates[] = {1.0/32.0, 1.0/16.0, 1.0/8.0, 1.0/4.0, 1.0/2.0, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0};
  double eps = 1.0e-16;

  for (unsigned int i=0; i<sizeof(candidates)/sizeof(candidates[0]); i++) {
    if (fabs(fac-candidates[i])<eps) return(true);
  }
  return(false);
}

bool splinefield::is_a_power_of_2(unsigned int fac) const
{
  if (fac >= 1) return(is_a_power_of_2(double(fac)));
  else return(false);
}

bool splinefield::are_a_power_of_2(const std::vector<double>& facs) const
{
  bool retval = true;
  for (unsigned int i=0; i<facs.size(); i++) if (!is_a_power_of_2(facs[i])) retval = false;
  return(retval);
}

bool splinefield::are_a_power_of_2(const std::vector<unsigned int>& facs) const
{
  bool retval = true;
  for (int unsigned i=0; i<facs.size(); i++) if (!is_a_power_of_2(facs[i])) retval = false;
  return(retval);
}

/////////////////////////////////////////////////////////////////////
//
// When zooming a field from one voxel-size to another we can only
// calculate the coefficients for the new field if there are some
// set of shared voxel-centres between the two representations.
// This will only be the case if the new voxel size is some integer
// factor of the old size (for downsampling) or a fraction
// (1.0/n, where n is integer) of the old size (for upsampling).
// The routines below ensure that is the case.
//
/////////////////////////////////////////////////////////////////////

bool splinefield::new_vxs_is_ok(const std::vector<double>& nvxs,
                                std::vector<double>        ovxs) const
{
  if (!ovxs.size()) {ovxs.resize(3); ovxs[0]=Vxs_x(); ovxs[1]=Vxs_y(); ovxs[2]=Vxs_z();}
  if (ovxs.size() != nvxs.size()) throw BasisfieldException("splinefield::new_vxs_is_ok: size mismatch between nvxs and ovxs");

  for (unsigned int i=0; i<nvxs.size(); i++) if (!new_vxs_is_ok(nvxs[i],ovxs[i])) return(false);
  return(true);
}
bool splinefield::new_vxs_is_ok(double nvxs,
                                double ovxs) const
{
  double eps = 1.0e-16;
  if (nvxs/ovxs < 1.0) {
    for (int i=32; i>1; i--) if (fabs((1.0/double(i))-(nvxs/ovxs)) < eps) return(true);
  }
  else if (fabs((nvxs/ovxs)-1.0) < eps) return(true);
  else {
    for (int i=2; i<33; i++) if (fabs(double(i)-(nvxs/ovxs)) < eps) return(true);
  }
  return(false);
}

/////////////////////////////////////////////////////////////////////
//
// These routines will provide a fake knot-spacing. When talking
// about "zooming" fields we may consider going from one voxel-size
// to another, or we may think of changing our parametrisation from
// one knot-spacing to another. In our main routine we treat these
// cases in an equivalent way by "transforming" the case where we
// go from one voxel size to another to a case where we change
// knot-spacing. For example if we have a knot-spacing of 3voxels
// and a voxel-size of 4mm and want to go to a knot-spacing of
// 3voxels for a 2mm we will "pretend" that we are in fact going
// from a knot-spacing of 6voxels to 3voxels in a 2mm voxel matrix.
//
/////////////////////////////////////////////////////////////////////

std::vector<unsigned int> splinefield::fake_old_ksp(const std::vector<double>&        nvxs,
                                                    const std::vector<unsigned int>&  nksp,
                                                    std::vector<double>               ovxs) const
{
  if (!ovxs.size()) {ovxs.resize(NDim()); ovxs[0]=Vxs_x(); if (NDim()>1) ovxs[1]=Vxs_y(); if (NDim()>2) ovxs[2]=Vxs_z();}
  if (ovxs.size()!=nvxs.size() || ovxs.size()!=nksp.size()) throw BasisfieldException("splinefield::fake_old_ksp: size mismatch between nvxs, ovxs and nksp");

  std::vector<unsigned int>   fksp(ovxs.size());
  for (unsigned int i=0; i<ovxs.size(); i++) fksp[i] = fake_old_ksp(nvxs[i],nksp[i],ovxs[i]);
  return(fksp);
}

unsigned int splinefield::fake_old_ksp(double        nvxs,
                                       unsigned int  nksp,
                                       double        ovxs) const
{
  return(static_cast<unsigned int>(MISCMATHS::round((ovxs/nvxs)*double(nksp))));
}

bool splinefield::faking_works(const std::vector<double>&        nvxs,
                               const std::vector<unsigned int>&  nksp,
                               std::vector<double>               ovxs) const
{
  if (!ovxs.size()) {ovxs.resize(NDim()); ovxs[0]=Vxs_x(); if (NDim()>1) ovxs[1]=Vxs_y(); if (NDim()>2) ovxs[2]=Vxs_z();}
  for (unsigned int i=0; i<ovxs.size(); i++) if (!faking_works(nvxs[i],nksp[i],ovxs[i])) return(false);
  return(true);
}

bool splinefield::faking_works(double        nvxs,
                               unsigned int  nksp,
                               double        ovxs) const
{
  double       dres = (ovxs/nvxs)*double(nksp);
  unsigned int ires = static_cast<unsigned int>(roundl(dres));
  if (fabs(double(ires)-dres) > 1e-6) return(false);
  return(true);
}

std::shared_ptr<BASISFIELD::basisfield> splinefield::zoom_field_in_stupid_way(const std::vector<unsigned int>&    nms,
                                                                                const std::vector<double>&          nvxs,
                                                                                const std::vector<unsigned int>&    nksp) const
{
  // Create new field
  std::shared_ptr<BASISFIELD::splinefield> tptr(new BASISFIELD::splinefield(nms,nvxs,nksp,this->Order()));

  // See if we have non-zero field that we need to resample.
  const std::shared_ptr<NEWMAT::ColumnVector> ocoef = GetCoef();
  NEWMAT::ColumnVector zerovec(ocoef->Nrows());
  zerovec = 0.0;
  if (*ocoef == zerovec) return(tptr);

  // Get image representation of current field in resolution of new field
  NEWIMAGE::volume<float>   nvol(nms[0],nms[1],nms[2]);
  nvol.setdims(nvxs[0],nvxs[1],nvxs[2]);
  double z=0.0;
  for (unsigned int k=0; k<nms[2] && z<double(FieldSz_z()); k++, z+=(nvxs[2]/Vxs_z())) {
    double y=0.0;
    for (unsigned int j=0; j<nms[1] && y<double(FieldSz_y()); j++, y+=(nvxs[1]/Vxs_y())) {
      double x=0.0;
      for (unsigned int i=0; i<nms[0] && x<double(FieldSz_x()); i++, x+=(nvxs[0]/Vxs_x())) {
        nvol(i,j,k) = Peek(x,y,z);
      }
    }
  }

  // Set this for new field, which means the coefficients will be calculated.
  tptr->Set(nvol);

  return(tptr);
}

helper_vol splinefield::get_coefs_one_dim(const helper_vol&                 in,
                                          unsigned int                      dim,
                                          unsigned int                      csz,
                                          unsigned int                      order,
                                          unsigned int                      ksp) const
{
  std::vector<unsigned int> nsz(3,0);
  for (unsigned int i=0; i<3; i++) {if (i==dim) nsz[i]=csz; else nsz[i]=in.Size(i); }
  helper_vol  out(nsz);

  unsigned int ii, jj;
  switch (dim) {
  case 0:
    ii=in.Size(1);
    jj=in.Size(2);
    break;
  case 1:
    ii=in.Size(0);
    jj=in.Size(2);
    break;
  case 2:
    ii=in.Size(0);
    jj=in.Size(1);
    break;
  default:
    throw ;
  }

  Spline1D<double>  sp(order,ksp);
  NEWMAT::Matrix A = sp.GetAMatrix(in.Size(dim),csz);
  NEWMAT::Matrix S = get_s_matrix(in.Size(dim),csz);
  NEWMAT::Matrix AS = A & 0.005*S;
  NEWMAT::Matrix M = (AS.t()*AS).i()*A.t();
  NEWMAT::ColumnVector y(in.Size(dim));
  for (unsigned int i=0; i<ii; i++) {
    for (unsigned int j=0; j<jj; j++) {
      in.GetColumn(i,j,dim,static_cast<double *>(y.Store()));
      NEWMAT::ColumnVector b = M*y;
      out.SetColumn(i,j,dim,static_cast<double *>(b.Store()));
    }
  }
  return(out);
}

NEWMAT::ReturnMatrix splinefield::get_s_matrix(unsigned int isz,
                                               unsigned int csz) const
{
  NEWMAT::Matrix S(isz,csz);
  S = 0.0;
  S(1,1) = 2.0; S(1,2) = -1.0; S(1,csz) = -1.0;
  S(isz,1) = -1.0; S(isz,csz-1) = -1.0; S(isz,csz) = 2.0;

  S.Release();
  return(S);
}

// This routine will return a value for the field (or a derivative of the field)
// from outside the "valid" FOV. Since the splines actually extend a bit outside
// the FOV there will be a gradual taper off to zero.
double splinefield::peek_outside_fov(int i, int j, int k, FieldIndex fi) const
{
  std::vector<double>  vox(3,0.0);
  vox[0] = static_cast<double>(i); vox[1] = static_cast<double>(j); vox[2] = static_cast<double>(k);
  std::vector<unsigned int> csz(3,0);
  csz[0] = CoefSz_x(); csz[1] = CoefSz_y(); csz[2] = CoefSz_z();
  std::vector<unsigned int> cindx(3,0);
  std::vector<unsigned int> first(3,0), last(3,0);

  double rval = 0.0;
  if (fi == FIELD) {
    _sp.RangeOfSplines(vox,csz,first,last);
    // cout << "vox = " << vox[0] << ",  " << vox[1] << ",  " << vox[2] << endl;
    // cout << "csz = " << csz[0] << ",  " << csz[1] << ",  " << csz[2] << endl;
    // cout << "first = " << first[0] << ",  " << first[1] << ",  " << first[2] << endl;
    // cout << "last = " << last[0] << ",  " << last[1] << ",  " << last[2] << endl;
    for (cindx[2]=first[2]; cindx[2]<last[2]; cindx[2]++) {
      for (cindx[1]=first[1]; cindx[1]<last[1]; cindx[1]++) {
        for (cindx[0]=first[0]; cindx[0]<last[0]; cindx[0]++) {
          // cout << "cindx = " << cindx[0] << ", " << cindx[1] << ", " << cindx[2] << endl;
          rval += GetCoef(cindx[0],cindx[1],cindx[2]) * _sp.SplineValueAtVoxel(vox,cindx);
	}
      }
    }
  }
  else {
    _dsp[fi-1]->RangeOfSplines(vox,csz,first,last);
    for (cindx[2]=first[2]; cindx[2]<last[2]; cindx[2]++) {
      for (cindx[1]=first[1]; cindx[1]<last[1]; cindx[1]++) {
        for (cindx[0]=first[0]; cindx[0]<last[0]; cindx[0]++) {
          rval += GetCoef(cindx[0],cindx[1],cindx[2]) * _dsp[fi-1]->SplineValueAtVoxel(vox,cindx);
	}
      }
    }
  }
  return(rval);
}

/////////////////////////////////////////////////////////////////////
//
// Member-functions for the ZeroSplineMap class. The class will
// keep track of splines for which the/an image is zero for all
// of their support, thereby avoiding unneccessary calculations.
//
/////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////
//
// Member-functions for the Memen_H_Helper class. The idea behind
// the class is that each column of H contains the same values
// spaced out in a particular patter (though some values might be
// "shifted out" of the volume and may be missing for a given
// column). A Memen_H_Helper object will calculate all unique values
// on construction and then one can use one of the access function
// (operator()(i,j,k) or Peek(i,j,k)) to populate H with this values.
//
/////////////////////////////////////////////////////////////////////

Memen_H_Helper::Memen_H_Helper(const Spline3D<double>&  sp) : _sz(3,0), _cntr(3,0), _data(NULL)
{
  // Fake a "really large" FOV
  std::vector<unsigned int>    isz(3,0);
  for (unsigned int i=0; i<3; i++) isz[i] = 1000*sp.KnotSpacing(i);

  // Pick an index somewhere in the centre
  std::vector<unsigned int>    cindx(3,0);
  for (unsigned int i=0; i<3; i++) cindx[i] = sp.NCoef(i,isz[i]) / 2;

  // Get indices of overlapping splines
  std::vector<unsigned int>    first(3,0), last(3,0);
  if (!sp.RangeOfOverlappingSplines(cindx,isz,first,last)) throw BasisfieldException("Memen_H_Helper::Memen_H_Helper: No overlapping splines");

  _sz[0] = last[0]-first[0]; _sz[1] = last[1]-first[1]; _sz[2] = last[2]-first[2];
  _cntr[0] = cindx[0]-first[0]; _cntr[1] = cindx[1]-first[1]; _cntr[2] = cindx[2]-first[2];
  _data = new double[_sz[0]*_sz[1]*_sz[2]];

  // Get values for all "overlaps" in "all positive" 1/8
  std::vector<unsigned int> cindx2(3,0);
  for (unsigned int ck=first[2]; ck<last[2]; ck++) {
    for (unsigned int cj=first[1]; cj<last[1]; cj++) {
      for (unsigned int ci=first[0]; ci<last[0]; ci++) {
        unsigned int li = (ck-cindx[2]+_cntr[2])*_sz[1]*_sz[0] + (cj-cindx[1]+_cntr[1])*_sz[0] + (ci-cindx[0]+_cntr[0]);
        cindx2[0] = ci; cindx2[1] = cj; cindx2[2] = ck;
        _data[li] = sp.MulByOther(cindx,cindx2,isz,sp);
      }
    }
  }
}

double Memen_H_Helper::operator()(int i, int j, int k) const
{
  if (i+_cntr[0] < 0 || i+_cntr[0] >= _sz[0] ||
      j+_cntr[1] < 0 || j+_cntr[1] >= _sz[1] ||
      k+_cntr[2] < 0 || k+_cntr[2] >= _sz[2]) {
    throw BasisfieldException("Memen_H_Helper::operator(): Index out of range");
  }
  return(Peek(i,j,k));
}

/////////////////////////////////////////////////////////////////////
//
// Member-functions for the helper_vol class. It is a little helper
// class that is used when deconvolving a field to obtain the spline
// coefficients.
//
/////////////////////////////////////////////////////////////////////

helper_vol::helper_vol(const std::vector<unsigned int>& sz)
{
  if (sz.size()!=3) throw BasisfieldException("helper_vol::helper_vol: s must have 3 elements");
  else {
    _sz[0]=sz[0]; _sz[1]=sz[1]; _sz[2]=sz[2];
    _data = new double[_sz[0]*_sz[1]*_sz[2]];
  }
}
helper_vol::helper_vol(const std::vector<unsigned int>& sz,
                       const NEWMAT::ColumnVector&      vec)
{
  if (sz.size()!=3) throw BasisfieldException("helper_vol::helper_vol: s must have 3 elements");
  else {
    _sz[0]=sz[0]; _sz[1]=sz[1]; _sz[2]=sz[2];
    _data = new double[_sz[0]*_sz[1]*_sz[2]];
    double *dp = static_cast<double *>(vec.Store());
    memcpy(_data,dp,_sz[0]*_sz[1]*_sz[2]*sizeof(double));
  }
}
helper_vol::helper_vol(const NEWIMAGE::volume<float>&   vol)
{
  _sz[0]=vol.xsize(); _sz[1]=vol.ysize(); _sz[2]=vol.zsize();
  _data = new double[_sz[0]*_sz[1]*_sz[2]];
  double *trgt = _data;
  for (NEWIMAGE::volume<float>::fast_const_iterator it=vol.fbegin(), it_end=vol.fend(); it!=it_end; it++, trgt++) {
    *trgt = static_cast<double>(*it);
  }
}
helper_vol::helper_vol(const helper_vol& in)
{
  _sz[0]=in._sz[0]; _sz[1]=in._sz[1]; _sz[2]=in._sz[2];
  _data = new double[_sz[0]*_sz[1]*_sz[2]];
  memcpy(_data,in._data,_sz[0]*_sz[1]*_sz[2]*sizeof(double));
}
helper_vol& helper_vol::operator=(const helper_vol& rhs)
{
  if (this == &rhs) return(*this);
  _sz[0]=rhs._sz[0]; _sz[1]=rhs._sz[1]; _sz[2]=rhs._sz[2];
  if (_data) delete [] _data;
  _data = new double[_sz[0]*_sz[1]*_sz[2]];
  memcpy(_data,rhs._data,_sz[0]*_sz[1]*_sz[2]*sizeof(double));
  return(*this);
}
NEWMAT::ReturnMatrix helper_vol::AsNewmatVector() const
{
  NEWMAT::ColumnVector  ovec(_sz[0]*_sz[1]*_sz[2]);
  double *dp = ovec.Store();
  memcpy(dp,_data,_sz[0]*_sz[1]*_sz[2]*sizeof(double));
  ovec.Release();
  return(ovec);
}
void helper_vol::GetColumn(unsigned int i, unsigned int j, unsigned int dim, double *col) const
{
  const double *endptr = end(i,j,dim);
  unsigned int stp = step(dim);
  for (const double *ptr=start(i,j,dim); ptr<endptr; ptr+=stp, col++) *col = *ptr;
}
void helper_vol::SetColumn(unsigned int i, unsigned int j, unsigned int dim, const double *col)
{
  const double *endptr = end(i,j,dim);
  unsigned int stp = step(dim);
  for (double *ptr=start(i,j,dim); ptr<endptr; ptr+=stp, col++) *ptr = *col;
}
double* helper_vol::start(unsigned int i, unsigned int j, unsigned int dim) const
{
  switch(dim) {
  case 0:
    return(&(_data[j*_sz[0]*_sz[1]+i*_sz[0]]));
    break;
  case 1:
    return(&(_data[j*_sz[0]*_sz[1]+i]));
    break;
  case 2:
    return(&(_data[j*_sz[0]+i]));
    break;
  default:
    throw BasisfieldException("helper_vol::start: dim must be 0, 1 or 2");
    break;
  }
}
const double* helper_vol::end(unsigned int i, unsigned int j, unsigned int dim) const
{
  const double *ptr = start(i,j,dim);
  ptr += _sz[dim]*step(dim);
  return(ptr);
}
unsigned int helper_vol::step(unsigned dim) const
{
  switch (dim) {
  case 0:
    return(1);
    break;
  case 1:
    return(_sz[0]);
    break;
  case 2:
    return(_sz[1]*_sz[0]);
  default:
    throw BasisfieldException("helper_vol::start: dim must be 0, 1 or 2");
    break;
  }
}


} // End namespace BASISFIELD

/*

////////////////////////////////////////////////////////////////
//
// Here is some decomissioned code that I don't dare to
// delete just yet.
//
////////////////////////////////////////////////////////////////

  double calculate_bender(const NEWMAT::ColumnVector&        b,
                          const std::vector<unsigned int>&   lksp,
                          const std::vector<unsigned int>&   csz) const;

  double calculate_memen(const NEWMAT::ColumnVector&        b,
                         const std::vector<unsigned int>&   lksp,
                         const std::vector<unsigned int>&   csz) const;

  void calculate_bender_grad(const NEWMAT::ColumnVector&       b,
                             const std::vector<unsigned int>&  lksp,
                             const std::vector<unsigned int>&  csz,
                             NEWMAT::ColumnVector&             grad) const;

  void calculate_memen_grad(const NEWMAT::ColumnVector&       b,
                            const std::vector<unsigned int>&  lksp,
                            const std::vector<unsigned int>&  csz,
                            NEWMAT::ColumnVector&             grad) const;


  NEWMAT::ReturnMatrix memen_HtHb_helper(const Spline3D<double>&            spd,
                                         const std::vector<unsigned int>&   csz,
                                         const NEWMAT::ColumnVector&        b,
                                         HtHbType                           what) const;

double splinefield::MemEnergy() const // Membrane energy of field
{
  const std::shared_ptr<NEWMAT::ColumnVector> lcoef = GetCoef();
  if (!lcoef) {throw BasisfieldException("splinefield::MemEnergy: No coefficients set yet");}

  std::vector<unsigned int>  csz(3,0);
  csz[0] = CoefSz_x(); csz[1] = CoefSz_y(); csz[2] = CoefSz_z();

  return(calculate_memen(*lcoef,_sp.KnotSpacing(),csz));
}

double splinefield::BendEnergy() const // Bending energy of field
{
  const std::shared_ptr<NEWMAT::ColumnVector> lcoef = GetCoef();
  if (!lcoef) {throw BasisfieldException("splinefield::BendEnergy: No coefficients set yet");}

  std::vector<unsigned int>  csz(3,0);
  csz[0] = CoefSz_x(); csz[1] = CoefSz_y(); csz[2] = CoefSz_z();

  return(calculate_bender(*lcoef,_sp.KnotSpacing(),csz));
}

NEWMAT::ReturnMatrix splinefield::MemEnergyGrad() const // Gradient of membrane energy of field
{
  const std::shared_ptr<NEWMAT::ColumnVector> lcoef = GetCoef();
  if (!lcoef) {throw BasisfieldException("splinefield::MemEnergyGrad: No coefficients set yet");}

  std::vector<unsigned int>  csz(3,0);
  csz[0] = CoefSz_x(); csz[1] = CoefSz_y(); csz[2] = CoefSz_z();
  NEWMAT::ColumnVector  grad(CoefSz());

  calculate_memen_grad(*lcoef,_sp.KnotSpacing(),csz,grad);

  grad.Release();
  return(grad);
}

NEWMAT::ReturnMatrix splinefield::BendEnergyGrad() const // Gradient of bending energy of field
{
  const std::shared_ptr<NEWMAT::ColumnVector> lcoef = GetCoef();
  if (!lcoef) {throw BasisfieldException("splinefield::BendEnergyGrad: No coefficients set yet");}

  std::vector<unsigned int>  csz(3,0);
  csz[0] = CoefSz_x(); csz[1] = CoefSz_y(); csz[2] = CoefSz_z();
  NEWMAT::ColumnVector  grad(CoefSz());

  calculate_bender_grad(*lcoef,_sp.KnotSpacing(),csz,grad);

  grad.Release();
  return(grad);
}


// Calculates bending energy

double splinefield::calculate_bender(const NEWMAT::ColumnVector&        b,
                                     const std::vector<unsigned int>&   lksp,
                                     const std::vector<unsigned int>&   csz)
const
{
  double memen = 0.0;
  double vxs[] = {Vxs_x(), Vxs_y(), Vxs_z()};

  // Sum over directions
  for (unsigned int d1=0; d1<3; d1++) {
    for (unsigned int d2=d1; d2<3; d2++) {
      std::vector<unsigned int> deriv(3,0);
      deriv[d1]++;
      deriv[d2]++;
      Spline3D<double>          spd(_sp.Order(),lksp,deriv);  // Spline twice differentiated
      spd /= (vxs[d1]*vxs[d2]);                     // Derivative in mm^{-1}
      NEWMAT::ColumnVector      hb = memen_HtHb_helper(spd,csz,b,Hb);
      if (d1==d2) memen += DotProduct(hb,hb);
      else memen += 2.0*DotProduct(hb,hb);
    }
  }
  return(memen);
}

// Calculates membrane energy

double splinefield::calculate_memen(const NEWMAT::ColumnVector&        b,
                                    const std::vector<unsigned int>&   lksp,
                                    const std::vector<unsigned int>&   csz)
const
{
  double memen = 0.0;
  double vxs[] = {Vxs_x(), Vxs_y(), Vxs_z()};

  // Sum over directions
  for (unsigned int d=0; d<3; d++) {
    std::vector<unsigned int> deriv(3,0);
    deriv[d] = 1;
    Spline3D<double>          spd(_sp.Order(),lksp,deriv);  // Spline differentiated in one direction
    spd /= vxs[d];                                // Derivative in mm^{-1}
    NEWMAT::ColumnVector      hb = memen_HtHb_helper(spd,csz,b,Hb);
    memen += DotProduct(hb,hb);
  }
  return(memen);
}

void splinefield::calculate_bender_grad(const NEWMAT::ColumnVector&       b,
                                        const std::vector<unsigned int>&  lksp,
                                        const std::vector<unsigned int>&  csz,
                                        NEWMAT::ColumnVector&             grad)
const
{
  if (static_cast<unsigned int>(grad.Nrows()) != csz[2]*csz[1]*csz[0]) grad.ReSize(csz[2]*csz[1]*csz[0]);
  grad = 0.0;

  double vxs[] = {Vxs_x(), Vxs_y(), Vxs_z()};

  // Sum over directions
  for (unsigned int d1=0; d1<3; d1++) {
    for (unsigned int d2=d1; d2<3; d2++) {
      std::vector<unsigned int> deriv(3,0);
      deriv[d1]++;
      deriv[d2]++;
      Spline3D<double>          spd(_sp.Order(),lksp,deriv);  // Spline twice differentiated
      spd /= (vxs[d1]*vxs[d2]);                     // Derivative in mm^{-1}
      if (d1==d2) grad += memen_HtHb_helper(spd,csz,b,HtHb);
      else grad += 2.0*memen_HtHb_helper(spd,csz,b,HtHb);
    }
  }
  grad *= 2.0;
}

void splinefield::calculate_memen_grad(const NEWMAT::ColumnVector&       b,
                                       const std::vector<unsigned int>&  lksp,
                                       const std::vector<unsigned int>&  csz,
                                       NEWMAT::ColumnVector&             grad)
const
{
  if (static_cast<unsigned int>(grad.Nrows()) != csz[2]*csz[1]*csz[0]) grad.ReSize(csz[2]*csz[1]*csz[0]);
  grad = 0.0;

  double vxs[] = {Vxs_x(), Vxs_y(), Vxs_z()};

  // Sum over directions
  for (unsigned int d=0; d<3; d++) {
    std::vector<unsigned int> deriv(3,0);
    deriv[d] = 1;
    Spline3D<double>          spd(_sp.Order(),lksp,deriv);  // Spline differentiated in one direction
    spd /= vxs[d];                                // Derivative in mm^{-1}
    grad += memen_HtHb_helper(spd,csz,b,HtHb);
  }
  grad *= 2.0;
}

/////////////////////////////////////////////////////////////////////
//
// This is a helper-routine that calculates H*b or H'*H*b depending
// on the switch "what". H in this case is the matrix such that
// b'*(H'*H)*b is the membrane energy for a field given by the
// coefficients b. It does so without explicitly representing H,
// which is good because it can be very large. It is a helper both
// for calculating the memebrane energy (as (H*b)'*(H*b)) and the
// gradient of the membrane energy (as H'*H*b).
// N.B. that we want to assess the energy for the "entire field",
// i.e. also the parts that extends beyond the FOV of the image.
// That means that (H*b).Nrows() > FieldSz().
//
/////////////////////////////////////////////////////////////////////

NEWMAT::ReturnMatrix splinefield::memen_HtHb_helper(const Spline3D<double>&            spd,
                                                    const std::vector<unsigned int>&   csz,
                                                    const NEWMAT::ColumnVector&        b,
                                                    HtHbType                           what)
const
{
  NEWMAT::ColumnVector  hb(spd.TotalFullFOV(csz));   // H*b
  hb = 0.0;
  double                *hbp = hb.Store();
  const double          *bp = b.Store();

  // Generate indicies of first spline into Hb
  unsigned int *indx = new unsigned int[spd.TotalKernelSize()];
  for (unsigned int k=0, li=0; k<spd.KernelSize(2); k++) {
    unsigned int b1 = k * spd.FullFOV(1,csz[1]) * spd.FullFOV(0,csz[0]);
    for (unsigned int j=0; j<spd.KernelSize(1); j++) {
      unsigned int b2 = j * spd.FullFOV(0,csz[0]);
      for (unsigned int i=0; i<spd.KernelSize(0); i++, li++) {
        indx[li] = b1 + b2 + i;
      }
    }
  }

  // Build H*b as linear combination of the columns of H
  unsigned int is = spd.KnotSpacing(0);
  unsigned int js = spd.KnotSpacing(1) * spd.FullFOV(0,csz[0]);
  unsigned int ks = spd.KnotSpacing(2) * spd.FullFOV(1,csz[1]) * spd.FullFOV(0,csz[0]);
  for (unsigned int ck=0, lci=0; ck<csz[2]; ck++) {
    for (unsigned int cj=0; cj<csz[1]; cj++) {
      for (unsigned int ci=0; ci<csz[0]; ci++, lci++) {
        unsigned int offset = ck*ks + cj*js + ci*is;
        for (unsigned int i=0; i<spd.TotalKernelSize(); i++) {
          hbp[indx[i]+offset] += bp[lci] * spd[i];
	}
      }
    }
  }
  if (what == Hb) { // Return if that is all we shall do
    delete [] indx;
    hb.Release();
    return(hb);
  }

  // Multiply by H'
  NEWMAT::ColumnVector   hthb(csz[2]*csz[1]*csz[0]);
  hthb = 0.0;
  double                 *hthbp = hthb.Store();
  for (unsigned int ck=0, lci=0; ck<csz[2]; ck++) {
    for (unsigned int cj=0; cj<csz[1]; cj++) {
      for (unsigned int ci=0; ci<csz[0]; ci++, lci++) {
        unsigned int offset = ck*ks + cj*js + ci*is;
        for (unsigned int i=0; i<spd.TotalKernelSize(); i++) {
          hthbp[lci] += spd[i] * hbp[indx[i]+offset];
	}
      }
    }
  }
  delete [] indx;
  hthb.Release();
  return(hthb);
}

*/
