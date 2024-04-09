//*****************************************************************
// Iterative template routine -- BiCGSTAB
//
// BiCGSTAB solves the unsymmetric linear system Ax = b 
// using the Preconditioned BiConjugate Gradient Stabilized method
//
// BiCGSTAB follows the algorithm described on p. 27 of the 
// SIAM Templates book.
//
// The return value indicates convergence within max_iter (input)
// iterations (0), or no convergence within max_iter iterations (1).
//
// Upon successful return, output arguments have the following values:
//  
//        x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//  
//*****************************************************************
//
// Slightly modified version of IML++ template. See ReadMe file.
//
// Jesper Andersson
//
/*    Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */
#ifndef bicgstab_h
#define bicgstab_h

namespace MISCMATHS {

template < class Matrix, class Vector, class Preconditioner, class Real >
int 
BiCGSTAB(const Matrix &A, Vector &x, const Vector &b,
         const Preconditioner &M, int &max_iter, Real &tol)
{
  Real resid;
  Vector rho_1(1), rho_2(1), alpha(1), beta(1), omega(1);
  Vector p, phat, s, shat, t, v;

  Real normb = b.NormFrobenius();
  Vector r = b - A * x;
  Vector rtilde = r;

  if (normb == 0.0)
    normb = 1;
  
  if ((resid = r.NormFrobenius() / normb) <= tol) {
    tol = resid;
    max_iter = 0;
    return 0;
  }

  for (int i = 1; i <= max_iter; i++) {
    rho_1(1) = DotProduct(rtilde, r);
    if (rho_1(1) == 0) {
      tol = r.NormFrobenius() / normb;
      return 2;
    }
    if (i == 1)
      p = r;
    else {
      beta(1) = (rho_1(1)/rho_2(1)) * (alpha(1)/omega(1));
      p = r + beta(1) * (p - omega(1) * v);
    }
    phat = M.solve(p);
    v = A * phat;
    alpha(1) = rho_1(1) / DotProduct(rtilde, v);
    s = r - alpha(1) * v;
    if ((resid = s.NormFrobenius()/normb) < tol) {
      x += alpha(1) * phat;
      tol = resid;
      return 0;
    }
    shat = M.solve(s);
    t = A * shat;
    omega = DotProduct(t,s) / DotProduct(t,t);
    x += alpha(1) * phat + omega(1) * shat;
    r = s - omega(1) * t;

    rho_2(1) = rho_1(1);
    if ((resid = r.NormFrobenius() / normb) < tol) {
      tol = resid;
      max_iter = i;
      return 0;
    }
    if (omega(1) == 0) {
      tol = r.NormFrobenius() / normb;
      return 3;
    }
  }

  tol = resid;
  return 1;
}

} // End namespace MISCMATHS

#endif // End #ifndef bicgstab_h
