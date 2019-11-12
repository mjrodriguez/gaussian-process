#include "ArrayUtils.h"
#include "LapackWrapper.h"
#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>

std::mt19937& re();
std::mt19937& re() {
  static std::mt19937 r;
  return r;
}

void matvec(DArray &y, const DArray &A, const DArray &x, bool trans,
            double alpha, double beta) {
  CBLAS_TRANSPOSE T = trans ? CblasTrans : CblasNoTrans;
  #ifdef kBoundsCheck
    if (trans) {
      assert(y.size(0) == A.size(1));
      assert(x.size(0) == A.size(0));
    } else {
      assert(y.size(0) == A.size(0));
      assert(x.size(0) == A.size(1));
    }
  #endif
  cblas_dgemv(CblasColMajor, T, A.size(0), A.size(1), alpha, A,
              A.size(0), x, x.strides[0], beta, y, y.strides[0]);
}

// C = alpha*A*B + beta*C
void matmat(DArray &C, const DArray &A, const DArray &B, bool TA, bool TB,
            double alpha, double beta) {
  CBLAS_TRANSPOSE BlTA = TA ? CblasTrans : CblasNoTrans;
  CBLAS_TRANSPOSE BlTB = TB ? CblasTrans : CblasNoTrans;
  int m = TA ? A.size(1) : A.size(0);
  int n = TB ? B.size(0) : B.nel/B.size(0);
  int k = TA ? A.size(0) : A.size(1);
  int lda = A.strides[1];
  int ldb = B.strides[1];
  int ldc = C.strides[1];
  cblas_dgemm(CblasColMajor, BlTA, BlTB, m, n, k, alpha, A, lda,
              B, ldb, beta, C, ldc);
}

void LU(DArray &A, IArray &P) {
  if (dgetrfCpp(A.size(0), A.size(1), A, A.size(0), P) != 0) {
    throw std::runtime_error("LAPACK DGETRF failed to factorize");
  }
}

void chol(DArray &A) {
  if (dpotrfCpp('U', A.size(0), A, A.size(0)) != 0) {
    throw std::runtime_error("LAPACK DPOTRF failed to factorize");
  }
}

void solve(DArray &B, const DArray &LUA, const IArray &P, bool trans) {
  char t = trans ? 'T' : 'N';
  int nrhs = (B.ndim > 1) ? B.nel/B.size(0) : 1;
  if (dgetrsCpp(t, LUA.size(0), nrhs, LUA, LUA.size(0), P, B, B.size(0)) != 0) {
    throw std::runtime_error("LAPACK DGETRS failed to solve system");
  }
}

void cholsolve(DArray &B, const DArray &cholA) {
  int nrhs = (B.ndim > 1) ? B.nel/B.size(0) : 1;
  int N = cholA.size(0);
  if (dpotrsCpp('U', N, nrhs, cholA, N, B, B.size(0)) != 0) {
    throw std::runtime_error("LAPACK DPOTRS failed to solve system");
  }
}

void inv(DArray &A, const IArray &P) {
  dgetriCpp(A.size(0), A, A.strides[1], P);
}

double dot(const DArray &x, const DArray &y) {
  #ifdef kBoundsCheck
    assert(x.nel == y.nel);
  #endif
  return cblas_ddot(x.nel, x, x.strides[0], y, y.strides[0]);
}

void copy(DArray &y, const DArray &x) {
  #ifdef kBoundsCheck
    assert(y.nel == x.nel);
  #endif
  if (x.strides[0] == 1 && y.strides[0] == 1) {
    std::copy(x.begin(), x.end(), y.begin());
  } else {
    cblas_dcopy(x.nel, x, x.strides[0], y, y.strides[0]);
  }
}

void axpy(DArray &y, const DArray &x, double a) {
  #ifdef kBoundsCheck
    assert(y.nel == x.nel);
  #endif
  cblas_daxpy(y.nel, a, x, x.strides[0], y, y.strides[0]);
}

// y = ax + by
void axpby(DArray &y, double b, const DArray &x, double a) {
  #ifdef kBoundsCheck
    assert(y.nel == x.nel);
  #endif

  if (b != 1.0) {
    cblas_dscal(y.nel, b, y, y.strides[0]);
  }

  if (a != 0.0) {
    cblas_daxpy(y.nel, a, x, x.strides[0], y, y.strides[0]);
  }
}

double norm(const DArray &x) {
  return cblas_dnrm2(x.nel, x, x.strides[0]);
}

double maxAbs(const DArray &x) {
  if (anynan(x)) return std::numeric_limits<double>::quiet_NaN();
  int idx = cblas_idamax(x.nel, x, x.strides[0]);
  return fabs(x[idx]);
}

int argMax(const DArray &x) {
  return cblas_idamax(x.nel, x, x.strides[0]);
}

void scale(DArray &x, double a) {
  cblas_dscal(x.nel, a, x, x.strides[0]);
}

void randu(DArray &x, double vmin, double vmax) {
  std::uniform_real_distribution<double> unif(vmin,vmax);

  for (int i = 0; i < x.nel; ++i) {
    x[i] = unif(re());
  }
}

// A = B \otimes C, A preallocated
void kron(DArray &A, const DArray &B, const DArray &C, double alpha,
          double beta) {
  #ifdef kBoundsCheck
    assert(A.size(0) == B.size(0)*C.size(0));
    assert(A.size(1) == B.size(1)*C.size(1));
  #endif

  int n = C.size(0);
  int m = C.size(1);

  for (int i=0; i<B.size(0); ++i) {
    for (int j=0; j<B.size(1); ++j) {
      for (int p=0; p<C.size(0); ++p) {
        for (int q=0; q<C.size(1); ++q) {
          A(i*n + p, j*m + q) = alpha*B(i,j)*C(p,q) + beta*A(i*n + p, j*m + q);
        }
      }
    }
  }
}

// Multiply by inverse (A \otimes B)^{-1} x
// by computing B^{-1}X A^{-T}
void kronSolve(DArray &x, const DArray &ALU, const IArray &AP,
               const DArray &BLU, const IArray &BP, DArray &tmp) {
  solve(x, BLU, BP);
  transp(tmp, x);
  solve(tmp, ALU, AP);
  transp(x, tmp);
}

void kronCholsolve(DArray &x, const DArray &A, const DArray &B, DArray &tmp) {
  cholsolve(x, B);
  transp(tmp, x);
  cholsolve(tmp, A);
  transp(x, tmp);
}

// Multiply (A \otimes B)X by computing BXA^T (or (A^T \otimes B^T)X)
void kronMult(DArray &x, const DArray &A, const DArray &B, DArray &tmp,
              bool trans) {
  matmat(tmp, B, x, trans, false);
  matmat(x, tmp, A, false, !trans);
}

// compute y = alpha (A \otimes B \otimes C) x + beta * y
// A, B, C are all (m x n) matrices (or transposed if trans==true)
// tmp1.size() = {m,n,n}
// tmp2.size() = {m*m,n}

// compute y = alpha(A \otimes B \ptimes C)x + beta*y
// sizes: x (n3,n2,n1), y(m3,m2,m1)
//        A (m1,n1), B (m2,n2), C (m3,n3)
//        tmp1 (m3,n2,n1), tmp2(m3*m2,n1)
void kronMult(DArray &y, DArray &x, const DArray &A, const DArray &B,
              const DArray &C, DArray &tmp1, DArray &tmp2, bool trans,
              double alpha, double beta) {

  int m1 = A.size(0);
  int n1 = A.size(1);
  int m2 = B.size(0);
  int n2 = B.size(1);
  int m3 = C.size(0);
  int n3 = C.size(1);

  int lda = m1;
  int ldb = m2;
  int ldc = m3;

  #ifdef kBoundsCheck
    assert(y.size(0) == m3 && y.size(1) == m2 && y.size(2) == m1);
    assert(x.size(0) == n3 && x.size(1) == n2 && x.size(2) == n1);
    assert(tmp1.size(0) == m3 && tmp1.size(1) == n2 && tmp1.size(2) == n1);
    assert(tmp2.size(0) == m3*m2 && tmp2.size(1) == n1);
  #endif

  CBLAS_TRANSPOSE  T =   trans  ? CblasTrans : CblasNoTrans;
  CBLAS_TRANSPOSE NT = (!trans) ? CblasTrans : CblasNoTrans;

  cblas_dgemm(CblasColMajor, T, CblasNoTrans, m3, n1*n2, n3, 1.0, C, ldc,
              x, n3, 0.0, tmp1, m3);
  for (int k=0; k<n1; ++k) {
    cblas_dgemm(CblasColMajor, CblasNoTrans, NT, m3, m2, n2, 1.0, &tmp1(0,0,k),
                m3, B, ldb, 0.0, &tmp2(0,k), m3);
  }
  cblas_dgemm(CblasColMajor, CblasNoTrans, NT, m2*m3, m1, n1, alpha, tmp2,
              m2*m3, A, lda, beta, y, m2*m3);
}

void kronSolve(DArray &x, const DArray &ALU, const IArray &AP,
               const DArray &BLU, const IArray &BP, const DArray &CLU,
               const IArray &CP, DArray &tmp1, DArray &tmp2,
               double alpha) {
  int n3 = ALU.size(0);
  int n2 = BLU.size(0);
  int n1 = CLU.size(0);
  dgetrsCpp('N', n1, n2*n3, CLU, n1, CP, x, n1);
  for (int k=0; k<n3; ++k) {
    for (int j=0; j<n1; ++j) {
      for (int i=0; i<n2; ++i) {
        tmp1(i,j,k) = x(j,i,k);
      }
    }
    dgetrsCpp('N', n2, n1, BLU, n2, BP, &tmp1(0,0,k), n2);
  }
  for (int k=0; k<n2; ++k)
    for (int j=0; j<n1; ++j)
      for (int i=0; i<n3; ++i)
        tmp2(i,j,k) = tmp1(k,j,i);
  dgetrsCpp('N', n3, n1*n2, ALU, n3, AP, tmp2, n3);
  for (int k=0; k<n3; ++k)
    for (int j=0; j<n2; ++j)
      for (int i=0; i<n1; ++i)
        x(i,j,k) = alpha*tmp2(k,i,j);
}

void solveSylvester(DArray &C, const DArray &A, const DArray &B, bool trana,
                    bool tranb,double scale) {
  int m = A.size(0);
  int n = B.size(0);
  char TA = trana ? 'T' : 'N';
  char TB = tranb ? 'T' : 'N';
  int info = dtrsylCpp(TA, TB, 1, m, n, A, m, B, n, C, m, scale);
  assert(info == 0);
}

void schurDecomp(DArray &A, DArray &QZ) {
  int n = A.size(0);
  DArray tau(n-1), wr(n), wi(n);
  int ilo=1, ihi=n;
  // it seems unnecessary to balance the matrix for better conditioning
  // transform to Hessenberg form
  dgehrdCpp(n, ilo, ihi, A, n, tau);
  copy(QZ, A);
  // generate the orthogonal similarity matrix
  dorghrCpp(n, ilo, ihi, QZ, n, tau);
  // transform to quasitriangular (i.e. complete Schur decomposition)
  dhseqrCpp('S', 'V', n, 1, n, A, n, wr, wi, QZ, n);
}

std::string formatElement(int x) {
  return stringFormat("%-6d  ", x);
}

std::string formatElement(double x) {
  return stringFormat("% 16.10e  ", x);
}

void print(const DArray &x) {
  print(x, "");
}

void print(const IArray &x) {
  print(x, "");
}
