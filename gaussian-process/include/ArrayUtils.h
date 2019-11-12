#pragma once

#include "Array.h"
#include <string>
#include <fstream>
#include <iostream>

// y = alpha*A*x + beta*y
void matvec(DArray &y, const DArray &A, const DArray &x, bool trans=false,
            double alpha=1.0, double beta=0.0);
// C = alpha*A*B + beta*C
void matmat(DArray &C, const DArray &A, const DArray &B, bool TA=false,
            bool TB=false, double alpha=1.0, double beta=0.0);
// compute permuted LU
void LU(DArray &A, IArray &P);
// compute Cholesky factorization
void chol(DArray &A);
// using LU of A
void solve(DArray &B, const DArray &LUA, const IArray &P, bool trans=false);
// using chol of A
void cholsolve(DArray &B, const DArray &cholA);
// invert using LU of A
void inv(DArray &A, const IArray &P);
// y = x
template <typename T>
void copy(TArray<T> &y, const TArray<T> &x) {
  #ifdef kBoundsCheck
    assert(y.nel == x.nel);
  #endif
  std::copy(x.begin(), x.end(), y.begin());
}
void copy(DArray &y, const DArray &x);
// B = A^T
template <typename T>
void transp(TArray<T> &B, TArray<T> &A) {
  #ifdef kBoundsCheck
    assert(B.ndim == 2);
    assert(A.ndim == 2);
    assert(A.size(0) == B.size(1) && A.size(1) == B.size(0));
  #endif

  for (int j=0; j<B.size(1); ++j) {
    for (int i=0; i<B.size(0); ++i) {
      B(i,j) = A(j,i);
    }
  }
}
// A = A^T
template <typename T>
void transp(TArray<T> &A) {
  #ifdef kBoundsCheck
    assert(A.size(0) == A.size(1));
  #endif

  for (int j=0; j<A.size(1); ++j) {
    for (int i=0; i<j; ++i) {
      std::swap(A(i,j), A(j,i));
    }
  }
}

// y = ax + y
void axpy(DArray &y, const DArray &x, double a=1.0);
// y = ax + by
void axpby(DArray &y, double b, const DArray &x, double a);
//
double dot(const DArray &x, const DArray &y);
// 2-norm
double norm(const DArray &x);

double maxAbs(const DArray &x);
int argMax(const DArray &x);

template <typename T>
bool anynan(const TArray<T> &x) {
  for (auto a : x)
    if (std::isnan(a)) return true;
  return false;
}

void scale(DArray &x, double a);
void randu(DArray &x, double vmin=0.0, double vmax=1.0);
void kron(DArray &A, const DArray &B, const DArray &C, double alpha=1.0,
          double beta=0.0);
void kronSolve(DArray &x, const DArray &ALU, const IArray &AP,
               const DArray &BLU, const IArray &BP, DArray &tmp);
void kronCholsolve(DArray &x, const DArray &A, const DArray &B, DArray &tmp);
void kronSolve(DArray &x, const DArray &ALU, const IArray &AP,
               const DArray &BLU, const IArray &BP, const DArray &CLU,
               const IArray &CP, DArray &tmp1, DArray &tmp2, double alpha=1.0);
void kronMult(DArray &x, const DArray &A, const DArray &B, DArray &tmp,
              bool trans=false);
void kronMult(DArray &y, DArray &x, const DArray &A, const DArray &B,
              const DArray &C, DArray &tmp1, DArray &tmp2, bool trans=false,
              double alpha=1.0, double beta=0.0);

void solveSylvester(DArray &C, const DArray &A, const DArray &B,
                    bool trana=false, bool tranb=true, double scale=1.0);
void schurDecomp(DArray &A, DArray &QZ);
template <typename T>
void writeRaw(TArray<T> &x, const std::string &filename) {
  std::ofstream f(filename);
  f.write((char *)x.data, sizeof(T)*x.nel);
  f.close();
}

// Based on this StackExchange answer: http://stackoverflow.com/a/26221725
template<typename ... Args>
std::string stringFormat(const std::string& format, Args ... args) {
  std::string result;
  // Extra space for the null character
  size_t size = std::snprintf(nullptr, 0, format.c_str(), args ...) + 1;
  result.resize(size);
  std::snprintf(&result[0], size, format.c_str(), args ...);
  result.resize(size-1);
  return result;
}

std::string formatElement(int x);
std::string formatElement(double x);

template <typename T>
void print(const TArray<T> &x, std::string name) {
  switch (x.ndim) {
    case 1:
      std::cout << stringFormat("Size %d vector ", x.dims[0])
                << name << std::endl;
      for (int i = 0; i < x.dims[0]; ++i) {
        //std::cout << "\t" << x(i) << std::endl;
        std::cout << "    " << formatElement(x(i)) << std::endl;
      }
      break;
    case 2:
      std::cout << stringFormat("Size %dx%d matrix ", x.dims[0], x.dims[1])
                << name << std::endl;
      for (int i = 0; i < x.dims[0]; ++i) {
        std::cout << "    ";
        for (int j = 0; j < x.dims[1]; ++j) {
          //std::cout << stringFormat("% 16.10e  ", x(i,j));
          std::cout << formatElement(x(i,j));
        }
        std::cout << std::endl;
      }
      break;
    case 3:
      std::cout << stringFormat("Size %dx%dx%d tensor ",
                                x.dims[0], x.dims[1], x.dims[2])
                << name << std::endl;
      for (int k = 0; k < x.dims[2]; ++k) {
        TArray<T> slice(&x(0,0,k), x.dims[0], x.dims[1]);
        print(slice, stringFormat("(:,:,%d)", k));
      }
      break;
    default:
      std::cout << "Cannot output arrays of dimension more than 2"
                << std::endl;
      break;
  }
}

template <typename T>
void print(const TArray<T> &x) {
  print(x, "");
}

void print(const DArray &x);
void print(const IArray &x);
