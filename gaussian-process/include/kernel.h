#ifndef KERNEL_H
#define KERNEL_H

#include "Vector.h"
#include "Matrix.h"

template<class T> class GenericKernel{
public:
    virtual T Kernel(const T* params, const T x, const T y) = 0;
    void ConstructKernel(const T* params, const Vector<T>& x, const Vector<T>& y, Matrix<T>& K);
};


template<class T> class SquaredExponential : public GenericKernel<T>{
public:
    T Kernel(const T* params, const T x, const T y);
};


#endif