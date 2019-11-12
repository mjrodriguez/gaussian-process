#include "../include/kernel.h"
#include "../include/Matrix.h"
#include<cmath>

template<class T> void GenericKernel<T>::ConstructKernel(const T* params, const Vector<T>& x, const Vector<T>& y, Matrix<T>& K){

    for (int i=0; i<K.GetNumberOfRows(); i++){
        for(int j=0; j<K.GetNumberOfCols(); j++){
            K(i,j) = Kernel(params, x.Read(i), y.Read(j));
        }
    }

}

template<class T> T SquaredExponential<T>::Kernel(const T* params, const T x, const T y){
    return params[0]*std::exp(-0.5*std::pow(x-y,2)/params[1]);
}

template class GenericKernel<double>;
template class GenericKernel<float>;
template class SquaredExponential<double>;
template class SquaredExponential<float>;