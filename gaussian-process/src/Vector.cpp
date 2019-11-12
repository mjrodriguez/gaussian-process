#include <cmath>
#include <iostream>
#include <cassert>
#include "../include/Vector.h"


// Initializes a Vector of length size and values from another Vector
template<class T> Vector<T>::Vector(const Vector<T>& otherVector){
    mSize = otherVector.GetSize();
    mData = new T [mSize];
    for (int i = 0; i < mSize; i++){
        mData[i] = otherVector.mData[i];
    }
}

// Initializes a Vector of length size with zeros
template<class T> Vector<T>::Vector(int size){
    assert(size > 0);
    mSize = size;
    mData = new T [mSize];
    for (int i = 0; i < mSize; i++){
        mData[i] = 0.0;
    }
}

// Initializes a Vector of length size with val
template<class T> Vector<T>::Vector(int size, T val){
    assert(size > 0);
    mSize = size;
    mData = new T [mSize];
    for (int i = 0; i < mSize; i++){
        mData[i] = val;
    }
}

// Desctructor fo Vector
template<class T> Vector<T>::~Vector(){
    delete[] mData;
}

// Returns the size of Vector
template<class T> int Vector<T>::GetSize() const {
    return mSize;
}

template<class T> T Vector<T>::Max() const {
    T val;
    val = mData[0];
    for (int i=0; i< mSize; i++){
        if (mData[i] > val) val = mData[i];
    }
    return val;
}

template<class T> T Vector<T>::Min() const {
    T val = 0;
    for (int i=0; i< mSize; i++){
        if (mData[i] < val) val = mData[i];
    }
    return val;
}


template<class T> T& Vector<T>::operator[](int i){
    assert(i > -1); assert(i < mSize);
    return mData[i];
}

// Equating two vectors together
template<class T> Vector<T>& Vector<T>::operator=(const Vector<T>& otherVector){
    assert(mSize == otherVector.mSize);
    for (int i=0; i<mSize; i++){
        mData[i] = otherVector.mData[i];
    }
    return *this;
}

template<class T> Vector<T> Vector<T>::operator+() const{
    Vector<T> v(mSize);
    for (int i=0; i<mSize; i++){
        v[i] = mData[i];
    }
    return v;
}

// Multiplying -1*V = -V
template<class T> Vector<T> Vector<T>::operator-() const{
    Vector<T> v(mSize);
    for (int i=0; i<mSize; i++){
        v[i] = -mData[i];
    }
    return v;
}

// Adding to Vectors together
template<class T> Vector<T> Vector<T>::operator+(const Vector<T>& v1) const{
    assert(mSize == v1.mSize);
    Vector<T> v(mSize);
    for (int i=0; i< mSize; i++){
        v[i] = mData[i] + v1.mData[i];
    }
    return v;
}

// Subtracting to Vectors together
template<class T> Vector<T> Vector<T>::operator-(const Vector<T>& v1) const{
    assert(mSize == v1.mSize);
    Vector<T> v(mSize);
    for (int i=0; i< mSize; i++){
        v[i] = mData[i] - v1.mData[i];
    }
    return v;
}

// Scalar multiplication
template<class T> Vector<T> Vector<T>::Scale(const T a) const{
    Vector<T> v(mSize);
    for (int i=0; i< mSize; i++){
        v[i] = a*mData[i];
    }
    return v;
}

template<class T> T* Vector<T>::RawPointer(){
    return mData;
}

// Dot Product
template<class T> T Vector<T>::Dot(const Vector<T> otherVec) const{
    assert(mSize == otherVec.GetSize());
    T dotResult = 0;

    for (int i=0; i<mSize; i++){
        dotResult += mData[i]*otherVec.mData[i];
    }

    return dotResult;
}

// Calculates norm for some int p. default p = 2;
template<class T> T Vector<T>:: CalculateNorm(int p) const{
    assert(p > 0);
    T norm_val, sum = 0.0;
    for (int i=0; i<mSize; i++){
        sum += pow(fabs(mData[i]),p);
    }
    norm_val = pow(sum,T(1.0)/ (T(p)));
    return norm_val;
}

template<class T> void Vector<T>::Print(std::string name) const{
    std::cout << mSize <<  "x1 Vector " << name << std::endl;
    for(int i = 0; i < mSize; i++){
        std::cout << "[ " << mData[i] << " ]" << std::endl;
    }
}

template<class T> void Vector<T>::Print() const {
    Print("");
}

template<class T> T Vector<T>::Read(int i) const {
    assert(i > -1 && i < mSize);
    return mData[i];
}

template class Vector<double>;
template class Vector<float>;
template class Vector<int>;
