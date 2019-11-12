#include <cmath>
#include <iostream>
#include <cassert>
#include "../include/Vector.h"
#include "../include/Matrix.h"

template<class T> void Matrix<T>::memalloc(){
    mData = new T* [mRows];
    for (int i=0; i<mRows; i++){
        mData[i] = new T [mCols];
    }
}

template<class T> Matrix<T>::Matrix(const Matrix<T>& otherMatrix){
    mRows = otherMatrix.mRows;
    mCols = otherMatrix.mCols;

    memalloc();

    for (int i = 0; i < mRows; i++){
        for (int j = 0; j < mCols; j++){
            mData[i][j] = otherMatrix.mData[i][j];
        }
    }
}

template<class T> Matrix<T>::Matrix(Vector<T>& vec){
    mRows = vec.GetSize(); mCols = vec.GetSize();
    memalloc();
    for (int i = 0; i < mRows; i++){
        for (int j = 0; j < mCols; j++){
            mData[i][j] = T(0);
            if ( i == j) mData[i][i] = vec[i];
        }
    }
}

template<class T> Matrix<T>::Matrix(int numRows, int numCols){
    assert(numRows > 0); assert(numCols > 0);

    mRows = numRows; mCols = numCols;
    memalloc();
    for (int i = 0; i < mRows; i++){
        for (int j = 0; j < mCols; j++){
            mData[i][j] = T(0);
        }
    }
}

template<class T> Matrix<T>::Matrix(int numRows, int numCols, T val){
    assert(numRows > 0); assert(numCols > 0);

    mRows = numRows; mCols = numCols;
    memalloc();
    for (int i = 0; i < mRows; i++){
        for (int j = 0; j < mCols; j++){
            mData[i][j] = val;
        }
    }
}


template<class T> Matrix<T>::~Matrix(){
    for (int i=0; i<mRows; i++){
        delete[] mData[i];
    }
    delete mData;
}

template<class T> int Matrix<T>::GetNumberOfRows() const {
    return mRows;
}
template<class T> int Matrix<T>::GetNumberOfCols() const {
    return mCols;
}


template<class T> Matrix<T> Matrix<T>::Transpose() const {
    Matrix<T> mat(mCols, mRows);
    for (int i=0; i<mRows; i++){
        for (int j=0; j<mCols; j++){
            mat(j,i) = mData[i][j];
        }
    }
    return mat;
}


template<class T> T& Matrix<T>::operator()(int i, int j)const{
    assert(i > -1 && i < mRows);
    assert(j > -1 && j < mCols);
    return mData[i][j];
}

template<class T> Matrix<T>& Matrix<T>::operator=(const Matrix<T>& otherMatrix){
    assert(mRows == otherMatrix.mRows);
    assert(mCols == otherMatrix.mCols);
    for (int i=0; i<mRows; i++){
        for (int j=0; j<mCols; j++){
            mData[i][j] = otherMatrix.mData[i][j];
        }
    }
    return *this;
}

template<class T> Matrix<T> Matrix<T>::operator+() const{
    Matrix<T> mat(mRows,mCols);
    for (int i=0; i<mRows; i++){
        for (int j=0; j<mCols; j++){
            mat(i,j) = mData[i][j];
        }
    }
    return mat;
}

template<class T> Matrix<T> Matrix<T>::operator-() const{
    Matrix<T> mat(mRows,mCols);
    for (int i=0; i<mRows; i++){
        for (int j=0; j<mCols; j++){
            mat(i,j) = -mData[i][j];
        }
    }
    return mat;
}

template<class T> Matrix<T> Matrix<T>::operator+(const Matrix<T>& m1) const{
    assert(mRows == m1.mRows && mCols == m1.mCols);

    Matrix<T> mat(mRows,mCols);
    for (int i=0; i<mRows; i++){
        for (int j=0; j<mCols; j++){
            mat(i,j) = mData[i][j] + m1.mData[i][j];
        }
    }
    return mat;
}

template<class T> Matrix<T> Matrix<T>::operator-(const Matrix<T>& m1) const{
    assert(mRows == m1.mRows && mCols == m1.mCols);

    Matrix<T> mat(mRows,mCols);
    for (int i=0; i<mRows; i++){
        for (int j=0; j<mCols; j++){
            mat(i,j) = mData[i][j] - m1.mData[i][j];
        }
    }
    return mat;
}

template<class T> Matrix<T> Matrix<T>::Scale(const T a) const{
    Matrix<T> mat(mRows,mCols);
    for (int i=0; i<mRows; i++){
        for (int j=0; j<mCols; j++){
            mat(i,j) = a*mData[i][j];
        }
    }
    return mat;
}

// Matrix-Matrix Multiplication
template<class T> Matrix<T> Matrix<T>::Multiply(const Matrix<T>& otherMat) const{
    assert(mCols == otherMat.GetNumberOfRows());
    int newNumRows = mRows;
    int newNumCols = otherMat.GetNumberOfCols();
    assert(newNumRows > 0 && newNumCols > 0);
    Matrix<T> resultMat(newNumRows,newNumCols);

    for(int i=0; i<mRows; i++){
        for(int j=0; j<otherMat.GetNumberOfRows(); j++){
            T sum = T(0);
            for(int k=0; k<mCols; k++){
                sum += mData[i][k]*otherMat(k,j);
            }
            resultMat(i,j) = sum;
        }
    }

    return resultMat;
}

// Matrix-Vector Multiplication
template<class T> Vector<T> Matrix<T>::Multiply(const Vector<T>& vec) const{
    assert(mCols == vec.GetSize());
    Vector<T> resultVec(mRows);
    for (int i=0; i<mRows; i++){
        for (int j=0; j<vec.GetSize(); j++){
            resultVec[i] += mData[i][j]*vec.Read(j);
        }
    }
    
    return resultVec;
}


template<class T> T** Matrix<T>::RawPointer(){
    return mData;
}

template<class T> T Matrix<T>::FrobeniusNorm() const{
    T sum = T(0);
    for (int i=0; i<mRows; i++){
        for (int j=0; j<mCols; j++){
            sum += abs(mData[i][j])*abs(mData[i][j]);
        }
    }
    return sqrt(sum);
}


template<class T> void Matrix<T>::Print() const {
    Print("");
}

template<class T> void Matrix<T>::Print(std::string name) const {
    std::cout << mRows << "x" << mCols << " Matrix " << name << std::endl;
    for (int i=0; i<mRows; i++){
        std::cout << "[ ";
        for (int j=0; j<mCols; j++){
            std::cout << mData[i][j] << " ";
        }
        std::cout << " ]" << std::endl;
    }
}





// template<class T> std::ostream& operator<<(std::ostream& output, Matrix<T>& matrix){
//     for (int i = 0; i < matrix.mRows; ++i) {
//         output << "[";
//         for (int j = 0; j < matrix.mCols; ++j) {
//             output << matrix[i][j];
//             if (j < matrix.mCols-1){
//                 output << " ";
//             }
//         }
//         output << " ]\n";
//     }
//     return output;
// }


template class Matrix<double>;
template class Matrix<float>;
template class Matrix<int>;
