#ifndef MATRIX_H
#define MATRIX_H

#include "Vector.h"
#include<string>

template<class T> class Matrix{
private:
    int mRows, mCols;
    T** mData;

public:
    // Constructors
    Matrix(const Matrix<T>& otherMatrix);
    Matrix(Vector<T>& vec); // Construct a square diag matrix from vector<T>
    Matrix(int numRows, int numCols);
    Matrix(int numRows, int numCols, const T val);
    ~Matrix();

    int GetNumberOfRows() const;
    int GetNumberOfCols() const;
    Matrix<T> Transpose() const;
    T& operator()(int i, int j) const;
    Matrix<T>& operator=(const Matrix<T>& otherMatrix);
    Matrix<T> operator+() const;
    Matrix<T> operator-() const;
    Matrix<T> operator+(const Matrix<T>& m1) const;
    Matrix<T> operator-(const Matrix<T>& m1) const;
    Matrix<T> Scale(const T a) const;
    Matrix<T> Multiply(const Matrix<T>& otherMat) const;
    Vector<T> Multiply(const Vector<T>& vec) const;
    // Matrix<T> operator*(T a) const;

    T** RawPointer();

    T FrobeniusNorm() const;
    void memalloc();
    void Print() const;
    void Print(std::string str) const;
    // friend Vector<T> operator*(const Matrix<T>& m, const Vector<T>& v);
    // friend std::ostream& operator<<(std::ostream& output, Matrix<T>& matrix);

};
// Vector<T> operator*(const Matrix<T>& m, const Vector<T>& v);

// Vector<double> operator*(const Matrix<double>& m, const Vector<double>& v);
// Vector<float> operator*(const Matrix<float>& m, const Vector<float>& v);
// Vector<int> operator*(const Matrix<int>& m, const Vector<int>& v);

#endif
