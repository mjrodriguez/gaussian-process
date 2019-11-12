#ifndef VECTOR_H
#define VECTOR_H

#include<string>

template<class T> class Vector{

private:
    int mSize;
    T* mData;

public:
    // Constructors
    Vector(const Vector<T>& otherVector);
    Vector(int size);
    Vector(int size, T val);

    // Destructors
    ~Vector();

    int GetSize() const;
    T Max() const;
    T Min() const;
    T Read(int i) const;
    T& operator[](int i); // zero based indexing


    Vector<T>& operator=(const Vector<T>& otherVector);
    Vector<T> operator+() const; // unary +
    Vector<T> operator-() const; // unary -
    Vector<T> operator+(const Vector<T>& v1) const; // binary +
    Vector<T> operator-(const Vector<T>& v1) const; // binary -
    // Vector<T> operator*(T a) const;
    Vector<T>  Scale(const T a) const;
    T* RawPointer();
    // Dot product
    T Dot(const Vector<T> otherVec) const;
    // Compute p-norm
    T CalculateNorm(int p=2) const;
    void Print() const;
    void Print(std::string name) const;

};


#endif
