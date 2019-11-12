#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string>
#include "include/Matrix.h"
#include "include/Vector.h"
#include "include/kernel.h"

using namespace std;

typedef double T;

int main(void){
    int cols, rows;
    rows = 4;
    cols = rows;
    T alpha = 2.0;
    T** Ar;

    Vector<T> v(rows,T(1));
    Matrix<T> A(rows,cols,T(1));
    Matrix<T> C(v);
    Matrix<T> B = C.Scale(alpha);
    Matrix<T> D = C.Multiply(B);
    Vector<T> matvec = C.Multiply(v);

    SquaredExponential<T> K;

    double params[2];
    params[0] = 1; params[1] = 1;

    v[1] = 2.4;

    v.Print();

    double res = K.Kernel(params, 1,1);
    cout << "result SE Kernel = " << res << endl;

    Matrix<T> B = A.Multiply(A);
    A.Print("A = ");
    B.Print("B = ");
    C.Print("C = ");
    D.Print("D = ");
    Ar = A.RawPointer();
    matvec.Print("Av =");

    cout << v.Dot(v) << endl;


    return 0;
}