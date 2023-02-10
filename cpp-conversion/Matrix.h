#include <iostream>
#include "Vector.h"

using namespace std;

class Matrix
{
    public:
        double r11, r12, r13;
        double r21, r22, r23;
        double r31, r32, r33;
        double elements[3][3];
        Matrix();
        Matrix
        (   double r11_in, double r12_in, double r13_in,
            double r21_in, double r22_in, double r23_in,
            double r31_in, double r32_in, double r33_in);

        // operators
};

Vector operator*(const Vector &v1, const Matrix &M);
Vector operator*(const Matrix &M, const Vector &v2);