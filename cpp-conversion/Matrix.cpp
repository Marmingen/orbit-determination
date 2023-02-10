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
};

Matrix::Matrix()
{
    r11 = 0, r12 = 0, r13 = 0;
    r21 = 0, r22 = 0, r23 = 0;
    r31 = 0, r32 = 0, r33 = 0;
    double elements = {0};
}

Vector operator*(const Vector &v1, const Matrix &M)
{
    Vector res;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            res.rowop[i] = v1.rowop[i]*M.elements[j][i];
        }
    }
    res.x = res.rowop[0];
    res.y = res.rowop[1];
    res.z = res.rowop[2];
    return res;
}

Vector operator*(const Matrix &M, const Vector &v2)
{
    Vector res;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            res.rowop[i] = v2.rowop[j]*M.elements[i][j];
        }
    }
    res.x = res.rowop[0];
    res.y = res.rowop[1];
    res.z = res.rowop[2];
    return res;
}

int main()
{
    Matrix M1;

    Vector v1(1,1,1);

    cout << (v1*M1).x;
}
