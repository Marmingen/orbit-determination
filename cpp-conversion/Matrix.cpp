#include "Matrix.h"

Matrix::Matrix()
{
    r11 = 0, r12 = 0, r13 = 0;
    r21 = 0, r22 = 0, r23 = 0;
    r31 = 0, r32 = 0, r33 = 0;
    elements[3][3];
}

Matrix::Matrix
(
    double r11_in, double r12_in, double r13_in,
    double r21_in, double r22_in, double r23_in,
    double r31_in, double r32_in, double r33_in
)
{
    r11 = r11_in, r12 = r12_in, r13 = r13_in,
    r21 = r21_in, r22 = r22_in, r23 = r23_in,
    r31 = r31_in, r32 = r32_in, r33 = r33_in;
    double elements[3][3] = {{r11,r12,r13},{r21,r22,r23},{r31,r32,r33}};
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

