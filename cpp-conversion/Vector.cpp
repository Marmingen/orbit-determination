#include "Vector.h"

Vector::Vector()
{   
    double x=0, y=0, z=0;
    double rowop = {0};
}

Vector::Vector(double in_x, double in_y, double in_z)
{
    x = in_x;
    y = in_y;
    z = in_z;
    rowop[0] = x;
    rowop[1] = y;
    rowop[2] = z;
}

Vector Vector::norm()
{
    Vector res;
    double a = abs();
    res.x = x/a;
    res.y = y/a;
    res.z = z/a;
    return res;
}

Vector Vector::cross(const Vector &v2)
{
    double x_temp = y*v2.z - z*v2.y;
    double y_temp = z*v2.x - x*v2.z;
    double z_temp = x*v2.y - y*v2.x;
    return Vector(x_temp, y_temp, z_temp);
}

double Vector::abs()
{
    return sqrt(pow(x,2)+pow(y,2)+pow(z,2));
}

Vector Vector::operator+(const Vector &v2)
{
    Vector res;
    res.x = x + v2.x;
    res.y = y + v2.y;
    res.z = z + v2.z;
    return res;
}

Vector Vector::operator-(const Vector &v2)
{
    Vector res;
    res.x = x - v2.x;
    res.y = y - v2.y;
    res.z = z - v2.z;
    return res;
}

Vector Vector::operator=(Vector &v2)
{
    v2.x = x;
    v2.y = y;
    v2.z = z;

    return v2;
}

// vector multiplication defined as scalar
double Vector::operator*(const Vector &v2)
{
    return x*v2.x+y*v2.y+z*v2.z;
}

Vector Vector::operator*(const double &a)
{
    return Vector(a*x, a*y, a*z);
}

std::ostream &operator<<(std::ostream &out, const Vector &v)
{
    out << "[" << floor(v.x*100.0)/100.0 << ", " << floor(v.y*100.0)/100.0 << ", " << floor(v.z*100.0)/100.0 << "]";
    return out;
}