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

Vector operator+(const Vector &v1, const Vector &v2)
{
    Vector res;
    res.x = v1.x + v2.x;
    res.y = v1.y + v2.y;
    res.z = v1.z + v2.z;
    return res;
}

Vector operator-(const Vector &v1, const Vector &v2)
{
    Vector res;
    res.x = v1.x - v2.x;
    res.y = v1.y - v2.y;
    res.z = v1.z - v2.z;
    return res;
}

Vector Vector::operator=(Vector &v1)
{
    v1.x = x;
    v1.y = y;
    v1.z = z;
}

// vector multiplication defined as scalar
double operator*(const Vector &v1, const Vector &v2)
{
    return v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;
}

Vector operator*(const double &a, const Vector &v2)
{
    return Vector(a*v2.x, a*v2.y, a*v2.z);
}

Vector operator*(const Vector &v1, const double &b)
{
    return Vector(b*v1.x, b*v1.y, b*v1.z);
}

ostream Vector::operator<<(ostream &out, const Vector &v)
{
    return out << "[" << floor(v.x*100.0)/100.0 << ", " << floor(v.y*100.0)/100.0 << ", " << floor(v.z*100.0)/100.0 << "]";
}



// int main() 
// {
//     Vector v1(12,10,10);

//     cout << v1;
// }