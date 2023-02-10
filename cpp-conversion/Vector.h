#include <iostream>
#include <cmath>
using namespace std;

class Vector 
{
    // private:
        
    public:
        double x, y, z;
        double rowop[3];
        Vector();
        Vector(double in_x, double in_y, double in_z);
        Vector norm();
        Vector cross(const Vector &v2);
        double abs();

        // operators
        Vector operator=(Vector &v2);
        Vector operator+(const Vector &v2);
        Vector operator-(const Vector &v2);
        double operator*(const Vector &v2);
        Vector operator*(const double &a);
};

ostream &operator<<(ostream &out, const Vector &v);