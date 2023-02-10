#include <iostream>
#include <cmath>
#include <math.h>
#include "Vector.h"
#include "Matrix.h"
#include "constants.h"


//########## angle calcs ############
double degToRad(const double deg, const double min, const double sec)
{
    int sgn = deg/abs(deg);
    double ang = 0;
    ang += sgn*sec/3600;
    ang += sgn*min/60;
    ang += deg;

    return ang*M_PI/180;
}

double timeToRad(const int h, const int m, const double s)
{
    double ang = 0;
    ang += s*2*M_PI/24/3600;
    ang += m*2*M_PI/24/60;
    ang += h*2*M_PI/24;

    return ang;
}

//########## vector calcs ############
Vector calcDeltaHat(const double eqCrds[2])
{
    double x = cos(eqCrds[0])*cos(eqCrds[1]);
    double y = sin(eqCrds[0])*cos(eqCrds[1]);
    double z = sin(eqCrds[1]);

    return Vector(x,y,z);
}

Vector calcDeltaHatP(const double deltaHat, const double xi, const double eta,  double zeta)
{
    return Vector(xi*deltaHat, eta*deltaHat, zeta*deltaHat);
}

Vector calcEta(Vector &delta1, Vector &delta3)
{
    Vector eta = (delta1.cross(delta3.cross(delta1)));
    return eta.norm();
}

Vector calcZeta(Vector &delta1, Vector &eta)
{
    return delta1.cross(eta);
}

//########## matrix calcs ############

Matrix ReqToC(const Vector &xi, const Vector &eta, const Vector &zeta)
{
    return Matrix
    (
        xi.x, xi.y, xi.z,
        eta.x, eta.y, eta.z,
        zeta.x, zeta.y, zeta.z
    );
}

Matrix ReqToec(const double deg)
{
    return Matrix
    (
        1, 0, 0,
        0, cos(deg), sin(deg),
        0, -sin(deg), cos(deg)
    );
}

Matrix RecToeq(const double deg)
{
    return Matrix
    (
        1, 0, 0,
        0, cos(deg), -sin(deg),
        0, sin(deg), cos(deg)
    );
}

//########## element calcs ############

void calcC(double cvals[], const float times[], const double y[])
{
    cvals[0] = (y[1]/y[0])*(times[2]-times[1])/(times[2]-times[0]);
    cvals[1] = (y[1]/y[2])*(times[1]-times[0])/(times[2]-times[0]);
}

void calcGeoDist(double deltas[3], const Vector transCrds2, const Vector transCrds3, const Vector solarDists[], const double c[])
{
    //           -c_1  R'_1z           +R'_2z         -c3   R'_3z             zeta2
    deltas[1] = (-c[0]*solarDists[0].z+solarDists[1].z-c[1]*solarDists[2].z)/transCrds2.z;

    //           D2        eta2         c_1  R'_1y            R'_2y          c_3  R'_3y              c_3  eta3
    deltas[2] = (deltas[1]*transCrds2.y+c[0]*solarDists[0].y-solarDists[1].y+c[1]*solarDists[2].y)/(c[1]*transCrds3.y);
    
    //           D2        xi2           c_3  D3       xi3           c_1  R'_1x            R'_2x       c_3  R'_3x             c_1
    deltas[0] = (deltas[1]*transCrds2.x-c[1]*deltas[2]*transCrds3.x+c[0]*solarDists[0].x-solarDists[1].x+c[1]*solarDists[2].x)/c[0];
}

void planetaryAbberation(double corrDeltas[3], const double deltas[3])
{
    for (int i = 0; i < 3; i++)
    {
        corrDeltas[i] = abs(deltas[i])*au/c/3600/24;
    }
}

void heliocentricDist(Vector r[3], double deltas[3], Vector dhats[3], const Vector solarPos[3])
{
    for (int i = 0; i < 3; i++)
    {
        r[i] = deltas[i]*dhats[i]-solarPos[i];
    }
}