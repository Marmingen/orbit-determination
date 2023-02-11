#include <iostream>
#include <cmath>
#include <math.h>
#include <string>
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

void calcC(double cvals[], const double times[], const double y[])
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


void calcYs(double y[3], int &iters, const double times[3],const Vector r[3])
{
    // y1 = y23
    double y1, y2, y3;
    int iters1, iters2, iters3;

    calcYab(iters1, y1, times[1], times[2], r[1], r[2]);
    calcYab(iters2, y2, times[0], times[2], r[0], r[2]);
    calcYab(iters3, y3, times[0], times[1], r[0], r[1]);

    // returns
    y[0] = y1;
    y[1] = y2;
    y[2] = y3;
    iters = iters1 + iters2 + iters3;
}

void calcYab(int &iters, double &y3, const double ta, const double tb, const Vector ra, const Vector rb)
{
    /*
    * calculates the y-value for yab
    * y1 = y23
    * y2 = y13
    * y3 = y12
    */

    double K = sqrt(2*(abs(ra)*abs(rb)+ra*rb));

    double m2 = pow((k*(tb-ta)),2)/pow(K,3);

    double l = (abs(ra)+abs(rb)-K)/(2*K);

    double y1 = 1;

    double y2;

    for (int i = 0; i < maxIters; i++)
    {
        steffensen(y1, y2, y3, m2, l);

        double ystar = y1 - 2*y2 + y3;

        if (abs(ystar) <= diffEps)
        {
            return;
        }

        y1 = y1 - pow(y2-y1,2)/ystar;

        iters++;
    }

    throw runtime_error("Not converging, too many iterations");

}

void steffensen(const double y1, double &y2, double &y3, const double m2, const double l)
{
    double x1 = m2/pow(y1,2)-l;

    y2 = 1 + (4/3) * (m2/pow(y1,2))*Q(x1);

    double x2 = m2/pow(y2,2) - l;

    y3 = 1 + (4/3)*(m2/pow(y2,2))*Q(x2);
}

double Q(const double x)
{
    if (0 < x <= 1/2)
    {
        return 3/16*(2*(2*x-1)*sqrt(x-pow(x,2))+asin(2*x-1)+M_PI/2)/(pow(x-pow(x,2),3/2));
    }
    else if (x < 0)
    {
        return 3/16*(2*(1-2*x)*sqrt(pow(x,2)-x)-log(1-2*x+2*sqrt(pow(x,2)-x)))/(pow(pow(x,2)-x,3/2));
    }
    else {return 1;}   
}

double calcV2(const Vector r[3], const double times[3])
{
    double f1 = 1-1/2*(pow(k,2)*pow(times[0]-times[1],2)*pow(abs(r[1]),-3));
    
    double f3 = 1-1/2*(pow(k,2)*pow(times[2]-times[1],2)*pow(abs(r[1]),-3));
    
    double g1 = times[0] - times[1] - 1/6*(pow(k,2)*pow(times[0]-times[1],3)*pow(abs(r[1]),-3));
    
    double g3 = times[2] - times[1] - 1/6*(pow(k,2)*pow(times[2]-times[1],3)*pow(abs(r[1]),-3));
    
    double v2_1 = 1/g1 * (r[0]-f1*r[1]);
    
    double v2_2 = 1/g3 * (r[2]-f3*r[1]);
    
    return (v2_1 + v2_2)*(1/2);
}

//########## element calcs ############

double calcSemimajor(const Vector r2, const Vector v2)
{
    double v = abs(v2);
    double r = abs(r2);

    return 1/pow(2/r-(v/k),2);
}

Vector calcEcc(const Vector r2, const Vector v2, const Vector h)
{
    return 1/pow(k,2)*v2.cross(h) - r2/abs(r2);
}

double calcInc(const Vector h)
{
    return acos(h.z/h.abs());
}

double calcNode(const Vector h)
{
    return atan2(h.x, -h.y);
}

double calcArg(const Vector h, const Vector e)
{
    return atan2(h.abs()*e.z, h.x*e.y-h.y*e.x);
}

double calcT(const double a, const double e, const double r, const double v, const double times[3])
{
    double x = (1-abs(r)/a)/e;
    double y = (r*v)/sqrt(a*pow(k,2))/e;

    double E = atan2(y,x);
    double M = E - e*sin(E);
    double n = sqrt(pow(k,2)/pow(a,3));

    return (times[1]-M/n);
}

void calcIters(const double ogTimes[3], const Vector dHats[3], const Vector transCrds2, const Vector transCrds3, const Vector solarDists[3], const Vector transSolPos[3])
{
    class RollingVals
    {
        public:
            double values[3];

            RollingVals()
            {
                double values[3] = {1,2,0};
            };

            double rollValues(const double newVal)
            {
                double temp[3] = {values[1], values[2], newVal};
                delete values;
                *values = *temp;
            };
    };

    RollingVals rollVals;

    double y[3] = {1,1,1};

    for (int i = 0; i < 100; i++)
    {
        // c = [c1, c3]
        double c[2];
        calcC(c, ogTimes, y);

        double deltas[3];
        calcGeoDist(deltas, transCrds2, transCrds3, transSolPos, c);

        double deltaTimes[3];
        planetaryAbberation(deltaTimes, deltas);

        double corrTimes[3];
        for (int i = 0; i < 3; i++)
        {
            corrTimes[i] = ogTimes[i] - deltaTimes[i];
        }

        // c = [c1, c3]
        calcC(c, corrTimes, y);

        Vector r[3];
        heliocentricDist(r, deltas, dHats, solarDists);

        int aIters = 0;
        calcYs(y, aIters, corrTimes, r);

        rollVals.rollValues(r[1].abs());
    }

    if (abs(rollVals.values[0]-2*rollVals.values[1]+rollVals.values[2]) < diffEps)
    {
        return;
    }

    throw runtime_error("Positions not converging, most likely lie on a great circle");
}

double nrmlToJDT(string date, const double corr)
{
    int Y = stoi(date.substr(0,4));
    int M = stoi(date.substr(4,6));
    int D = stoi(date.substr(6));

    double JD = 367*Y - int(7/4*(Y + int((M+9)/12)));
    D += int(275*M/9);
    JD += D + 1721013.5;
    JD -=corr/24;
}