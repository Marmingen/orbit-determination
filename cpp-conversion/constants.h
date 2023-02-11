#include <cmath>

// celestial instance constants
// 23.439 at J2000
// 23.436 in 2022
const double eps = 23.439*M_PI/180;     //rad
const double J2000 = 2451545.0;

// physical constants
const double k = 0.01720209895;   // rad/day
const double au = 1.495978707e11; // m
const double c =  2.99792458e8; // m/s

// computer stuff
const double diffEps = 1e-12;  // min diff
const int maxIters = 100;