#################################################################################
## IMPORTS

import math
from mathematics import *
from constants import *

#########################################
## calculations

def deg_to_rad(deg, min, sec):
    
    sgn = deg//abs(deg)
    
    ang = 0
    
    ang += sgn*sec/3600
    
    ang += sgn*min/60
    
    ang += deg
      
    return ang/180*math.pi

"""
converts radians to hours, mins, and secs
(one revolution is 24 hours)
"""
def time_to_rad(h, m, s):

    ang = 0
    
    ang += s*2*math.pi/24/3600
    
    ang += m*2*math.pi/24/60
    
    ang += h*2*math.pi/24
    
    return ang


def calc_delta_hat(eq_crds):
    """
    calculates the unit position vectors of the object in the 
    eq-system as seen from Earth
    
    the input will be of the format:
    eq_crds = [alpha, delta]
            = [right ascension, declination]
    """
    x = math.cos(eq_crds[0])*math.cos(eq_crds[1])
    y = math.sin(eq_crds[0])*math.cos(eq_crds[1])
    z = math.sin(eq_crds[1])
    
    return Vector(x,y,z)

# calculate specific xi, eta, zeta values
def calc_delta_hatp(delta_hat, xi, eta, zeta):
    return Vector(xi*delta_hat, eta*delta_hat, zeta*delta_hat)


def calc_eta(delta1=Vector, delta3=Vector):
    eta = delta1.cross(delta3.cross(delta1))
    eta = eta/abs(eta)
    
    return eta

# calculate 
def calc_zeta(delta1, eta):
    return delta1.cross(eta)

#########################################
## matrices

def Req2C(xi, eta, zeta):
    """
    calculates the rotational matrix from eq-system to
    C-system, used for solar and unit pos
    """
    return Matrix("Req2C",\
        xi.x, xi.y, xi.z,\
        eta.x, eta.y, eta.z,\
        zeta.x, zeta.y, zeta.z)
    
    
def Req2ec():
    return Matrix("Req2ec",\
        1,0,0,\
        0, math.cos(eps), math.sin(eps),\
        0, -math.sin(eps), math.cos(eps))


def R_ec2eq(deg):
    return Matrix("R_ec2eq",\
        1, 0, 0,\
        0, math.cos(deg),-math.sin(deg),\
        0, math.sin(deg), math.cos(deg))
    
    
def calc_c(times, y = [1,1,1]):
    # constants = [c1, c3]
    
    #     y2/y1(t3-t2)/(t3-t1)
    c1 = (y[1]/y[0])*(times[2]-times[1])/(times[2]-times[0])
    
    #     y2/y3(t2-t1)/(t3-t1)
    c3 = (y[1]/y[2])*(times[1]-times[0])/(times[2]-times[0])
    
    return [c1, c3]


def calc_geo_dist(trans_crds2, trans_crds3, solar_dists, c):
    """
    calculating the geocentric distances,
    
    the input has format:
    trans_crds2 = [xi2, eta2, zeta2]
    trans_crds3 = [xi3, eta3, 0]
    """
    
    #         -c_1  R'_1z           +R'_2z           -c3   R'_3z             zeta2
    delta2 = (-c[0]*solar_dists[0].z+solar_dists[1].z-c[1]*solar_dists[2].z)/trans_crds2.z
    
    #         D2     eta2          c_1  R'_1y            R'_2y            c_3  R'_3y              c_3  eta3
    delta3 = (delta2*trans_crds2.y+c[0]*solar_dists[0].y-solar_dists[1].y+c[1]*solar_dists[2].y)/(c[1]*trans_crds3.y)
    
    #         D2     xi2           c_3  D3     xi3           c_1  R'_1x            R'_2x            c_3  R'_3x             c_1
    delta1 = (delta2*trans_crds2.x-c[1]*delta3*trans_crds3.x+c[0]*solar_dists[0].x-solar_dists[1].x+c[1]*solar_dists[2].x)/c[0]
    
    return [delta1, delta2, delta3]

# due to special relativity
def planetary_abberation(deltas):
    correction = [abs(delta)*au/c/3600/24 for delta in deltas]
    return correction

# estimates the heliocentric distance
def heliocentric_dist(deltas, dhats, solar_pos):
    r = [delta*dhat-solpos for delta, dhat, solpos in zip(deltas, dhats, solar_pos)]
    return r


def calc_ys(times, r):
    """
    this handles the calculations for step 7
    i.e. estimating the swept areas
    """
    # y1 = y23
    y1 = calc_yab(times[1], times[2], r[1], r[2])
    
    # y2 = y13
    y2 = calc_yab(times[0], times[2], r[0], r[2])
    
    # y3 = y12
    y3 = calc_yab(times[0], times[1], r[0], r[1])
    
    return [y1, y2, y3]


def calc_yab(ta, tb, ra, rb):
    """
    calculates the y-value for yab
    y1 = y23
    y2 = y13
    y3 = y12
    """
        
    K = math.sqrt(2*(abs(ra)*abs(rb) + ra*rb))
    
    m2 = (k*(tb-ta))**2/(K**3)
    
    l = (abs(ra) + abs(rb) - K)/(2*K)
    
    y1 = 1
    
    iters = 0
    
    for _ in range(100):
        
        [y2, y3] = steffensen(m2, l, y1)
        
        ystar = y1 - 2*y2 + y3 
        
        if abs(ystar) <= diff_eps:
            return y3
        
        y1 = y1 - (y2-y1)**2/ystar
        
        iters+=1
    
    raise Exception(f"Not converging, {iters} iterations")


def steffensen(m2, l, y1):
    x1 = m2/y1**2-l
    
    y2 = 1 + (m2/y1**2)*4/3*Q(x1)
    
    x2 = m2/y2**2 - l
    
    y3 = 1 + (m2/y2**2)*4/3*Q(x2)
    
    return [y2, y3]
    
    
def Q(x):
    if 0 < x <= 1/2:
        return 3/16*(2*(2*x-1)*math.sqrt(x-x**2)+math.asin(2*x-1)+math.pi/2)/((x-x**2)**(3/2))
    elif x < 0:
        return 3/16*(2*(1-2*x)*math.sqrt(x**2-x)-math.log(1-2*x+2*math.sqrt(x**2-x)))/((x**2-x)**(3/2))
    else:
        return 1  
    
    
def calc_v2(r, times):
    f1 = 1-1/2*(k**2*(times[0]-times[1])**2*abs(r[1])**(-3))
    
    f3 = 1-1/2*(k**2*(times[2]-times[1])**2*abs(r[1])**(-3))
    
    g1 = times[0] - times[1] - 1/6*(k**2*(times[0]-times[1])**3*abs(r[1])**(-3))
    
    g3 = times[2] - times[1] - 1/6*(k**2*(times[2]-times[1])**3*abs(r[1])**(-3))
    
    v2_1 = 1/g1 * (r[0]-f1*r[1])
    
    v2_2 = 1/g3 * (r[2]-f3*r[1])
    
    return (v2_1 + v2_2)*(1/2)

# calculate the semimajor axis
def calc_a(r2,v2):
    
    v = abs(v2)
    r = abs(r2)

    return 1/(2/r-(v/k)**2)
    
# calculate the eccentricity vector
def calc_e(r2, v2, h):
    return 1/k**2 * v2.cross(h) - r2/abs(r2)

# calculate the inclination (h must be in equatorial)
def calc_i(h):
    return math.acos(h.z/abs(h))

# calculate the longitude of the ascending node
def calc_W(h):
    return math.atan2(h.x,-h.y)

# calculate the argument of perihelion (h,e must be in ecliptic)
def calc_w(h, e):    
    return math.atan2(abs(h)*e.z, h.x*e.y-h.y*e.x)

# calculate the time since last perihelion
def calc_T(a,e,r,v,times):
    
    x = (1-abs(r)/a)/e
    
    y = (r*v)/math.sqrt(a*k**2)/e
    
    E = math.atan2(y,x)
    
    M = E - e*math.sin(E)
    
    n = math.sqrt(k**2/a**3)
    
    T = times[1]-M/n
    
    return T

def step45678_iteration(ogtimes, dhats, trans_crds2, trans_crds3, solar_dists, trans_solpos):    
    
    def roll(lis, val):
        lis.pop(0)
        return lis.append(val)
    
    rolling_values = [1,2,0]
    
    y = [1,1,1]
    
    for iters in range(100):
    
        # c = [c1, c3]
        c = calc_c(ogtimes, y)

        deltas = calc_geo_dist(trans_crds2, trans_crds3, trans_solpos, c)
        
        delta_times = planetary_abberation(deltas)
        
        corr_times = [ogtime - corr for ogtime,corr in zip(ogtimes, delta_times)]
        
        r = heliocentric_dist(deltas, dhats, solar_dists)
        
        y = calc_ys(corr_times, r)
        
        roll(rolling_values, abs(r[1]))
        
        if abs(rolling_values[0]-2*rolling_values[1]+rolling_values[2]) < diff_eps:            
            return [r, corr_times, iters]
        
    raise Exception("Selected positions are probably on a great circle")
        
    
def nrml_to_JDT(date, corr=0):
    Y = int(date[0:4])
    M = int(date[4:6])
    D = int(date[6:].strip())
    # TT is assumed to be 00 hrs
    
    JD = 367*Y - int(7/4*(Y + int((M+9)/12)))
    JD += int(275*M/9)
    JD += D + 1721013.5
    JD -=corr/24

    return JD