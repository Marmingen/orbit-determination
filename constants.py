import math
import os

#################################################################################
## GLOBAL CONSTANTS

# Meta usage
bar = "**********************************************"

# Celestial instance constants
# 23.439 at J2000
# 23.436 at 2022
eps = 23.439*math.pi/180    #rad
J2000 = 2451545.0

# Physical constants
G = 6.674e-11 # m**3/kg/s**2
k = 0.01720209895   # rad/day
au = 1.495978707e11 # m
c =  2.99792458e8 # m/s

# Lambda funcs
clear = lambda: os.system('cls')

# error
# diff_eps = 0.0000001

diff_eps = 1e-10