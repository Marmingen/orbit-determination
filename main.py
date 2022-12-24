#################################################################################
## IMPORTS

import math
from mathematics import *
from funcs import *

#################################################################################
## FUNCTIONS


def handle_orbital_calcs(r2, v2, times):
    h = r2.cross(v2)
    
    a = calc_a(r2, v2)
    e = calc_e(r2, v2, h)
    
    hec = Req2ec()*h
    eec = Req2ec()*e
    
    i = calc_i(hec)
    w = calc_w(hec, eec)
    W = calc_W(hec)
    T = calc_T(a, abs(e), r2, v2, times)
    
    return [a, abs(e), i, w, W, T]


def algorithm_orbital_elements(eq_coords, times, solar_pos):
    
    # step 1  || use terrestrial observations to calculate unit pos vec
    
    dhats = [calc_delta_hat(eq_coords[0]), calc_delta_hat(eq_coords[1]),\
                calc_delta_hat(eq_coords[2])]
    
    # step 2  || calculate the transformation matrix from equatorial to C
    
    xi = dhats[0]
    
    eta = calc_eta(xi,dhats[2])
    
    zeta = calc_zeta(xi,eta)
    
    R_e2C = Req2C(xi, eta, zeta)
    
    # step 3  || calculate some quantities and transform solar position
    
    transformed_solar_pos = [R_e2C*solpos for solpos in solar_pos]
    
    transformed_delta_hat2 = calc_delta_hatp(dhats[1], xi, eta, zeta)
    
    transformed_delta_hat3 = calc_delta_hatp(dhats[2], xi, eta, zeta)
    
    # START OF ITERATION
    # step 4  || | initial guess of the c-parameters
    # ---------- |
    # step 5  || | calculate the geocentric distances
    # ---------- |
    # step 6  || | evaluate the heliocentric positions
    # ---------- |
    # step 7  || | calculate the swept areas
    # ---------- |
    # step 8  || | recalculate the c-parameters
    # ---------- |
    # step 9  || | iterate 5-8 to convergence (if converges)
    
    [r, corr_times, iters] = step45678_iteration(times, dhats,\
        transformed_delta_hat2, transformed_delta_hat3, solar_pos, transformed_solar_pos)
    
    # step 10 || estimate the velocity vector for r2
    
    v2 = calc_v2(r, corr_times)
    
    # step 11 || calculate the orbital elements
    
    orbital_elements = handle_orbital_calcs(r[1], v2, times)
    
    for i in range(2,5):
        orbital_elements[i] = math.degrees(orbital_elements[i])
        if i > 2 and orbital_elements[i] < 0:
            orbital_elements[i] = orbital_elements[i] + 360
    
    print_results(orbital_elements, iters)


#########################################
## printing

def print_results(orbital_elements, iters):
    
    clear()
    print("CALCULATION INFORMATION")
    print(bar)
    print(f"{'differential epsilon: ':<15}{diff_eps}")
    print(f"{'area iterations: ':<15}{iters}")
    print()
    print("CALCULATED ORBITAL ELEMENTS")
    print(bar)
    for name, element, unit in zip(["a: ", "e: ", \
        "i: ", "\u03C9: ",\
        "\u03A9: ", "T: "],orbital_elements,\
        ["au", "", "◦", "◦", "◦", " JD"]):
        
        print(f"{name:<15}{round(element,5)}{unit}")


#################################################################################
## MAIN

def main():
    
    choice = input_selection()
    
    algorithm_orbital_elements(choice[0], choice[1], choice[2])


#################################################################################
## RUN CODE

if __name__ == "__main__":
    main()