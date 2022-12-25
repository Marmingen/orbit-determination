from constants import *
from mathematics import *
from .calculations import *
from time import sleep

#################################################################################
## input

def input_selection():  
    choices = {"dav":_davidsson, "man": _manaslu, "data": _input_custom,\
        "quit": _quit}
    
    while True:
    
        clear()
        print("SELECTION")
        print(bar)
        print(f"{'use data from 11798 Davidsson: ':{'.'}<42}{' dav'}")
        print(f"{'use data from 6918 Manaslu: ':{'.'}<42}{' man'}")
        print(f"{'enter data manually: ':{'.'}<41}{' data'}")
        print(f"{'quit the program: ':{'.'}<41}{' quit'}")
        print("")
        choice = input("input: ")
        
        if choice in choices.keys():
            return choices[choice]()
    
        else:
            print("incorrect format")
            sleep(1)

def _quit():
    return False

def _input_custom():
    clear()
    print("manual data entry")
    print(bar)
    
    name = input("name of object: ")

    print("enter the sun\'s geocentric position in the eccliptic")
    print(bar)
    
    solpos = [_solar_pos(i) for i in range(1,4)]
    
    print("enter the observational dates")
    print(bar)
    times = [float(input(f"date {i} [JD]: ")) for i in range(1,4)]
    
    
    print()
    print("enter the observed angular coordinates")
    print(bar)
    eq_coords = [_angular_coords(i) for i in range(1,4)]
    
    return [eq_coords, times, solpos, name]

def _solar_pos(i):
    print(f"enter position {i}")
    print(bar)
    x = float(input("x-coordinate [au]: "))
    y = float(input("y-coordinate [au]: "))
    z = float(input("z-coordinate [au]: "))
    print()
    return Vector(x,y,z)

def _angular_coords(i):
    print(f"coordinate set {i}")
    print(bar)
    alpha = input("right ascension [h m s]: ").split(" ")
    delta = input("declination [degrees arcminutes arcseconds]: ").split(" ")
    print()
    
    alpha = [float(alp) for alp in alpha]
    delta = [float(delt) for delt in delta]
    
    return [time_to_rad(alpha[0], alpha[1], alpha[2]),\
        deg_to_rad(delta[0], delta[1], delta[2])]
    

#################################################################################

def _davidsson():
    
    # sun geocentric position in eccliptic
    sol1 = Vector(0.1765528,-0.8875529,-0.3847901)
    sol2 = Vector(0.2786189,-0.8652295,-0.3751088)
    sol3 = Vector(0.5026422,-0.7762153,-0.3365189)
    
    # observation times
    times = [2453736.5, 2453742.5, 2453756.5]
    
    # Terrestrial observations as angular coordinates
    first = [time_to_rad(h=15, m=58, s=50.5), deg_to_rad(deg=-29, min=55, sec=38)]
    
    second = [time_to_rad(h=16, m=8, s=17.1), deg_to_rad(deg=-30, min=35, sec=20)]
    
    third = [time_to_rad(h=16, m=29, s=46.3), deg_to_rad(deg=-32, min=2, sec=34)]
    
    eq_coords = [first, second, third]
    
    return [eq_coords, times, [sol1, sol2, sol3], "Davidsson"]

def _manaslu():
    
    # sun geocentric position in eccliptic
    sol1 = Vector(0.8814882, -0.41157715846813614, -0.17841808859776345)    
    sol2 = Vector(0.9953944, 0.048270359200408924, 0.020917519100428113)
    sol3 = Vector(0.8475964, 0.4959487672703848, 0.21497992247295952)
    
    # observation times
    times = [2459632.5, 2459662.5, 2459692.5]
    
    # Terrestrial observations as angular coordinates
    first = [time_to_rad(h=3, m=18, s=34), deg_to_rad(deg=16, min=26, sec=10)]
    
    second = [time_to_rad(h=4, m=1, s=51), deg_to_rad(deg=19, min=1, sec=58)]
    
    third = [time_to_rad(h=4, m=55, s=19.1), deg_to_rad(deg=21, min=8, sec=14)]
    
    eq_coords = [first, second, third]
    
    return [eq_coords, times, [sol1, sol2, sol3], "Manaslu"]