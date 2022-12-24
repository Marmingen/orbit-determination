from constants import *
from mathematics import *
from .calculations import *

#################################################################################
## input

def input_selection():  
    choices = {"dav":_davidsson, "man": _manaslu}
    
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
    

#################################################################################

def _davidsson():
    sol1 = Vector(0.1765528,-0.8875529,-0.3847901)
    sol2 = Vector(0.2786189,-0.8652295,-0.3751088)
    sol3 = Vector(0.5026422,-0.7762153,-0.3365189)
    
    times = [2453736.5, 2453742.5, 2453756.5]
    
    first = [time_to_rad(h=15, m=58, s=50.5), deg_to_rad(deg=-29, min=55, sec=38)]
    
    second = [time_to_rad(h=16, m=8, s=17.1), deg_to_rad(deg=-30, min=35, sec=20)]
    
    third = [time_to_rad(h=16, m=29, s=46.3), deg_to_rad(deg=-32, min=2, sec=34)]
    
    eq_coords = [first, second, third]
    
    return [eq_coords, times, [sol1, sol2, sol3]]

def _manaslu():
    sol1 = Vector(0.8814882, -0.41157715846813614, -0.17841808859776345)    
    sol2 = Vector(0.9953944, 0.048270359200408924, 0.020917519100428113)
    sol3 = Vector(0.8475964, 0.4959487672703848, 0.21497992247295952)

    # sol1 = Vector(0.8815, -0.44859, 0.0000183)
    
    # sol3 = Vector(0.9954, 0.0526, 0)

    times = [2459632.5, 2459662.5, 2459692.5]
    
    first = [time_to_rad(h=3, m=18, s=34), deg_to_rad(deg=16, min=26, sec=10)]
    
    second = [time_to_rad(h=4, m=1, s=51), deg_to_rad(deg=19, min=1, sec=58)]
    
    # second = [time_to_rad(h=3, m=40, s=8.6), deg_to_rad(deg=17, min=50, sec=1)]
    
    third = [time_to_rad(h=4, m=55, s=19.1), deg_to_rad(deg=21, min=8, sec=14)]
    
    # times[1] = 2459648.5
    
    # sol2 = Vector(9.750096102174219E-01,-1.879539702586240E-01,5.023082231012434E-06)
    
    # sol2 = Vector(0.975, -0.18795, 0)
    
    eq_coords = [first, second, third]
    
    return [eq_coords, times, [sol1, sol2, sol3]]