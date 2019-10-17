# Simple calculator to convert betwen RA, DEC and l,m
import sys
import math
import numpy as np

def eq_dcc(d,ra,d_0,ra_0):
    d_ra = ra-ra_0
    l = math.cos(d) * math.sin(d_ra)
    m = math.sin(d) * math.cos(d_0) - math.cos(d) * math.sin(d_0) * math.cos(d_ra)
    return (l,m)

def dcc_eq(l,m,d_0,ra_0):
    d = math.asin(m * math.cos(d_0) + math.sin(d_0) * math.sqrt(1 - l**2 - m**2))
    ra = ra_0 + math.atan(1/(math.acos(d_0) * math.sqrt(1 - l ** 2 - m ** 2) - m * math.sin(d_0)))
    return (d,ra)

if __name__ == "__main__":
    if sys.argv[1] == "eq":
        print(dcc_eq(float(sys.argv[2]) * np.pi/180, float(sys.argv[3]) * np.pi/180, float(sys.argv[4]) * np.pi/180, float(sys.argv[5]) * np.pi/180))
    elif sys.argv[1] == "dcc":
        print(eq_dcc(float(sys.argv[2]) * np.pi/180, float(sys.argv[3]) * np.pi/12, float(sys.argv[4]) * np.pi/180, float(sys.argv[5]) * np.pi/12))
