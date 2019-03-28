import numpy as np
import matplotlib.pyplot as plt
import scipy.constants
import math
from mpl_toolkits.mplot3d import Axes3D
import plotBL
from matplotlib.patches import Ellipse

def plot_baseline(b_ENU, L, f, ant1, ant2):
    D = math.sqrt(np.sum((b_ENU)**2))
    A = np.arctan2(b_ENU[0],b_ENU[1])
    E = np.arcsin(b_ENU[2]/D)
    # plotBL.sphere(ant1,ant2,A,E,D,L)
    B = np.array([D * (math.cos(L)*math.sin(E) - math.sin(L) * math.cos(E)*math.cos(A)),
                D * (math.cos(E)*math.sin(A)),
                D * (math.sin(L)*math.sin(E) + math.cos(L) * math.cos(E)*math.cos(A))])
    c = scipy.constants.c                                        # Speed of light
    lam = c/f
    # dec = ?
    time_steps = 600
    h = np.linspace(-12,12,num=time_steps)*np.pi/12
    X = B[0]
    Y = B[1]
    Z = B[2]
    u = lam**(-1)*(np.sin(h)*X+np.cos(h)*Y)
    v = lam**(-1)*(-np.sin(dec)*np.cos(h)*X+np.sin(dec)*np.sin(h)*Y+np.cos(dec)*Z)
    w = lam**(-1)*(np.cos(dec)*np.cos(h)*X-np.cos(dec)*np.sin(h)*Y+np.sin(dec)*Z)
    # plotBL.UV(u,v,w)
    a=np.sqrt(X**2+Y**2)/lam # major axis
    b=a*np.sin(dec)              # minor axis
    v0=(Z/lam)*np.cos(dec)  # center of ellipse
    plotBL.UVellipse(u,v,w,a,b,v0)

def plot_array(antennas):
    plt.scatter(antennas[:,0], antennas[:,1])
    plt.grid(True)
    plt.xlabel('E-W [m]')
    plt.ylabel('N-S [m]')
    plt.title('TART Array Layout')
    plt.show()
