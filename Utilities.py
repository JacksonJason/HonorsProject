import numpy as np
import matplotlib.pyplot as plt
import scipy.constants
import math
from mpl_toolkits.mplot3d import Axes3D
import plotBL
from matplotlib.patches import Ellipse
import matplotlib.pylab as pl

h = np.linspace(-12,12,num=600)*np.pi/12
dec = 0

# def draw_matrix(matrix):
#     print(matrix)
#     # plt.figure()
#     tb = pl.table(cellText=matrix, loc=(0,0), cellLoc='center')
#     tc = tb.properties()['child_artists']
#     # for cell in tc:
#     #     cell.set_height(1/len(matrix))
#     #     cell.set_width(1/len(matrix))
#
#     ax = pl.gca()
#     ax.set_xticks([])
#     ax.set_yticks([])
#     plt.show()

def get_B(b_ENU, L):
    D = math.sqrt(np.sum((b_ENU)**2))
    A = np.arctan2(b_ENU[0],b_ENU[1])
    E = np.arcsin(b_ENU[2]/D)
    # plotBL.sphere(ant1,ant2,A,E,D,L)
    B = np.array([D * (math.cos(L)*math.sin(E) - math.sin(L) * math.cos(E)*math.cos(A)),
                D * (math.cos(E)*math.sin(A)),
                D * (math.sin(L)*math.sin(E) + math.cos(L) * math.cos(E)*math.cos(A))])
    return B

def get_lambda(f):
    c = scipy.constants.c                                        # Speed of light
    lam = c/f
    return lam

def plot_baseline(b_ENU, L, f, ant1, ant2):
    B = get_B(b_ENU, L)

    lam = get_lambda(f)
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

def plot_visibilities(u, v, b_ENU, L, f):
    B = get_B(b_ENU, L)
    lam = get_lambda(f)
    d = 1/lam * B[1]
    u_d = d*np.cos(h)
    v_d = d*np.sin(h)*np.sin(dec)
    # X = B[0]
    # Y = B[1]
    # Z = B[2]
    # u = lam**(-1)*(np.sin(h)*X+np.cos(h)*Y)
    # v = lam**(-1)*(-np.sin(dec)*np.cos(h)*X+np.sin(dec)*np.sin(h)*Y+np.cos(dec)*Z)
    # w = lam**(-1)*(np.cos(dec)*np.cos(h)*X-np.cos(dec)*np.sin(h)*Y+np.sin(dec)*Z)
    uu, vv = np.meshgrid(u, v)
    zz = np.zeros(uu.shape).astype(complex)
    s = point_sources.shape
    for counter in range(1, s[0]+1):
        A_i = point_sources[counter-1,0]
        l_i = point_sources[counter-1,1]
        m_i = point_sources[counter-1,2]
        zz += A_i*np.exp(-2*np.pi*1j*(uu*l_i+vv*m_i))
    zz = zz[:,::-1]

    plt.figure()
    plt.subplot(121)
    plt.imshow(zz.real,extent=[-1*(np.amax(np.abs(u_d)))-10, np.amax(np.abs(u_d))+10,-1*(np.amax(abs(v_d)))-10, \
                               np.amax(abs(v_d))+10])
    plt.plot(u_d,v_d,"k")
    plt.xlim([-1*(np.amax(np.abs(u_d)))-10, np.amax(np.abs(u_d))+10])
    plt.ylim(-1*(np.amax(abs(v_d)))-10, np.amax(abs(v_d))+10)
    plt.xlabel("u")
    plt.ylabel("v")
    plt.title("Real part of visibilities")

    plt.subplot(122)
    plt.imshow(zz.imag,extent=[-1*(np.amax(np.abs(u_d)))-10, np.amax(np.abs(u_d))+10,-1*(np.amax(abs(v_d)))-10, \
                               np.amax(abs(v_d))+10])
    plt.plot(u_d,v_d,"k")
    plt.xlim([-1*(np.amax(np.abs(u_d)))-10, np.amax(np.abs(u_d))+10])
    plt.ylim(-1*(np.amax(abs(v_d)))-10, np.amax(abs(v_d))+10)
    plt.xlabel("u")
    plt.ylabel("v")
    plt.title("Imaginary part of visibilities")
