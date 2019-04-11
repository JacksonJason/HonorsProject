import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.constants
import math
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Ellipse
import matplotlib.pylab as pl

h = np.linspace(-12,12,num=600)*np.pi/12
dec = 0

def draw_matrix(matrix):
    plt.figure()
    plt.subplot(121)
    plt.imshow(matrix.real)
    # plt.xlabel("Timeslots")
    # plt.ylabel("Jy")
    plt.title("Real: visibilities")

    plt.subplot(122)
    plt.imshow(matrix.imag)
    # plt.xlabel("Timeslots")
    # plt.ylabel("Jy")
    plt.title("Imag: visibilities")
    plt.savefig('Plots/Antenna_Visibilities.png')
    plt.close()

def tabulate_matrix(matrix):
    n = []
    for i in range(matrix.shape[0]):
        m = []
        for j in range(matrix.shape[1]):
            r = '%.6f' % matrix[i][j].real
            im = '%.6f' % matrix[i][j].imag
            if "-" in im:
                m.append(r + im + "i")
            else:
                m.append(r + "+" + im + "i")
        n.append(m)
    m = np.array(n)

    fig=plt.figure(figsize=(m.shape[0]*1, m.shape[1]*1))
    ax = fig.add_subplot(111)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    the_table = ax.table(cellText=m,
              loc='center')
    the_table.set_zorder(10)
    plt.title("Visibility Matrix")
    plt.savefig("Plots/Matrix.svg")

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
    c = scipy.constants.c
    lam = c/f
    return lam

# This code is taken from the fundamentals of interferometry notebook.
def UVellipse(u,v,w,a,b,v0):
    fig=plt.figure(0, figsize=(8,8))

    e1=Ellipse(xy=np.array([0,v0]),width=2*a,height=2*b,angle=0)
    e2=Ellipse(xy=np.array([0,-v0]),width=2*a,height=2*b,angle=0)

    ax=fig.add_subplot(111,aspect="equal")

    ax.plot([0],[v0],"go")
    ax.plot([0],[-v0],"go")
    ax.plot(u[0],v[0],"bo")
    ax.plot(u[-1],v[-1],"bo")

    ax.plot(-u[0],-v[0],"ro")
    ax.plot(-u[-1],-v[-1],"ro")

    ax.add_artist(e1)
    e1.set_lw(1)
    e1.set_ls("--")
    e1.set_facecolor("w")
    e1.set_edgecolor("b")
    e1.set_alpha(0.5)
    ax.add_artist(e2)

    e2.set_lw(1)
    e2.set_ls("--")
    e2.set_facecolor("w")
    e2.set_edgecolor("r")
    e2.set_alpha(0.5)
    ax.plot(u,v,"b")
    ax.plot(-u,-v,"r")
    ax.grid(True)
    plt.savefig('Plots/UVCoverage.png')

def plot_baseline(b_ENU, L, f, ant1, ant2):
    B = get_B(b_ENU, L)

    lam = get_lambda(f)
    # h = np.linspace(-12,12,num=time_steps)*np.pi/12
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
    UVellipse(u,v,w,a,b,v0)

def plot_array(antennas, name):
    plt.figure()
    plt.scatter(antennas[:,0], antennas[:,1])
    plt.grid(True)
    plt.xlabel('E-W [m]')
    plt.ylabel('N-S [m]')
    plt.title(name + ' Array Layout')
    plt.savefig('Plots/AntennaLayout.png')

def plot_visibilities(u, v, b_ENU, L, f):

    #################################################################################################
    RA_sources = np.array([-4 + 44/60.0 + 6.686/3600, -4 + 44/60.0 + 6.686/3600])
    DEC_sources = np.array([-74 + 39/60.0 + 37.481/3600, -73 + 39.0/60.0 + 37.298/3600])
    Flux_sources_labels = np.array(["1 Jy","0.2 Jy"])
    Flux_sources = np.array([1,0.2]) #in Jy
    step_size = 200
    print("Phase center     Source 1")
    print("RA="+str(RA_sources))
    print("DEC="+str(DEC_sources))
    RA_rad = np.array(RA_sources)*(np.pi/12)
    DEC_rad = np.array(DEC_sources)*(np.pi/180)
    RA_delta_rad = RA_rad-RA_rad[0]

    l = np.cos(DEC_rad)*np.sin(RA_delta_rad)
    m = (np.sin(DEC_rad)*np.cos(DEC_rad[0])-np.cos(DEC_rad)*np.sin(DEC_rad[0])*np.cos(RA_delta_rad))
    print("l=",l)
    print("m=",m)

    point_sources = np.zeros((len(RA_sources),3))
    point_sources[:,0] = Flux_sources
    point_sources[:,1] = l[0:]
    point_sources[:,2] = m[0:]
    dec = DEC_sources[0]
    #################################################################################################

    B = get_B(b_ENU, L)
    lam = get_lambda(f)

    X = B[0]
    Y = B[1]
    Z = B[2]
    u_d = lam**(-1)*(np.sin(h)*X+np.cos(h)*Y)
    v_d = lam**(-1)*(-np.sin(dec)*np.cos(h)*X+np.sin(dec)*np.sin(h)*Y+np.cos(dec)*Z)
    # delta = 60*(np.pi/180) #Declination in degrees
    u = np.linspace(-1*(np.amax(np.abs(u_d)))-10, np.amax(np.abs(u_d))+10, num=step_size, endpoint=True)
    v = np.linspace(-1*(np.amax(abs(v_d)))-10, np.amax(abs(v_d))+10, num=step_size, endpoint=True)
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
    plt.savefig('Plots/Visibilities.png')
    plt.close()

    u_track = u_d
    v_track = v_d
    z = np.zeros(u_track.shape).astype(complex)
    plt.figure()

    s = point_sources.shape
    for counter in range(1, s[0]+1):
        A_i = point_sources[counter-1,0]
        l_i = point_sources[counter-1,1]
        m_i = point_sources[counter-1,2]
        z += A_i*np.exp(-1*2*np.pi*1j*(u_track*l_i+v_track*m_i))
    plt.subplot(121)
    plt.plot(z.real)
    plt.xlabel("Timeslots")
    plt.ylabel("Jy")
    plt.title("Real: sampled visibilities")

    plt.subplot(122)
    plt.plot(z.imag)
    plt.xlabel("Timeslots")
    plt.ylabel("Jy")
    plt.title("Imag: sampled visibilities")
    plt.savefig('Plots/SampledVisibilities.png')
    plt.close()
