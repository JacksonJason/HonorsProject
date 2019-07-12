import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.constants
import math
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Ellipse
import matplotlib.pylab as pl
import Tigger

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
    plt.savefig('Plots/Antenna_Visibilities.png', transparent=True)
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
    plt.savefig("Plots/Matrix.svg", transparent=True)
    plt.close()

def get_B(b_ENU, L):
    D = math.sqrt(np.sum((b_ENU)**2))
    A = np.arctan2(b_ENU[0],b_ENU[1])
    E = np.arcsin(b_ENU[2]/D)
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
    plt.title("UV Coverage")
    plt.savefig('Plots/UVCoverage.png', transparent=True)
    plt.close()

def plot_baseline(b_ENU, L, f, h0, h1, dec):
    # dec = model.dec0
    B = get_B(b_ENU, L)

    lam = get_lambda(f)
    h = np.linspace(h0,h1,num=600)*np.pi/12
    X = B[0]
    Y = B[1]
    Z = B[2]
    u = lam**(-1)*(np.sin(h)*X+np.cos(h)*Y)
    v = lam**(-1)*(-np.sin(dec)*np.cos(h)*X+np.sin(dec)*np.sin(h)*Y+np.cos(dec)*Z)
    w = lam**(-1)*(np.cos(dec)*np.cos(h)*X-np.cos(dec)*np.sin(h)*Y+np.sin(dec)*Z)
    # plotBL.UV(u,v,w)
    a = np.sqrt(X**2+Y**2)/lam # major axis
    b = a*np.sin(dec)              # minor axis
    v0 = (Z/lam)*np.cos(dec)  # center of ellipse
    UVellipse(u,v,w,a,b,v0)

def plot_array(antennas, name):
    plt.figure()
    plt.scatter(antennas[:,0], antennas[:,1])
    plt.grid(True)
    plt.xlabel('E-W [m]')
    plt.ylabel('N-S [m]')
    plt.title(name + ' Array Layout')
    plt.savefig('Plots/' + name + 'AntennaLayout.png', transparent=True)
    plt.close()

def plot_visibilities(b_ENU, L, f, h0, h1, model_name, cos):
    h = np.linspace(h0,h1,num=600)*np.pi/12
    model = Tigger.load(model_name)
############################################################
    RA_sources = []
    DEC_sources = []
    # Flux_sources_labels = []
    Flux_sources = []
    for val in model.sources:
        RA_sources.append(val.pos.ra)
        DEC_sources.append(val.pos.dec)
        # Flux_sources_labels.append(str(val.flux.I))
        Flux_sources.append(val.flux.I)
    RA_sources = np.array(RA_sources)
    DEC_sources = np.array(DEC_sources)
    # Flux_sources_labels = np.array(Flux_sources_labels)
    Flux_sources = np.array(Flux_sources)

    ra_0 = model.ra0
    dec_0 = model.dec0
    ra_0_rad = ra_0 * (np.pi/12)
    dec_0_rad = dec_0 * (np.pi/180)

    step_size = 200
    RA_rad = RA_sources*(np.pi/12)
    DEC_rad = DEC_sources*(np.pi/180)
    RA_delta_rad = RA_rad-ra_0_rad


    l = np.cos(DEC_rad)*np.sin(RA_delta_rad)
    m = (np.sin(DEC_rad)*np.cos(dec_0_rad)-np.cos(DEC_rad)*np.sin(dec_0_rad)*np.cos(RA_delta_rad))

    point_sources = np.zeros((len(RA_sources),3))
    point_sources[:,0] = Flux_sources
    point_sources[:,1] = l[0:]
    point_sources[:,2] = m[0:]
    dec = dec_0
    #################################################################################################


    B = get_B(b_ENU, L)
    lam = get_lambda(f)

    X = B[0]
    Y = B[1]
    Z = B[2]
    u_d = lam**(-1)*(np.sin(h)*X+np.cos(h)*Y)
    v_d = lam**(-1)*(-np.sin(dec)*np.cos(h)*X+np.sin(dec)*np.sin(h)*Y+np.cos(dec)*Z)
    u = np.linspace(-1*(np.amax(np.abs(u_d)))-10, np.amax(np.abs(u_d))+10, num=step_size, endpoint=True)
    v = np.linspace(-1*(np.amax(abs(v_d)))-10, np.amax(abs(v_d))+10, num=step_size, endpoint=True)
    uu, vv = np.meshgrid(u, v)
    zz = np.zeros(uu.shape).astype(complex)
    s = point_sources.shape
    for counter in range(0, s[0]):
        A_i = point_sources[counter,0]
        l_i = point_sources[counter,1]
        m_i = point_sources[counter,2]
        zz += A_i*np.exp(-2*np.pi*1j*(uu*l_i+vv*m_i))
    zz = zz[:,::-1]

    plt.figure()
    plt.set_cmap('viridis')
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
    plt.title("Imaginary part of visibilities")
    plt.savefig('Plots/Visibilities.png', transparent=True)
    plt.close()

    uv_tracks = plot_sampled_visibilities(point_sources, u_d, v_d)
    uv = []
    for i in range(len(u_d)):
        uv.append([u_d[i], v_d[i]])
    uv = np.array(uv)
    cell_size_l = 0.1 #should be manually input
    cell_size_m = 0.1
    degrees_l = 4 #should determine through sky model
    degrees_m = 4
    Nl = int(np.round(degrees_l / cell_size_l))
    Nm = int(np.round(degrees_m / cell_size_m))
    Nl = find_closest_power_of_two(Nl * 5)
    Nm = find_closest_power_of_two(Nm * 5)
    rad_d_l = cell_size_l * (np.pi/180)
    rad_d_m = cell_size_m * (np.pi/180)
    cell_size_u = 1 / (2 * Nl * rad_d_l)
    cell_size_v = 1 / (2 * Nm * rad_d_m)

    scaled_uv = np.copy(uv)
    scaled_uv[:,0] *= np.deg2rad(cell_size_l * Nl)
    scaled_uv[:,1] *= np.deg2rad(cell_size_m * Nm)

    gridded = grid(Nl, Nm, uv_tracks, cell_size_u, cell_size_v, scaled_uv)

    # plot_sky_model(l, m, Flux_sources, Nl, Nm)
    if cos == "1":
        plot_sky_model(l*(180/np.pi), m*(180/np.pi), Flux_sources, Nl, Nm, "l [degrees]", "m [degrees]")

        print(RA_delta_rad)
        L = np.cos(dec_0) * np.sin(0)
        M = np.sin(dec_0) * np.cos(dec_0) - np.cos(dec_0) * np.sin(dec_0) * np.cos(0)
        image_visibilities(gridded, Nl, Nm, cell_size_l, cell_size_m, L, M)
    else:
        plot_sky_model(RA_sources, DEC_sources, Flux_sources, Nl, Nm, "RA", "DEC")
        # convert grid to RA/DEC
        # dec = math.asin(m * np.cos(dec_0) + np.sin(dec_0) * math.sqrt(1 - l**2 - m**2))
        # ra = ra_0 + np.atan(1 / (np.cos(dec_0) * math.sqrt(1 - l**2 - m**2) - m * np.sin(dec_0)))
        image_visibilities(gridded, Nl, Nm, cell_size_l, cell_size_m, dec_0, ra_0)


def find_closest_power_of_two(number):
    s = 2
    for i in range(0, 15):
        if number < s:
            return s
        else:
            s *= 2

def grid(Nl, Nm, uv_tracks, d_u, d_v, uv):
    vis = np.zeros((Nl, Nm), dtype=complex)
    counter = np.zeros((Nl, Nm))
    positions = np.empty((Nl, Nm), dtype=object)
    half_l = int(Nl / 2)
    half_m = int(Nm / 2)
    o_d_u = d_u
    o_d_v = d_v
    for i in range(len(positions)):
        for j in range(len(positions[i])):
            positions[i][j] = (i - half_l, j - half_m)
    for i in range(len(uv)):
        y,x = int(np.round(uv[i][0])), int(np.round(uv[i][1])) # why is it flipped, a question to ask dr Grobler perhaps
        vis[x][y] += uv_tracks[i]
        counter[x][y] += 1
    for i in range(len(vis)):
        for j in range(len(vis[i])):
            if not counter[i][j] == 0:
                vis[i][j] = vis[i][j] / counter[i][j]
    return vis

def plot_sampled_visibilities(point_sources, u_d, v_d):
    u_track = u_d
    v_track = v_d
    z = np.zeros(u_track.shape).astype(complex)
    plt.figure()

    s = point_sources.shape
    for counter in range(0, s[0]):
        A_i = point_sources[counter,0]
        l_i = point_sources[counter,1]
        m_i = point_sources[counter,2]
        z += A_i*np.exp(-1*2*np.pi*1j*((u_track*l_i)+(v_track*m_i)))
    plt.subplot(121)
    plt.plot(z.real)
    plt.xlabel("Timeslots")
    plt.ylabel("Jy")
    plt.title("Real: sampled visibilities")

    plt.subplot(122)
    plt.plot(z.imag)
    plt.xlabel("Timeslots")
    plt.title("Imag: sampled visibilities")
    plt.savefig('Plots/SampledVisibilities.png', transparent=True)
    plt.close()
    return z

def plot_sky_model(l, m, Flux_sources, Nl, Nm, x, y):
    # Plot sky model in L and M
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)
    plt.xlabel(x)
    plt.ylabel(y)
    max_flux = max(Flux_sources)
    if max_flux > 1:
        col = (Flux_sources/max_flux)
    else:
        col = Flux_sources
    colour = []
    for i in col:
        colour.append((i,i,i))
    # l_d = (l*(180/np.pi)) * Nl/3
    # m_d = (m*(180/np.pi)) * Nm/3
    # im = np.zeros((Nl,Nm))
    # for i in range(len(l_d)):
    #     im[int(np.round(l_d[i]))][int(np.round(m_d[i]))] = Flux_sources[i]
    #     print(im[int(np.round(l_d[i]))][int(np.round(m_d[i]))])
    # plt.imshow(im, cmap="gray",)
    plt.scatter(l,m,c=colour,s=8)
    # plt.scatter(l,m,c=colour,s=8)
    ax.set_facecolor('xkcd:black')
    fig.patch.set_alpha(0)
    plt.title("Sky Model")
    plt.savefig("Plots/SkyModel.png", transparent=False)
    plt.close()

def image_visibilities(grid, Nl, Nm, cell_size_l, cell_size_m, RA, DECLINATION):
    image = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(grid)))
    image = np.abs(image)
    img = plt.figure(figsize=(10,10))
    plt.title("Reconstructed Sky Model")
    plt.set_cmap('gray')
    plt.imshow(np.real(image), origin='lower', extent=[RA - Nl / 2 * cell_size_l, RA + Nl / 2 * cell_size_l,
                                                            DECLINATION - Nm / 2 * cell_size_m, DECLINATION + Nm / 2 * cell_size_m])
    # plt.imshow(np.real(image), origin='lower')
    plt.savefig('Plots/ReconstructedSkyModel.png', transparent=True)
    plt.close()
