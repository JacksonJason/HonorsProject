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
    plt.set_cmap('viridis')
    plt.imshow(matrix.real)
    plt.title("Real: visibilities")

    plt.subplot(122)
    plt.imshow(matrix.imag)
    plt.title("Imag: visibilities")
    plt.savefig('Plots/Antenna_Visibilities.png', transparent=True)
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
def UVellipse(u,v,w,a,b,v0, name):
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
    plt.title("UV Coverage", size=20)
    plt.xlabel("u", size=18)
    plt.ylabel("v", size=18)
    plt.savefig('Plots/' + name + 'UVCoverage.png', transparent=True)
    plt.close()

def plot_baseline(b_ENU, L, f, h0, h1, dec, name):
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
    UVellipse(u,v,w,a,b,v0, name)

def plot_array(antennas, name):
    plt.figure()
    plt.scatter(antennas[:,0], antennas[:,1])
    plt.grid(True)
    plt.xlabel('E-W [m]')
    plt.ylabel('N-S [m]')
    plt.title(name + ' Array Layout')
    plt.savefig('Plots/' + name + 'AntennaLayout.png', transparent=True)
    plt.close()

def get_uv_tracks(b_ENU, L, f, h, dec):
    B = get_B(b_ENU, L)
    lam = get_lambda(f)

    X = B[0]
    Y = B[1]
    Z = B[2]
    u_d = lam**(-1)*(np.sin(h)*X+np.cos(h)*Y)
    v_d = lam**(-1)*(-np.sin(dec)*np.cos(h)*X+np.sin(dec)*np.sin(h)*Y+np.cos(dec)*Z)
    return u_d, v_d

def get_uv_and_tracks(b_ENU, L, f, h, dec, point_sources):
    u_d, v_d = get_uv_tracks(b_ENU, L, f, h, dec)
    step_size = 200
    u = np.linspace(-1*(np.amax(np.abs(u_d)))-10, np.amax(np.abs(u_d))+10, num=step_size, endpoint=True)
    v = np.linspace(-1*(np.amax(abs(v_d)))-10, np.amax(abs(v_d))+10, num=step_size, endpoint=True)
    uu, vv = np.meshgrid(u, v)
    uv_tracks = plot_sampled_visibilities(point_sources, u_d, v_d)
    uv = []
    for i in range(len(u_d)):
        uv.append([u_d[i], v_d[i]])
    uv = np.array(uv)
    return uv, u_d, v_d, uu, vv, uv_tracks

def plot_visibilities(b_ENU, L, f, h0, h1, model_name, cos, layout):
    h = np.linspace(h0,h1,num=600)*np.pi/12
    model = Tigger.load(model_name)
############################################################
    RA_sources = []
    DEC_sources = []
    Flux_sources = []
    for val in model.sources:
        RA_sources.append(val.pos.ra)
        DEC_sources.append(val.pos.dec)
        Flux_sources.append(val.flux.I)
    RA_sources = np.array(RA_sources)
    DEC_sources = np.array(DEC_sources)
    Flux_sources = np.array(Flux_sources)

    ra_0 = model.ra0
    dec_0 = model.dec0
    ra_0_rad = ra_0 * (np.pi/12)
    dec_0_rad = dec_0 * (np.pi/180)

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

    if cos == "1":
        plot_sky_model(l*(180/np.pi), m*(180/np.pi), Flux_sources, "l [degrees]", "m [degrees]")
    else:
        plot_sky_model(RA_sources, DEC_sources, Flux_sources, "RA", "DEC")

    uv, u_d, v_d, uu, vv, uv_tracks = get_uv_and_tracks(b_ENU, L, f, h, dec, point_sources)

    plt.subplot(121)
    plt.plot(uv_tracks.real)
    plt.xlabel("Timeslots")
    plt.ylabel("Jy")
    plt.title("Real: sampled visibilities")

    plt.subplot(122)
    plt.plot(uv_tracks.imag)
    plt.xlabel("Timeslots")
    plt.title("Imag: sampled visibilities")
    plt.savefig('Plots/SampledVisibilities.png', transparent=True)
    plt.close()

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

    all_uv_tracks = []
    all_uv = []

    for i in range(len(layout)):
        for j in range(i+1, len(layout)):
            b = layout[j] - layout[i]
            uv, u_d, v_d, uu, vv, uv_tracks = get_uv_and_tracks(b, L, f, h, dec, point_sources)
            all_uv.append(uv)
            all_uv_tracks.append(uv_tracks)
            b = layout[i] - layout[j]
            uv, u_d, v_d, uu, vv, uv_tracks = get_uv_and_tracks(b, L, f, h, dec, point_sources)
            all_uv.append(uv)
            all_uv_tracks.append(uv_tracks)

    return all_uv, all_uv_tracks, dec_0


def image(uv, uv_tracks, cell_size, cos, dec_0, res, name):
    c_s = float(cell_size)
    cell_size_l = c_s
    cell_size_m = c_s
    degrees_l = float(res)
    degrees_m = float(res)
    Nl = int(np.round(degrees_l / cell_size_l))
    Nm = int(np.round(degrees_m / cell_size_m))
    Nl = find_closest_power_of_two(Nl)
    Nm = find_closest_power_of_two(Nm)
    rad_d_l = cell_size_l * (np.pi/180)
    rad_d_m = cell_size_m * (np.pi/180)
    cell_size_u = 1 / (2 * Nl * rad_d_l)
    cell_size_v = 1 / (2 * Nm * rad_d_m)

    gridded, cell_size_error = grid(Nl, Nm, uv_tracks, cell_size_u, cell_size_v, uv, cell_size_l, cell_size_m)
    img = plt.figure(figsize=(10,10))
    plt.title("Baseline Grid", size=20)
    plt.set_cmap('nipy_spectral')
    im = plt.imshow(np.real(np.abs(gridded)), origin='lower')
    plt.axis('off')
    # plt.xlabel("l", size=18)
    # plt.ylabel("m", size=18)
    plt.savefig('Plots/' + name + 'grid.png', transparent=True)
    if cos == "1":
        L = np.cos(dec_0) * np.sin(0)
        M = np.sin(dec_0) * np.cos(dec_0) - np.cos(dec_0) * np.sin(dec_0) * np.cos(0)
        image = image_visibilities(gridded)
        psf = np.ones ((np.array(uv_tracks).shape), dtype=complex)
        psf_grid, cell_size_error = grid(Nl, Nm, psf, cell_size_u, cell_size_v, uv, cell_size_l, cell_size_m)
        psf_image = image_visibilities(psf_grid)
        scale_factor = psf_image[int(psf_image.shape[0]/2)][int(psf_image.shape[1]/2)]
        image /= scale_factor
        psf_image /= scale_factor
        draw_image(image, Nl, Nm, cell_size_l, cell_size_m, L, M, name+"SkyModel", "l", "m", cell_size_error)
        draw_image(psf_image, Nl, Nm, cell_size_l, cell_size_m, L, M, name+"PSF", "l", "m", cell_size_error)

    else:
        # convert grid to RA/DEC
        # dec = math.asin(m * np.cos(dec_0) + np.sin(dec_0) * math.sqrt(1 - l**2 - m**2))
        # ra = ra_0 + np.atan(1 / (np.cos(dec_0) * math.sqrt(1 - l**2 - m**2) - m * np.sin(dec_0)))
        image_visibilities(gridded, Nl, Nm, cell_size_l, cell_size_m, dec_0, ra_0, "SkyModel")
        # psf = np.ones((Nl, Nm), dtype=complex)
        # image_visibilities(psf, Nl, Nm, cell_size_l, cell_size_m, L, M, "PSF")


def find_closest_power_of_two(number):
    s = 2
    for i in range(0, 15):
        if number < s:
            return s
        else:
            s *= 2

def grid(Nl, Nm, uv_tracks, d_u, d_v, uv, cell_size_l, cell_size_m):
    vis = np.zeros((Nl, Nm), dtype=complex)
    counter = np.zeros((Nl, Nm))
    half_l = int(Nl / 2)
    half_m = int(Nm / 2)
    cell_size_error = False
    for i in range(len(uv)):
        scaled_uv = np.copy(uv[i])
        scaled_uv[:,0] *= np.deg2rad(cell_size_l * Nl)
        scaled_uv[:,1] *= np.deg2rad(cell_size_m * Nm)
        for j in range(len(scaled_uv)):
            y,x = int(np.round(scaled_uv[j][0])), int(np.round(scaled_uv[j][1]))
            x += half_l
            y += half_m
            if not x >= vis.shape[0] and not y >= vis.shape[1] and not x < -vis.shape[0] and not y < -vis.shape[1]:
                vis[x][y] += uv_tracks[i][j]
                counter[x][y] += 1
            else:
                cell_size_error = True

    for i in range(len(vis)):
        for j in range(len(vis[i])):
            if not counter[i][j] == 0:
                vis[i][j] = vis[i][j] / counter[i][j]
    return vis, cell_size_error

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

    return z

def plot_sky_model(l, m, Flux_sources, x, y):
    # Plot sky model in L and M
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)
    plt.xlabel(x, size=18)
    plt.ylabel(y, size=18)
    max_flux = max(Flux_sources)
    if max_flux > 1:
        col = (Flux_sources/max_flux)
    else:
        col = Flux_sources
    colour = []
    for i in col:
        colour.append((i,i,i))
    plt.scatter(l,m,c=colour,s=8)
    ax.set_facecolor('xkcd:black')
    fig.patch.set_alpha(0)
    plt.title("Sky Model", size=20)
    plt.savefig("Plots/SkyModel.png", transparent=False)
    plt.close()

def image_visibilities(grid):
    image = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(grid)))
    image = np.abs(image)
    return np.real(image)


def draw_image(image, Nl, Nm, cell_size_l, cell_size_m, RA, DECLINATION, name, x_title, y_title, cell_size_error):
    img = plt.figure(figsize=(10,10))
    plt.title("Reconstructed" + name,size=20)
    plt.set_cmap('nipy_spectral')
    im_vis = plt.imshow(image, origin='lower', extent=[RA - Nl / 2 * cell_size_l, RA + Nl / 2 * cell_size_l,
                                                        DECLINATION - Nm / 2 * cell_size_m, DECLINATION + Nm / 2 * cell_size_m])
    cbr = img.colorbar(im_vis)
    cbr.set_label('Jy per Beam',size=18)
    plt.xlabel(x_title,size=18)
    plt.ylabel(y_title,size=18)
    if cell_size_error:
        txt = "INVALID CELL SIZE"
        plt.figtext(0.5, 0.1, txt, wrap=True, horizontalalignment='center', fontsize=25, color='red')
    plt.savefig('Plots/Reconstructed' + name + '.png', transparent=True)
    plt.close()
