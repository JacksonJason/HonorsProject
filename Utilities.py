import numpy as np
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants
import math
from mpl_toolkits.mplot3d import Axes3D
import plotBL
from matplotlib.patches import Ellipse

def plot_baseline():
    pass

def plot_array(antennas):
    plt.scatter(antennas[:,0], antennas[:,1])
    plt.grid(True)
    plt.xlabel('E-W [m]')
    plt.ylabel('N-S [m]')
    plt.title('TART Array Layout')
    plt.show()
