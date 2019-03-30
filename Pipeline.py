import TARTRequests as TR
import numpy as np
import Utilities as ut

def get_antenna_layout():
    layoutJSON = TR.antenna_layout()
    return np.array(layoutJSON)

def make_vis_matrix():
    vis = TR.get_visibilities()
    vis = np.array(vis["data"])
    i,j = parse_vis(vis)
    vis_matrix = np.zeros((i, j)).astype(complex)
    for v in vis:
        vis_matrix[v["i"]][v["j"]-1] = v["re"] + 1j*v["im"]
    # ut.draw_matrix(vis_matrix)
    np.savetxt("visibilities.csv", vis_matrix, delimiter=",")
    return vis_matrix

def parse_vis(vis):
    i = 0
    j = 0
    for v in vis:
        if v['i'] >= i:
            i = v['i']
        if v['j'] >= j:
            j = v['j']
    return i+1,j

if __name__ == "__main__":
    layout = get_antenna_layout()
    # ut.plot_array(layout)
    b12 = layout[1] - layout[2]
    L,f = TR.get_latitude_and_frequency()
    ## Need Declination before I can continue propely. Need to figure out how to extract that from the
    ## JSON data.
    # ut.plot_baseline(b12, L, f, layout[0], layout[1])
    visibilities = make_vis_matrix()
    # ut.plot_visibilities(0,0, b12, L, f)
