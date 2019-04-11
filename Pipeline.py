import TARTRequests as TR
import numpy as np
import Utilities as ut
import cherrypy
import os
from jinja2 import Environment, FileSystemLoader
env = Environment(loader=FileSystemLoader('html'))

class pipeline(object):
    global layout
    global L, f

    def get_antenna_layout(self):
        layoutJSON = TR.antenna_layout()
        return np.array(layoutJSON)

    def make_vis_matrix(self):
        vis = TR.get_visibilities()
        vis = np.array(vis["data"])
        i,j = self.parse_vis(vis)
        vis_matrix = np.zeros((i+1, j+1)).astype(complex)
        for v in vis:
            vis_matrix[v["i"]][v["j"]] = v["re"] + 1j*v["im"]
        cc = vis_matrix.T
        cc = np.conj(cc)
        for i in range(cc.shape[0]):
            cc[i][i] = 0 + 0j
        vis_matrix = vis_matrix + cc
        ut.draw_matrix(vis_matrix)
        # ut.tabulate_matrix(vis_matrix)
        np.savetxt("Plots/visibilities.csv", vis_matrix, delimiter=",")
        return vis_matrix

    def parse_vis(self, vis):
        i = 0
        j = 0
        for v in vis:
            if v['i'] >= i:
                i = v['i']
            if v['j'] >= j:
                j = v['j']
        return i+1,j

    @cherrypy.expose
    def generate_graphs(self):
        visibilities = self.make_vis_matrix()
        ut.plot_array(layout, "TART")
        b12 = layout[1] - layout[2]
        ## Need Declination before I can continue propely. Need to figure out how to extract that from the
        ## JSON data.
        ut.plot_baseline(b12, L, f, layout[0], layout[1])
        ut.plot_visibilities(0,0, b12, L, f)

    @cherrypy.expose
    def index(self):
        tmpl = env.get_template('index.html')
        global layout
        layout = self.get_antenna_layout()
        global L, f
        L,f = TR.get_latitude_and_frequency()
        return tmpl.render(target='Imaging pipeline')

if __name__ == '__main__':
    conf = {
        '/': {
            'tools.sessions.on': True,
            'tools.staticdir.root': os.path.abspath(os.getcwd())
        },
        '/static': {
            'tools.staticdir.on': True,
            'tools.staticdir.dir': './Plots'
        }
    }
    cherrypy.quickstart(pipeline(), '/', conf)
