import TARTRequests as TR
import numpy as np
import Utilities as ut
import cherrypy
import os
import json
from jinja2 import Environment, FileSystemLoader
env = Environment(loader=FileSystemLoader('html'))

class pipeline(object):

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
        # np.savetxt("Plots/visibilities.csv", vis_matrix, delimiter=",")
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
    def generate_custom_graphs(self, input_file=None, lsm_file=None, baseline=None, cos=None):
        upload_path = os.path.dirname(__file__)
        if input_file is not "" and lsm_file is not "" and baseline is not "" and cos is not "":
            bl = baseline.split(" ")
            bl_1 = int(bl[0]) - 1
            bl_2 = int(bl[1]) - 1
            input_file = input_file.split("\\")[-1]
            lsm_file = lsm_file.split("\\")[-1]
            input_file = os.path.normpath(
            os.path.join(upload_path, input_file))
            with open("Antenna_Layouts/" + input_file) as outfile:
                json_antenna = json.load(outfile)
            custom_layout = np.array(json_antenna['antennas'])
            ut.plot_array(custom_layout, "Custom")
            b12 = custom_layout[bl_2] - custom_layout[bl_1]
            custom_L = json_antenna['latitude']
            custom_L = (np.pi/180)* (custom_L[0] + custom_L[1]/60. + custom_L[2]/3600.)
            custom_f = json_antenna['frequency']
            custom_f = custom_f * 10**9
            sha = json_antenna['sha']
            eha = json_antenna['eha']
            dec = json_antenna['center_dec']
            dec = dec[0] + dec[1]/60. + dec[2]/3600.
            # asc = json_antenna['center_asc']
            ut.plot_baseline(b12, custom_L, custom_f, sha, eha, dec)
            ut.plot_visibilities(b12, custom_L, custom_f, sha, eha, "Sky_Models/" + lsm_file, cos)

    @cherrypy.expose
    def generate_graphs(self):
        layout = self.get_antenna_layout()
        L,f = TR.get_latitude_and_frequency()
        visibilities = self.make_vis_matrix()
        ut.plot_array(layout, "TART")
        # b12 = layout[1] - layout[2]
        # ut.plot_baseline(b12, L, f, layout[0], layout[1])
        # ut.plot_visibilities(0,0, b12, L, f)

    @cherrypy.expose
    def index(self):
        tmpl = env.get_template('index.html')
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
