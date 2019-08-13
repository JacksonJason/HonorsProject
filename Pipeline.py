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
    def generate_custom_graphs(self, input_file=None, lsm_file=None, baseline=None, cos=None, cell_size=None, res=None):
        upload_path = os.path.dirname(__file__)
        if res is not "" and input_file is not "" and lsm_file is not "" and baseline is not "" and cos is not "" and cell_size is not "":
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
            b = custom_layout[bl_2] - custom_layout[bl_1]
            custom_L = json_antenna['latitude']
            custom_L = (np.pi/180)* (custom_L[0] + custom_L[1]/60. + custom_L[2]/3600.)
            custom_f = json_antenna['frequency']
            custom_f = custom_f * 10**9
            sha = json_antenna['sha']
            eha = json_antenna['eha']
            dec = json_antenna['center_dec']
            dec = dec[0] + dec[1]/60. + dec[2]/3600.
            ut.plot_baseline(b, custom_L, custom_f, sha, eha, dec, "CUSTOM")
            uv, uv_tracks, dec_0 = ut.plot_visibilities(b, custom_L, custom_f, sha, eha, "Sky_Models/" + lsm_file, cos, custom_layout)
            ut.image(uv, uv_tracks, cell_size, cos, dec_0, res, "CUSTOM")

    @cherrypy.expose
    def generate_graphs(self, cos=None, cell_size=None, res=None):
        if res is not "" and cos is not "" and cell_size is not "":
            layout = self.get_antenna_layout()
            L,f = TR.get_latitude_and_frequency()
            visibilities = self.make_vis_matrix()
            ut.plot_array(layout, "TART")
            #dec_0 is 0
            # h = np.linspace(h0,h1,num=600)*np.pi/12
            all_uv = []
            all_uv_tracks = []
            for i in range(len(layout)):
                for j in range(i+1, len(layout)):
                    b = layout[j] - layout[i]
                    u_d, v_d =  ut.get_uv_tracks(b, L, f, 0, L)
                    uv = []
                    uv.append([[u_d, v_d]])
                    uv = np.array(uv)
                    all_uv.append(uv)
                    uv_tracks = [[1]]
                    all_uv_tracks.append(uv_tracks)
                    b = layout[i] - layout[j]
                    u_d, v_d =  ut.get_uv_tracks(b, L, f, 0, L)
                    uv = []
                    uv.append([[u_d, v_d]])
                    uv = np.array(uv)
                    all_uv.append(uv)
                    uv_tracks = [[1]]
                    all_uv_tracks.append(uv_tracks)
            ut.image(uv, uv_tracks, cell_size, cos, 0, res, "TART")



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
