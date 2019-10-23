import TARTRequests as TR
import numpy as np
import Utilities as ut
import cherrypy
import os
import json
from jinja2 import Environment, FileSystemLoader
env = Environment(loader=FileSystemLoader('html'))

class pipeline(object):

    def make_vis_matrix(self, loc):
        """
        Requests the visibilities from the API, then converts that into a 2d
        array of visibilities.

        :param loc: The location of the interferometer, NZ or SA
        :type loc: str

        :returns: A matrix of visibilities
        """
        vis = TR.get_visibilities(loc)
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

        return vis_matrix

    def parse_vis(self, vis):
        """
        Finds the maximum i and j postions in the vis array

        :param vis: The visibility to convert
        :type vis: numpy array

        :returns: the x and y max coordinates
        """
        i = 0
        j = 0

        for v in vis:
            if v['i'] >= i:
                i = v['i']

            if v['j'] >= j:
                j = v['j']

        return i+1,j

    def get_antenna_layout(self, loc):
        """
        Gets the antenna layout and returns it in JSON format.

        :returns: The antenna layout in JSON format.
        """
        layoutJSON = TR.antenna_layout(loc)
        return np.array(layoutJSON)

    @cherrypy.expose
    def generate_custom_graphs(self, input_file=None, lsm_file=None,
                                baseline=None, cell_size=None, res=None,
                                showGrid=False):
        """
        Called by the HTML page when the generate Custom graphs button is pressed.
        The HTML page sends the information entered in the fields to the function and it is
        executed.

        :param input_file: The file that contains the telescope information
        :type input_file: str

        :param lsm_file: The file that contains the sky model
        :type lsm_file: str

        :param baseline: The baseline to plot, for information purposes, in the
                         form "x y" where x and y are numbers
        :type baseline: str

        :param cell_size: The cell size for the imaging process
        :type cell_size: str

        :param res: The resolution for the imaging process
        :type: str

        :param showGrid: Whether or not to show the grid on the image.
        :type: str

        :returns: nothing
        """
        upload_path = os.path.dirname(__file__)
        if (res is not "" and input_file is not "" and lsm_file is not ""
            and baseline is not "" and cell_size is not ""):

            bl = baseline.split(" ")
            bl_1 = int(bl[0]) - 1
            bl_2 = int(bl[1]) - 1
            input_file = input_file.split("\\")[-1]
            lsm_file = lsm_file.split("\\")[-1]
            input_file = os.path.normpath(os.path.join(upload_path, input_file))

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
            uv, uv_tracks, dec_0, ra_0 = ut.get_visibilities(b, custom_L, custom_f, sha, eha, "Sky_Models/" + lsm_file, custom_layout)
            if showGrid == "true":
                showGrid = True
            else:
                showGrid = False
            ut.image(uv, uv_tracks, cell_size, dec_0, res, "CUSTOM", showGrid, ra_0)

    @cherrypy.expose
    def generate_graphs(self, cell_size=None, loc=None, showGrid=False):
        """
        Called by the HTML page when the generate Tart graphs button is pressed.
        The HTML page sends the information entered in the fields to the function and it is
        executed.

        :param cell_size: The cell size for the imaging process
        :type cell_size: str

        :param loc: location, which telescope to use, New Zealand or South Africa
        :type: str

        :param showGrid: Whether or not to show the grid on the image.
        :type: str

        :returns: nothing
        """
        if cell_size is not "" and loc is not "":
            layout = self.get_antenna_layout(loc)
            L,f = TR.get_latitude_and_frequency(loc)
            visibilities = self.make_vis_matrix(loc)
            ut.plot_array(layout, "TART")

            # dec_0 is 0 as it is TART
            all_uv = []
            all_uv_tracks = []

            for i in range(len(layout)):
                for j in range(i+1, len(layout)):
                    b = layout[j] - layout[i]
                    u_d, v_d =  ut.get_uv_tracks(b, L, f, 0, L)
                    uv = []
                    uv.append([u_d, v_d])
                    uv = np.array(uv)
                    all_uv.append(uv)
                    uv_tracks = [visibilities[i][j]]
                    all_uv_tracks.append(uv_tracks)

                    b = layout[i] - layout[j]
                    u_d, v_d =  ut.get_uv_tracks(b, L, f, 0, L)
                    uv = []
                    uv.append([u_d, v_d])
                    uv = np.array(uv)
                    all_uv.append(uv)
                    uv_tracks = [visibilities[j][i]]
                    all_uv_tracks.append(uv_tracks)

            res = 2 * 180/np.pi
            if showGrid == "true":
                showGrid = True
            else:
                showGrid = False
            ut.image(all_uv, all_uv_tracks, cell_size, 0, res, "TART", showGrid)


    @cherrypy.expose
    def index(self):
        """
        Tells Cherrypy where to start
        :returns: the rendered template from Jinga
        """
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
        } # This sets the locations of images and such for cherrypy
    }
    cherrypy.quickstart(pipeline(), '/', conf)
