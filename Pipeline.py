import TARTRequests as TR
import numpy as np
import Utilities as ut

def get_antenna_layout():
    layoutJSON = TR.antenna_layout()
    return np.array(layoutJSON)

if __name__ == "__main__":
    layout = get_antenna_layout()
    ut.plot_array(layout)
