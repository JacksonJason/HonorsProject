import TARTRequests as TR
import numpy as np
import Utilities as ut

def get_antenna_layout():
    layoutJSON = TR.antenna_layout()
    return np.array(layoutJSON)

if __name__ == "__main__":
    layout = get_antenna_layout()
    # ut.plot_array(layout)
    b12 = layout[1] - layout[2]
    L,f = TR.get_latitude_and_frequency()
    # ut.plot_baseline(b12, L, f, layout[0], layout[1])
    print(TR.get_visibilities())
