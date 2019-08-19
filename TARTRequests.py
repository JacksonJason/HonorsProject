import requests

def antenna_layout(loc):
    if loc == "1":
        # SU
        r = requests.get("http://146.232.222.105/api/v1/imaging/antenna_positions")
    else:
        # NZ
        r = requests.get("https://tart.elec.ac.nz/signal/api/v1/imaging/antenna_positions")
    if not r.status_code == 200:
        return "Error retreiving antenna positions"
    return r.json()

def get_visibilities(loc):
    if loc == "1":
        # SU
        r = requests.get("http://146.232.222.105/api/v1/imaging/vis")
    else:
        # NZ
        r = requests.get("https://tart.elec.ac.nz/signal/api/v1/imaging/vis")

    if not r.status_code == 200:
        return "Error retreiving visibilities"
    return r.json()

def get_latitude_and_frequency(loc):
    if loc == "1":
        # SU
        r = requests.get("http://146.232.222.105/api/v1/info")
    else:
        # NZ
        r = requests.get("https://tart.elec.ac.nz/signal/api/v1/info")
    if not r.status_code == 200:
        return "Error retreiving latitude"
    info = r.json()["info"]
    return(info["location"]["lat"], info["operating_frequency"])
