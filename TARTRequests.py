import requests

def antenna_layout():
    r = requests.get("http://146.232.222.105/api/v1/imaging/antenna_positions")
    if not r.status_code == 200:
        return "Error retreiving antenna positions"
    return r.json()

def get_visibilities():
    r = requests.get("http://146.232.222.105/api/v1/imaging/vis")
    if not r.status_code == 200:
        return "Error retreiving visibilities"
    return r.json()
