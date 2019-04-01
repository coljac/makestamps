import requests
import shutil
from requests.auth import HTTPBasicAuth

username = "cjacobs"
# password = "cja70chips"
password = "qEVgduvtARj"
prefix = "https://desar2.cosmology.illinois.edu/DESFiles/desarchive/"

def download_file(url, user=username, password=password):
    url = prefix + url
    local_filename = url.split('/')[-1]
    r = requests.get(url, stream=True, auth=HTTPBasicAuth(user, password))
    with open(local_filename, 'wb') as f:
        shutil.copyfileobj(r.raw, f)

    # print "Getting " + url
    return local_filename

# url = "https://desar2.cosmology.illinois.edu/DESFiles/desarchive/OPS/" +\
        # "multiepoch/Y3A1/r2689/DES0530-5248/p01/coadd/DES0530-5248_r2689p01_g.fits.fz"
if __name__ == "__main__":
    import sys
    fileurl = sys.argv[1]
    if fileurl.startswith(prefix):
        fileurl = fileurl[len(prefix):]
    else:
        pass
    download_file(fileurl)
