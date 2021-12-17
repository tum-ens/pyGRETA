# How to create authentication file for downloading data with GES DISC account (https://disc.gsfc.nasa.gov/data-access)
# On Mac/Linux:
# 1) cd ~ or cd $HOME
# 2) touch .netrc
# 3) echo "machine urs.earthdata.nasa.gov login <uid> password <password>" >> .netrc (where <uid> is your user name and <password> is your Earthdata Login password without the brackets)
# 4) chmod 0600 .netrc (so only you can access it)
# On Windows:
# 1) Open Notepad
# 2) Enter (without quotes): machine urs.earthdata.nasa.gov login <uid> password <password>
# 3) Save as: C:\Users\<username>\.netrc

# Generate link file
# Go to GES disc data collection and search for "tavg1_2d_slv_NX"

import os
import requests

# Parameters
file_links = "subset_M2I1NXLFO_5.12.4_20210930_073038.txt"  # Downloaded links from GES DISC
file_links = 'subset_MAT1NXSLV_5.2.0_20211217_094609.txt'
path_downloadedFiles = "download/"

# Create folder for downloaded files
if not os.path.exists(path_downloadedFiles):
    os.makedirs(path_downloadedFiles)

# Download each link of text file
file = open(file_links, 'r')
Lines = file.readlines()

for line in Lines:
    url = line.strip()
    filename = url.split("/")[-1]

    if not os.path.isfile(path_downloadedFiles + filename):
        result = requests.get(url)
        try:
            result.raise_for_status()
            f = open(path_downloadedFiles + filename, 'wb')
            f.write(result.content)
            f.close()
            print('contents of URL written to ' + filename)
        except:
            print('requests.get() returned an error code ' + str(result.status_code))