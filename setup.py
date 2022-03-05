##########################################
##########################################
####
####  Set up data/ folder to contain all data sources for hourly eGrid calculation.
####  Data sources:
####     - PUDL (CEMS, EIA 860, EIA 923)
####     - eGRID
####  EIA-930 data is used, but is availible as needed via API
####
####  Check https://zenodo.org/record/5701406#.Yh-3kN9OmLc for new PUDL releases
####
####  These sources can also be downloaded manually from the listed links.
####  If downloaded manually, the expected structure is:
####     hourly-egrid/data/pudl-v0.5.0-2021-11-14 contains unzipped pudl folder from zenodo
####     hourly-egrid/data/egrid contains egrid excel spreadsheets
####     hourly-egrid/data/eia930/  empty dir to hold cached 930 data
####
##########################################
##########################################

import os
import requests
import tarfile

if os.getcwd().split("/")[-1] != "hourly-egrid":
    print("You are not running this script in the top-level directory of the hourly egrid project.")
    throw(AssertionError)

os.mkdir("data")
os.chdir("data")

############### EIA-930 data ########################

# This will be downloaded and cached as-needed, so for now just make directory
os.mkdir("eia930")


############### PUDL data ###########################

print("Downloading pudl (be patient, this is a large file).....")
r = requests.get("https://zenodo.org/record/5701406/files/pudl-v0.5.0-2021-11-14.tgz", params={"download":"1"}, stream=True)
with open("pudl.tgz", 'wb') as fd:
    for chunk in r.iter_content(chunk_size=1024):
        fd.write(chunk)

with tarfile.open("pudl.tgz") as tar:
    tar.extractall()

################# eGRID data #########################

egrids = ["2021-02/egrid2019_data.xlsx", "2021-02/egrid2018_data_v2.xlsx"]
os.mkdir("egrid")
os.chdir("egrid")

print("Downloading eGRID....")
for filename in egrids:
    r = requests.get("https://www.epa.gov/sites/default/files/"+filename, stream=True)
    localname = filename.split("/")[-1]
    with open(localname, 'wb') as fd:
        for chunk in r.iter_content(chunk_size=1024):
            fd.write(chunk)

os.chdir("../..")
