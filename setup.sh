##########################################
##########################################
####
####  Set up data/ folder to contain all data sources for hourly eGrid calculation.
####  Data sources:
####     - PUDL (CEMS, EIA 860, EIA 923)
####     - eGRID
####     - EIA 930 <- can access via API on as-needed basis.
####
####  These sources can also be downloaded manually from the listed links.
####
####
##########################################
##########################################

mkdir data
cd data
wget https://sandbox.zenodo.org/record/764696/files/databeta-2021-11-14.tgz?download=1
tar -xzf databeta-2021-11-14.tgz

mkdir egrid
cd egrid
wget https://www.epa.gov/sites/default/files/2021-02/egrid2019_data.xlsx
wget https://www.epa.gov/sites/default/files/2020-03/egrid2018_data_v2.xlsx
cd ..
