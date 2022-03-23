import pandas as pd
import os

year = 2019
technology = "OpenFieldPV,Biomass,WindOn"
path = "../configs-list/"

ISO = pd.read_csv("../assumptions/configs_list.csv", sep=";")

# Remove all old files
for f in os.listdir(path):
    os.remove(os.path.join(path, f))

# Iterate over each country
for i in range(len(ISO)):
    country = ISO.iloc[i]
    f = open(path + str(i) + "_" + country['Country Name']+ ".txt", "w")
    f.write("regions:gadm40_" + str(country['Country Code']) + "_" + country['ID Level'][-1] +".shp" + "\n")
    f.write("region_name:" + country['Country Name'] + "\n")
    f.write("subregions_name:" + country['Country Name'] + "_level" + country['ID Level'][-1] + "\n")
    f.write("country_code:" + str(country['Country Code']) + "\n")
    f.write("year:" + str(year) + "\n")
    f.write("technology:" + technology + "\n")
    f.write("gid:" + country['ID Level'] + "\n")
    f.close()
