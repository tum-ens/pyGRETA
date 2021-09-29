import pandas as pd
import os

year = 2019
technology = "OpenFieldPV,RoofTopPV,WindOn"
path = "../configs/"

ISO = pd.read_csv("../assumptions/ISO3166-1.csv", sep=";")

# Remove all old files
for f in os.listdir(path):
    os.remove(os.path.join(path, f))

# Iterate over each country
for i in range(len(ISO)):
    country = ISO.iloc[i]
    f = open(path + str(i) + "_" + country['Name']+ ".txt", "w")
    f.write("regions:gadm36_" + str(country['Alpha-3']) + "_0.shp" + "\n")
    f.write("region_name:" + country['Name'] + "\n")
    f.write("subregions_name:" + country['Name'] + "_level0" + "\n")
    f.write("country_code:" + str(country['Alpha-3']) + "\n")
    f.write("year:" + str(year) + "\n")
    f.write("technology:" + technology + "\n")
    f.close()