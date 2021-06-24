import pandas as pd
import numpy as np
import os

year = 2019
technology = "PV,WindOn"
path = "configs/"

df = pd.read_excel('World-Timeplan.xlsx', sheet_name='Time Plan', engine='openpyxl')
# print(df)
# print(df['Country Code'])

for f in os.listdir(path):
    os.remove(os.path.join(path, f))    # Remove all old files

for i in range(len(df)):
    if (not np.isnan(df['ID'][i]) and str(df['Country Name'][i])!='nan'):
        f = open(path + str(int(df['ID'][i])) + "_" + df['Country Name'][i] + ".txt", "w")
        f.write("spatial_scope:gadm36_" + str(df['Country Code'][i]) + "_0.shp" + "\n")
        f.write("subregions:gadm36_" + str(df['Country Code'][i]) + "_0.shp" + "\n")
        f.write("region_name:" + df['Country Name'][i] + "_gwa" + "\n")
        f.write("subregions_name:" + df['Country Name'][i] + "_level0" + "\n")
        f.write("country_code:" + str(df['Country Code'][i]) + "\n")
        f.write("year:" + str(year) + "\n")
        f.write("technology:" + technology + "\n")
        f.close()