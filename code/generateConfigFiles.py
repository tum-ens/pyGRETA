import pandas as pd

year = 2019
technology = "WindOn,PV"
path = "configs/"

df = pd.read_excel('World-Timeplan.xlsx', sheet_name='Time Plan', engine='openpyxl')
# print(df)
# print(df['Country Code'])

for i in range(len(df)):
    f = open(path + str(i) + "_" + df['Country Name'][i] + ".txt", "w")
    f.write("spatial_scope:gadm36_" + str(df['Country Code'][i]) + "_0.shp" + "\n")
    f.write("subregions:gadm36_" + str(df['Country Code'][i]) + "_0.shp" + "\n")
    f.write("region_name:" + df['Country Name'][i] + "_gwa" + "\n")
    f.write("subregions_name:" + df['Country Name'][i] + "_level0" + "\n")
    f.write("country_code:" + str(df['Country Code'][i]) + "\n")
    f.write("year:" + str(year) + "\n")
    f.write("technology:" + technology + "\n")
    f.close()