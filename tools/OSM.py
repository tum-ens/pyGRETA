import os
import zipfile
import urllib.request
import pandas as pd
import geopandas as gpd
import numpy as np

filter = True  # Activate saving reduced/filtered shapefiles


def extractShapefile(file_path, shapefile_name, output_filename):
    with zipfile.ZipFile(file_path, "r") as zipObj:
        output_file = output_filename + ".shp"
        if not os.path.isfile(output_file):
            with open(output_file, "wb") as f:
                f.write(zipObj.read(shapefile_name + ".shp"))

        output_file = output_filename + ".shx"
        if not os.path.isfile(output_file):
            with open(output_file, "wb") as f:
                f.write(zipObj.read(shapefile_name + ".shx"))

        output_file = output_filename + ".dbf"
        if not os.path.isfile(output_file):
            with open(output_file, "wb") as f:
                f.write(zipObj.read(shapefile_name + ".dbf"))

        output_file = output_filename + ".prj"
        if not os.path.isfile(output_file):
            with open(output_file, "wb") as f:
                f.write(zipObj.read(shapefile_name + ".prj"))


def filterShapefile(input_file, output_file, classes):
    if not os.path.isfile(output_file):
        gdf = gpd.read_file(input_file)

        selection = np.zeros((len(gdf),), dtype=bool)
        for element in classes:
            selection = selection | (gdf["fclass"].to_numpy() == element)

        if selection.sum():
            gdf_filtered = gdf[selection]
            gdf_filtered.to_file(output_file)
        else:
            print("filtering unnecessary")


def downloadShapefile(filename):
    if not os.path.exists("zip"):
        os.mkdir("zip")

    if not os.path.isfile("zip/" + filename):
        try:
            urllib.request.urlretrieve(url, "zip/" + filename)
            print("Successfully downloaded: " + filename)
        except:
            print(print("wrong filename: " + filename))


print("Started ... - Version 6")

# Import links to Geofabrik data
ISO = pd.read_csv("../assumptions/ISO3166-1.csv", sep=";")

# Iterate over each country
for country in ISO["Country"]:
    # Download shapefiles from Geofabrik
    url = ISO[ISO["Country"] == country]["shp.zip"].iloc[0]
    if type(url) == str:  # download only countries with url
        urls = url.split(",")
        if len(urls) == 1:

            country_letterCode = ISO[ISO["Country"] == country]["Alpha-3"].iloc[0]
            print('Started: ' + country + ", " + country_letterCode)

            # Download shapefile as *.zip if not already available
            filename = country_letterCode + '-latest-free.shp.zip'
            downloadShapefile(filename)

            # Unzip
            file_path = os.path.join("zip/", filename)
            if not os.path.exists("shapefiles"):
                os.mkdir("shapefiles")
            if not os.path.exists("shapefiles/railways"):
                os.mkdir("shapefiles/railways")

            shapefile_name = "gis_osm_landuse_a_free_1"
            output_filename = "shapefiles/" + country_letterCode + "-landuse"
            extractShapefile(file_path, shapefile_name, output_filename)

            shapefile_name = "gis_osm_roads_free_1"
            output_filename = "shapefiles/" + country_letterCode + "-roads"
            extractShapefile(file_path, shapefile_name, output_filename)

            shapefile_name = "gis_osm_railways_free_1"
            output_filename = "shapefiles/railways/" + country_letterCode + "-railways"
            extractShapefile(file_path, shapefile_name, output_filename)

            # Filter classes
            if filter:
                if not os.path.exists("shapefiles/filtered"):
                    os.mkdir("shapefiles/filtered")

                # Landuse
                # print('landuse')
                input_file = "shapefiles/" + country_letterCode + "-landuse.shp"
                output_file = "shapefiles/filtered/" + country_letterCode + "-landuse.shp"
                # "fclass" ILIKE 'commercial' OR "fclass" ILIKE 'industrial' OR "fclass" ILIKE 'quarry' OR "fclass" ILIKE 'military' OR "fclass" ILIKE 'park' OR "fclass" ILIKE 'recreation_ground'
                classes = ['commercial', 'industrial', 'quarry', 'military', 'park', 'recreation_ground']
                filterShapefile(input_file, output_file, classes)

                # Roads
                # print('roads')
                input_file = "shapefiles/" + country_letterCode + "-roads.shp"
                output_file = "shapefiles/filtered/" + country_letterCode + "-roads.shp"
                # "fclass" ILIKE 'primary' OR "fclass" ILIKE 'secondary'
                classes = ['primary', 'secondary']
                filterShapefile(input_file, output_file, classes)

print("done!")
