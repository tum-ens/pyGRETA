import os
import zipfile
import urllib.request
import pandas as pd
import geopandas as gpd
import numpy as np

filter = True  # Activate saving reduced/filtered shapefiles


def downloadShapefile(url, filename):
    if not os.path.exists("zip"):
        os.mkdir("zip")

    if not os.path.isfile("zip/" + filename):
        try:
            urllib.request.urlretrieve(url, "zip/" + filename)
            print("Successfully downloaded: " + filename)
        except:
            print(print("wrong filename: " + filename))


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

        output_file = output_filename + ".cpg"
        if not os.path.isfile(output_file):
            with open(output_file, "wb") as f:
                f.write(zipObj.read(shapefile_name + ".cpg"))


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


def mergeShapefiles(input_files, output_file):
    if not os.path.isfile(output_file):
        gdf_merged = gpd.GeoDataFrame(pd.concat([gpd.read_file(i) for i in input_files], ignore_index=True), crs=gpd.read_file(input_files[0]).crs)
        gdf_merged.to_file(output_file)


print("Started ... - Version 7")

# Import links to Geofabrik data
ISO = pd.read_csv("../assumptions/ISO3166-1.csv", sep=";")

# Iterate over each country
for country in ISO["Country"]:
    # Download shapefiles from Geofabrik
    url = ISO[ISO["Country"] == country]["shp.zip"].iloc[0]
    if type(url) == str:  # download only countries with url
        country_letterCode = ISO[ISO["Country"] == country]["Alpha-3"].iloc[0]
        urls = url.split(",")
        numberOfparts = len(urls)
        print('Started: ' + country + ", " + country_letterCode + " - parts: " + numberOfparts)

        for i, link in enumerate(urls):

            if len(urls) == 1:
                i = ''  # remove additional number in shapefile names for unsplitted countries

            # Download shapefile as *.zip if not already available
            filename = country_letterCode + str(i) + '-latest-free.shp.zip'
            downloadShapefile(link, filename)

            # Unzip downloaded shapefiles
            file_path = os.path.join("zip/", filename)
            if not os.path.exists("shapefiles"):
                os.mkdir("shapefiles")

            shapefile_name = "gis_osm_landuse_a_free_1"
            output_filename = "shapefiles/" + country_letterCode + str(i) + "-landuse"
            extractShapefile(file_path, shapefile_name, output_filename)

            shapefile_name = "gis_osm_roads_free_1"
            output_filename = "shapefiles/" + country_letterCode + str(i) + "-roads"
            extractShapefile(file_path, shapefile_name, output_filename)

            if not os.path.exists("shapefiles/railways"):
                os.mkdir("shapefiles/railways")
            shapefile_name = "gis_osm_railways_free_1"
            output_filename = "shapefiles/railways/" + country_letterCode + str(i) + "-railways"
            extractShapefile(file_path, shapefile_name, output_filename)


            # Remove unused classes in shapefiles
            if filter:
                print("filtering")
                if not os.path.exists("shapefiles/filtered"):
                    os.mkdir("shapefiles/filtered")

                # Landuse
                input_file = "shapefiles/" + country_letterCode + str(i) + "-landuse.shp"
                output_file = "shapefiles/filtered/" + country_letterCode + str(i) + "-landuse.shp"
                # "fclass" ILIKE 'commercial' OR "fclass" ILIKE 'industrial' OR "fclass" ILIKE 'quarry' OR "fclass" ILIKE 'military' OR "fclass" ILIKE 'park' OR "fclass" ILIKE 'recreation_ground'
                classes = ['commercial', 'industrial', 'quarry', 'military', 'park', 'recreation_ground']
                filterShapefile(input_file, output_file, classes)

                # Roads
                input_file = "shapefiles/" + country_letterCode + str(i) + "-roads.shp"
                output_file = "shapefiles/filtered/" + country_letterCode + str(i) + "-roads.shp"
                # "fclass" ILIKE 'primary' OR "fclass" ILIKE 'secondary'
                classes = ['primary', 'secondary']
                filterShapefile(input_file, output_file, classes)

        # Merge if country is splitted into multiple shapefiles
        if numberOfparts >= 1:
            print("merging")
            input_file = ["shapefiles/filtered/" + country_letterCode + str(i) + "-landuse.shp" for i in range(numberOfparts)]
            output_file = "shapefiles/filtered/" + country_letterCode + "-landuse.shp"
            mergeShapefiles(input_file, output_file)

            input_file = ["shapefiles/filtered/" + country_letterCode + str(i) + "-roads.shp" for i in range(numberOfparts)]
            output_file = "shapefiles/filtered/" + country_letterCode + "-roads.shp"
            mergeShapefiles(input_file, output_file)

            input_file = ["shapefiles/railways/" + country_letterCode + str(i) + "-railways.shp" for i in range(numberOfparts)]
            output_file = "shapefiles/railways/" + country_letterCode + "-railways.shp"
            mergeShapefiles(input_file, output_file)


print("done!")
