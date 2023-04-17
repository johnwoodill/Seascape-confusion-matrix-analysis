import numpy as np
import pandas as pd
import netCDF4 as nc
import geopandas as gpd
import xarray as xr
import rasterio.mask
import rioxarray as rxr
import glob

from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, classification_report
from scipy.stats import mode
from shapely.geometry import Polygon


# Define a function to download data from a URL
def download_CW(url, filename):
    from urllib.request import urlretrieve
    urlretrieve(url, filename)


def get_geom(xmin, ymin, xmax, ymax):
    # Create the polygon
    polygon = Polygon([(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)])

    # Create a GeoDataFrame with the polygon as the geometry
    gdf = gpd.GeoDataFrame(geometry=[polygon])

    return gdf


def download_seascapes():
    # Download the data
    download_CW("https://cwcgom.aoml.noaa.gov/thredds/ncss/SEASCAPE_MONTH/SEASCAPES.nc?var=CLASS&var=P&north=90&west=-180&east=180&south=-90&disableProjSubset=on&horizStride=1&time_start=2017-01-01T12%3A00%3A00Z&time_end=2017-12-15T12%3A00%3A00Z&timeStride=1&addLatLon=true&accept=netcdf", "data/seascapes/MODISSS_2017.nc")
    print("Saving: data/seascapes/MODISSS_2017.nc")

    download_CW("https://cwcgom.aoml.noaa.gov/thredds/ncss/SEASCAPE_MONTH/SEASCAPES.nc?var=CLASS&var=P&north=90&west=-180&east=180&south=-90&disableProjSubset=on&horizStride=1&time_start=2018-01-01T12%3A00%3A00Z&time_end=2018-12-15T12%3A00%3A00Z&timeStride=1&addLatLon=true&accept=netcdf", "data/seascapes/MODISSS_2018.nc")
    print("Saving: data/seascapes/MODISSS_2018.nc")
    
    download_CW("https://cwcgom.aoml.noaa.gov/thredds/ncss/SEASCAPE_MONTH/SEASCAPES.nc?var=CLASS&var=P&north=90&west=-180&east=180&south=-90&disableProjSubset=on&horizStride=1&time_start=2019-01-01T12%3A00%3A00Z&time_end=2019-12-15T12%3A00%3A00Z&timeStride=1&addLatLon=true&accept=netcdf", "data/seascapes/MODISSS_2019.nc")
    print(f"Saving: data/seascapes/MODISSS_2019.nc")
    
    download_CW("https://cwcgom.aoml.noaa.gov/thredds/ncss/SEASCAPE_MONTH/SEASCAPES.nc?var=CLASS&var=P&north=90&west=-180&east=180&south=-90&disableProjSubset=on&horizStride=1&time_start=2020-01-01T12%3A00%3A00Z&time_end=2020-12-15T12%3A00%3A00Z&timeStride=1&addLatLon=true&accept=netcdf", "data/seascapes/MODISSS_2020.nc")
    print(f"Saving: data/seascapes/MODISSS_2020.nc")
    
    download_CW("https://cwcgom.aoml.noaa.gov/thredds/ncss/SEASCAPE_MONTH/SEASCAPES.nc?var=CLASS&var=P&north=90&west=-180&east=180&south=-90&disableProjSubset=on&horizStride=1&time_start=2021-01-01T12%3A00%3A00Z&time_end=2021-12-15T12%3A00%3A00Z&timeStride=1&addLatLon=true&accept=netcdf", "data/seascapes/MODISSS_2021.nc")
    print(f"Saving: data/seascapes/MODISSS_2021.nc")

    download_CW("https://cwcgom.aoml.noaa.gov/thredds/ncss/SEASCAPE_MONTH/SEASCAPES.nc?var=CLASS&var=P&north=90&west=-180&east=180&south=-90&disableProjSubset=on&horizStride=1&time_start=2022-01-01T12%3A00%3A00Z&time_end=2022-12-15T12%3A00%3A00Z&timeStride=1&addLatLon=true&accept=netcdf", "data/seascapes/MODISSS_2022.nc")
    print(f"Saving: data/seascapes/MODISSS_2022.nc")

    download_CW("https://cwcgom.aoml.noaa.gov/thredds/ncss/SEASCAPE_MONTH_VIIRS_OLCI/SEASCAPES.nc?var=CLASS&var=P&north=90&west=-180&east=180&south=-90&disableProjSubset=on&horizStride=1&time_start=2017-01-01T12%3A00%3A00Z&time_end=2017-12-15T12%3A00%3A00Z&timeStride=1&addLatLon=true&accept=netcdf", "data/seascapes/VIIRSSS_2017.nc")
    print(f"Saving: data/seascapes/VIIRSSS_2017.nc")
    
    download_CW("https://cwcgom.aoml.noaa.gov/thredds/ncss/SEASCAPE_MONTH_VIIRS_OLCI/SEASCAPES.nc?var=CLASS&var=P&north=90&west=-180&east=180&south=-90&disableProjSubset=on&horizStride=1&time_start=2018-01-01T12%3A00%3A00Z&time_end=2018-12-15T12%3A00%3A00Z&timeStride=1&addLatLon=true&accept=netcdf", "data/seascapes/VIIRSSS_2018.nc")
    print(f"Saving: data/seascapes/VIIRSSS_2018.nc")

    download_CW("https://cwcgom.aoml.noaa.gov/thredds/ncss/SEASCAPE_MONTH_VIIRS_OLCI/SEASCAPES.nc?var=CLASS&var=P&north=90&west=-180&east=180&south=-90&disableProjSubset=on&horizStride=1&time_start=2019-01-01T12%3A00%3A00Z&time_end=2019-12-15T12%3A00%3A00Z&timeStride=1&addLatLon=true&accept=netcdf", "data/seascapes/VIIRSSS_2019.nc")
    print(f"Saving: data/seascapes/VIIRSSS_2019.nc")

    download_CW("https://cwcgom.aoml.noaa.gov/thredds/ncss/SEASCAPE_MONTH_VIIRS_OLCI/SEASCAPES.nc?var=CLASS&var=P&north=90&west=-180&east=180&south=-90&disableProjSubset=on&horizStride=1&time_start=2020-01-01T12%3A00%3A00Z&time_end=2020-12-15T12%3A00%3A00Z&timeStride=1&addLatLon=true&accept=netcdf", "data/seascapes/VIIRSSS_2020.nc")
    print(f"Saving: data/seascapes/VIIRSSS_2020.nc")

    download_CW("https://cwcgom.aoml.noaa.gov/thredds/ncss/SEASCAPE_MONTH_VIIRS_OLCI/SEASCAPES.nc?var=CLASS&var=P&north=90&west=-180&east=180&south=-90&disableProjSubset=on&horizStride=1&time_start=2021-01-01T12%3A00%3A00Z&time_end=2021-12-15T12%3A00%3A00Z&timeStride=1&addLatLon=true&accept=netcdf", "data/seascapes/VIIRSSS_2021.nc")
    print(f"Saving: data/seascapes/VIIRSSS_2021.nc")

    download_CW("https://cwcgom.aoml.noaa.gov/thredds/ncss/SEASCAPE_MONTH_VIIRS_OLCI/SEASCAPES.nc?var=CLASS&var=P&north=90&west=-180&east=180&south=-90&disableProjSubset=on&horizStride=1&time_start=2022-01-01T12%3A00%3A00Z&time_end=2022-12-15T12%3A00%3A00Z&timeStride=1&addLatLon=true&accept=netcdf", "data/seascapes/VIIRSSS_2022.nc")
    print(f"Saving: data/seascapes/VIIRSSS_2022.nc")
    

    
def proc_seascape_data(year):

    # Get shapefiles for mask
    # -----------------------------------------------------------------------------
    # Load main ocean boundary data
    oceans = gpd.read_file('data/shapefiles/world_oceans/World_Seas_IHO_v3.shp')

    # ca_coast = [-130, 32.5, -114.5, 42]
    ca_coast = get_geom(-130, 32.5, -114.5, 42)

    # east_coast = [-82, 25, -65, 45]
    east_coast = get_geom(-82, 25, -65, 45)

    # Additional oceans
    gulf_of_mexico = oceans[oceans['NAME'] == 'Gulf of Mexico']
    north_pacific = oceans[oceans['NAME'] == 'North Pacific Ocean']
    south_pacific =oceans[oceans['NAME'] == 'South Pacific Ocean']
    north_atlantic = oceans[oceans['NAME'] == 'North Atlantic Ocean']
    south_atlantic = oceans[oceans['NAME'] == 'South Atlantic Ocean']
    indian_ocean = oceans[oceans['NAME'] == 'Indian Ocean']
    
    # Load data
    mdat = rxr.open_rasterio(f"data/seascapes/MODISSS_{year}.nc", crs='EPSG:4326')
    vdat = rxr.open_rasterio(f"data/seascapes/VIIRSSS_{year}.nc", crs='EPSG:4326')

    # Rename columns
    mdat = mdat.rename({'CLASS': 'MCLASS', 'P': 'MP'})
    vdat = vdat.rename({'CLASS': 'VCLASS', 'P': 'VP'})

    # Set crs
    mdat = mdat.rio.write_crs("EPSG:4326")
    vdat = vdat.rio.write_crs("EPSG:4326")

    # Check crs
    # print(mdat.rio.crs)
    # print(vdat.rio.crs)

    # Merge two netcdf files and then get difference
    gdf = xr.merge([mdat, vdat])
    
    # Mask CA Coast
    geom = ca_coast
    outdat = gdf.rio.clip(geom.geometry, all_touched=False, drop=True, invert=False, from_disk=True) 
    outdat.to_netcdf(f"data/processed/oceans/ca_coast_{year}.nc")
    print(f"Saving: data/processed/oceans/ca_coast_{year}.nc")
    outdat2 = outdat.to_dataframe().reset_index()
    outdat2.to_csv(f"data/processed/oceans/ca_coast_{year}.csv", index=False)
    print(f"Saving: data/processed/oceans/ca_coast_{year}.csv")
    

    # Mask East Coast
    geom = east_coast
    outdat = gdf.rio.clip(geom.geometry, all_touched=False, drop=True, invert=False, from_disk=True) 
    outdat.to_netcdf(f"data/processed/oceans/east_coast_{year}.nc")
    print(f"Saving: data/processed/oceans/east_coast_{year}.nc")
    outdat2 = outdat.to_dataframe().reset_index()
    outdat2.to_csv(f"data/processed/oceans/east_coast_{year}.csv", index=False)
    print(f"Saving: data/processed/oceans/east_coast_{year}.csv")

    # Mask gulf of mexico
    geom = gulf_of_mexico
    outdat = gdf.rio.clip(geom.geometry, all_touched=False, drop=True, invert=False, from_disk=True) 
    outdat.to_netcdf(f"data/processed/oceans/gulf_of_mexico_{year}.nc")
    print(f"Saving: data/processed/oceans/gulf_of_mexico_{year}.nc")
    outdat2 = outdat.to_dataframe().reset_index()
    outdat2.to_csv(f"data/processed/oceans/gulf_of_mexico_{year}.csv", index=False)
    print(f"Saving: data/processed/oceans/gulf_of_mexico_{year}.csv")


    # Mask North Pacific
    geom = north_pacific
    outdat = gdf.rio.clip(geom.geometry, all_touched=False, drop=True, invert=False, from_disk=True) 
    outdat.to_netcdf(f"data/processed/oceans/north_pacific_{year}.nc")
    print(f"Saving: data/processed/oceans/north_pacific_{year}.nc")
    outdat2 = outdat.to_dataframe().reset_index()
    outdat2.to_csv(f"data/processed/oceans/north_pacific_{year}.csv", index=False)
    print(f"Saving: data/processed/oceans/north_pacific_{year}.csv")

    # Mask South Pacific
    geom = south_pacific
    outdat = gdf.rio.clip(geom.geometry, all_touched=False, drop=True, invert=False, from_disk=True) 
    outdat.to_netcdf(f"data/processed/oceans/south_pacific_{year}.nc")
    print(f"Saving: data/processed/oceans/south_pacific_{year}.nc")
    outdat2 = outdat.to_dataframe().reset_index()
    outdat2.to_csv(f"data/processed/oceans/south_pacific_{year}.csv", index=False)
    print(f"Saving: data/processed/oceans/south_pacific_{year}.csv")

    # Mask North Atlantic
    geom = north_atlantic
    outdat = gdf.rio.clip(geom.geometry, all_touched=False, drop=True, invert=False, from_disk=True) 
    outdat.to_netcdf(f"data/processed/oceans/north_atlantic_{year}.nc")
    print(f"Saving: data/processed/oceans/north_atlantic_{year}.nc")
    outdat2 = outdat.to_dataframe().reset_index()
    outdat2.to_csv(f"data/processed/oceans/north_atlantic_{year}.csv", index=False)
    print(f"Saving: data/processed/oceans/north_atlantic_{year}.csv")

    # Mask South Atlantic
    geom = south_atlantic
    outdat = gdf.rio.clip(geom.geometry, all_touched=False, drop=True, invert=False, from_disk=True) 
    outdat.to_netcdf(f"data/processed/oceans/south_atlantic_{year}.nc")
    print(f"Saving: data/processed/oceans/south_atlantic_{year}.nc")
    outdat2 = outdat.to_dataframe().reset_index()
    outdat2.to_csv(f"data/processed/oceans/south_atlantic_{year}.csv", index=False)
    print(f"Saving: data/processed/oceans/south_atlantic_{year}.csv")


    # Mask Indian Ocean
    geom = indian_ocean.geometry
    outdat = gdf.rio.clip(geom.geometry, all_touched=False, drop=True, invert=False, from_disk=True) 
    outdat.to_netcdf(f"data/processed/oceans/indian_ocean_{year}.nc")
    print(f"Saving: data/processed/oceans/indian_ocean_{year}.nc")
    outdat2 = outdat.to_dataframe().reset_index()
    outdat2.to_csv(f"data/processed/oceans/indian_ocean_{year}.csv", index=False)
    print(f"Saving: data/processed/oceans/indian_ocean_{year}.csv")






file_ = "data/processed/oceans/ca_coast_2018.nc"

def proc_conf_mats(file_):

        filename = file_.split("/")[-1].split(".")[0]
        region = "_".join(file_.split("/")[-1].split(".")[0].split("_")[0:-1])
        year = file_.split("/")[-1].split(".")[0].split("_")[-1]

        dat = rxr.open_rasterio(file_)

        indat = dat.to_dataframe().reset_index()

        true_labels = indat['MCLASS']
        pred_labels = indat['VCLASS']

        # create a list of all possible labels
        all_labels = np.unique(np.concatenate((true_labels, pred_labels)))

        # create the confusion matrix using the true and predicted labels
        cm = confusion_matrix(true_labels, pred_labels, labels=all_labels)

        # calculate the accuracy for each class
        accuracies = {}
        for i, label in enumerate(all_labels):
            total = sum(cm[i, :])
            correct = cm[i, i]
            accuracy = correct / total
            accuracies[label] = accuracy

        df_cm = pd.DataFrame(cm, index=all_labels, columns=all_labels).reset_index()
        df_cm = df_cm.rename(columns={'index': 'seascape'})
        df_cm.insert(1, 'acc', pd.DataFrame(accuracies, index=[0]).T)
        df_cm.insert(2, 'region', region)
        df_cm.insert(3, 'year', year)
        df_cm = df_cm.sort_values('acc', ascending=False).reset_index(drop=True)

        full_accuracy = accuracy_score(true_labels, pred_labels)
        precision = precision_score(true_labels, pred_labels, average=None, labels=np.unique(true_labels))
        recall = recall_score(true_labels, pred_labels, average=None, labels=np.unique(true_labels))
        f1 = f1_score(true_labels, pred_labels, average=None, labels=np.unique(true_labels))

        amat = pd.DataFrame({'seascape': np.unique(true_labels)})
        amat['full_accuracy'] = full_accuracy
        amat['precision'] = precision
        amat['recall'] = recall
        amat['f1'] = f1

        df_cm.merge(amat, on='seascape', how='left')

        df_cm.to_csv(f"data/processed/conf_mats/{filename}.csv", index=False)


[proc_conf_mats(x) for x in glob.glob("data/processed/oceans/ca_coast*")]




if __name__ == "__main__":

    # download_seascapes()

    years_ = ["2017", "2018", "2019", "2020", "2021", "2022"]

    # proc_seascape_data("2017")

    results = [proc_seascape_data(years) for years in years_]









