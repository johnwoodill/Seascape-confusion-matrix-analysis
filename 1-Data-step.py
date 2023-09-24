import numpy as np
import pandas as pd
# import netCDF4 as nc
import geopandas as gpd
import xarray as xr
import rasterio.mask
import rioxarray as rxr
import glob

from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, classification_report
from scipy.stats import mode
from urllib.request import urlretrieve

from shapely.geometry import Polygon
from shapely.ops import unary_union


# Define a function to download data from a URL
def download_CW(url, filename):
    urlretrieve(url, filename)


def get_geom(xmin, ymin, xmax, ymax):
    # Create the polygon
    polygon = Polygon([(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)])

    # Create a GeoDataFrame with the polygon as the geometry
    gdf = gpd.GeoDataFrame(geometry=[polygon])

    return gdf


def download_monthly_seascapes():
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



    
def download_8day_seascapes():
    for year_ in np.arange(2017, 2022):
        print(year_)

        download_CW(f"https://cwcgom.aoml.noaa.gov/erddap/griddap/noaa_aoml_seascapes_8day.nc?CLASS[({year_}-01-01):1:({year_}-04-30)][(-80.00):1:(80.00)][(-179.975):1:(179.975)],P[({year_}-01-01):1:({year_}-04-30)][(-80.00):1:(80.00)][(-179.975):1:(179.975)]", f"data/seascapes/MODISSS_8DAY_0104-{year_}.nc")
        print(f"Saving: data/seascapes/MODISSS_8DAY_0104-{year_}.nc")

        download_CW(f"https://cwcgom.aoml.noaa.gov/erddap/griddap/noaa_aoml_seascapes_8day.nc?CLASS[({year_}-05-01):1:({year_}-08-31)][(-80.00):1:(80.00)][(-179.975):1:(179.975)],P[({year_}-05-01):1:({year_}-08-31)][(-80.00):1:(80.00)][(-179.975):1:(179.975)]", f"data/seascapes/MODISSS_8DAY_0508-{year_}.nc")
        print(f"Saving: data/seascapes/MODISSS_8DAY_0508-{year_}.nc")

        download_CW(f"https://cwcgom.aoml.noaa.gov/erddap/griddap/noaa_aoml_seascapes_8day.nc?CLASS[({year_}-09-01):1:({year_}-12-31)][(-80.00):1:(80.00)][(-179.975):1:(179.975)],P[({year_}-09-01):1:({year_}-12-31)][(-80.00):1:(80.00)][(-179.975):1:(179.975)]", f"data/seascapes/MODISSS_8DAY_0912-{year_}.nc")
        print(f"Saving: data/seascapes/MODISSS_8DAY_0912-{year_}.nc")

        download_CW(f"https://cwcgom.aoml.noaa.gov/erddap/griddap/noaa_aoml_seascapes_8day_viirs_olci.nc?CLASS[({year_}-01-01):1:({year_}-04-30)][(-80.00):1:(80.00)][(-179.975):1:(179.975)],P[({year_}-01-01):1:({year_}-04-30)][(-80.00):1:(80.00)][(-179.975):1:(179.975)]", f"data/seascapes/VIIRSSS_8DAY_0104-{year_}.nc")
        print(f"Saving: data/seascapes/VIIRSSS_8DAY_0104-{year_}.nc")

        download_CW(f"https://cwcgom.aoml.noaa.gov/erddap/griddap/noaa_aoml_seascapes_8day_viirs_olci.nc?CLASS[({year_}-05-01):1:({year_}-08-31)][(-80.00):1:(80.00)][(-179.975):1:(179.975)],P[({year_}-05-01):1:({year_}-08-31)][(-80.00):1:(80.00)][(-179.975):1:(179.975)]", f"data/seascapes/VIIRSSS_8DAY_0508-{year_}.nc")
        print(f"Saving: data/seascapes/VIIRSSS_8DAY_0508-{year_}.nc")

        download_CW(f"https://cwcgom.aoml.noaa.gov/erddap/griddap/noaa_aoml_seascapes_8day_viirs_olci.nc?CLASS[({year_}-09-01):1:({year_}-12-31)][(-80.00):1:(80.00)][(-179.975):1:(179.975)],P[({year_}-09-01):1:({year_}-12-31)][(-80.00):1:(80.00)][(-179.975):1:(179.975)]", f"data/seascapes/VIIRSSS_8DAY_0912-{year_}.nc")
        print(f"Saving: data/seascapes/VIIRSSS_8DAY_0912-{year_}.nc")

    
def proc_monthly_seascape_data(year):

    # Get shapefiles for mask
    # -----------------------------------------------------------------------------
    # Load main ocean boundary data
    oceans = gpd.read_file('data/shapefiles/world_oceans/World_Seas_IHO_v3.shp')

    # ca_coast = [-130, 32.5, -114.5, 42]
    ca_coast = get_geom(-130, 32.5, -114.5, 42)

    # east_coast = [-82, 25, -65, 45]
    east_coast = get_geom(-82, 25, -65, 45)

    # Oregon Coast
    pnw_coast = get_geom(-130, 42, -114.5, 49)

    # Additional oceans
    gulf_of_mexico = oceans[oceans['NAME'] == 'Gulf of Mexico']
    north_pacific = oceans[oceans['NAME'] == 'North Pacific Ocean']
    south_pacific =oceans[oceans['NAME'] == 'South Pacific Ocean']
    north_atlantic = oceans[oceans['NAME'] == 'North Atlantic Ocean']
    south_atlantic = oceans[oceans['NAME'] == 'South Atlantic Ocean']
    indian_ocean = oceans[oceans['NAME'] == 'Indian Ocean']
    chukchi_sea = oceans[oceans['NAME'] == 'Chukchi Sea']
    beaufort_sea = oceans[oceans['NAME'] == 'Beaufort Sea']
    
    chukchi_sea = chukchi_sea.explode()
    chukchi_sea = chukchi_sea.iloc[[0], :]

    chukchi_sea.bounds.minx.min()
    chukchi_sea.bounds.maxx.max()

    chukchi_sea.bounds.miny.min()
    chukchi_sea.bounds.maxy.max()

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
    
    # Mask PNW Coast
    geom = pnw_coast
    outdat = gdf.rio.clip(geom.geometry, all_touched=False, drop=True, invert=False, from_disk=True) 
    outdat.to_netcdf(f"data/processed/oceans/pnw_coast_{year}.nc")
    print(f"Saving: data/processed/oceans/pnw_coast_{year}.nc")
    outdat2 = outdat.to_dataframe().reset_index()
    outdat2.to_csv(f"data/processed/oceans/pnw_coast_{year}.csv", index=False)
    print(f"Saving: data/processed/oceans/pnw_coast_{year}.csv")


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
    geom = indian_ocean
    outdat = gdf.rio.clip(geom.geometry, all_touched=False, drop=True, invert=False, from_disk=True) 
    outdat.to_netcdf(f"data/processed/oceans/indian_ocean_{year}.nc")
    print(f"Saving: data/processed/oceans/indian_ocean_{year}.nc")
    outdat2 = outdat.to_dataframe().reset_index()
    outdat2.to_csv(f"data/processed/oceans/indian_ocean_{year}.csv", index=False)
    print(f"Saving: data/processed/oceans/indian_ocean_{year}.csv")

    # Mask Chukchi Beaufort sea Ocean
    geom = chukchi_sea
    outdat = gdf.rio.clip(geom.geometry, all_touched=False, drop=True, invert=False, from_disk=True) 
    outdat.to_netcdf(f"data/processed/oceans/chukchi_sea_{year}.nc")
    print(f"Saving: data/processed/oceans/chukchi_sea_{year}.nc")
    outdat2 = outdat.to_dataframe().reset_index()
    outdat2.to_csv(f"data/processed/oceans/chukchi_sea_{year}.csv", index=False)
    print(f"Saving: data/processed/oceans/chukchi_sea_{year}.csv")

    geom = beaufort_sea
    outdat = gdf.rio.clip(geom.geometry, all_touched=False, drop=True, invert=False, from_disk=True) 
    outdat.to_netcdf(f"data/processed/oceans/beaufort_sea_{year}.nc")
    print(f"Saving: data/processed/oceans/beaufort_sea_{year}.nc")
    outdat2 = outdat.to_dataframe().reset_index()
    outdat2.to_csv(f"data/processed/oceans/beaufort_sea_{year}.csv", index=False)
    print(f"Saving: data/processed/oceans/beaufort_sea_{year}.csv")


year = "2017"
year = "2018"
proc_8day_seascape_data("2019")
proc_8day_seascape_data("2020")
proc_8day_seascape_data("2021")

def proc_8day_seascape_data(year):

    # Get shapefiles for mask
    # -----------------------------------------------------------------------------
    # Load main ocean boundary data
    oceans = gpd.read_file('data/shapefiles/world_oceans/World_Seas_IHO_v3.shp')

    # ca_coast = [-130, 32.5, -114.5, 42]
    ca_coast = get_geom(-130, 32.5, -114.5, 42)

    # east_coast = [-82, 25, -65, 45]
    east_coast = get_geom(-82, 25, -65, 45)

    # Oregon Coast
    pnw_coast = get_geom(-130, 42, -114.5, 49)

    # Additional oceans
    gulf_of_mexico = oceans[oceans['NAME'] == 'Gulf of Mexico']
    north_pacific = oceans[oceans['NAME'] == 'North Pacific Ocean']
    south_pacific =oceans[oceans['NAME'] == 'South Pacific Ocean']
    north_atlantic = oceans[oceans['NAME'] == 'North Atlantic Ocean']
    south_atlantic = oceans[oceans['NAME'] == 'South Atlantic Ocean']
    indian_ocean = oceans[oceans['NAME'] == 'Indian Ocean']
    chukchi_sea = oceans[oceans['NAME'] == 'Chukchi Sea']
    beaufort_sea = oceans[oceans['NAME'] == 'Beaufort Sea']
    
    chukchi_sea = chukchi_sea.explode()
    chukchi_sea = chukchi_sea.iloc[[0], :]

    chukchi_sea.bounds.minx.min()
    chukchi_sea.bounds.maxx.max()

    chukchi_sea.bounds.miny.min()
    chukchi_sea.bounds.maxy.max()

    # Load data
    mdat1 = rxr.open_rasterio(f"data/seascapes/MODISSS_8DAY_0104-{year}.nc", crs='EPSG:4326')
    mdat2 = rxr.open_rasterio(f"data/seascapes/MODISSS_8DAY_0508-{year}.nc", crs='EPSG:4326')
    mdat3 = rxr.open_rasterio(f"data/seascapes/MODISSS_8DAY_0912-{year}.nc", crs='EPSG:4326')

    vdat1 = rxr.open_rasterio(f"data/seascapes/VIIRSSS_8DAY_0104-{year}.nc", crs='EPSG:4326')
    vdat2 = rxr.open_rasterio(f"data/seascapes/VIIRSSS_8DAY_0508-{year}.nc", crs='EPSG:4326')
    vdat3 = rxr.open_rasterio(f"data/seascapes/VIIRSSS_8DAY_0912-{year}.nc", crs='EPSG:4326')

    # Rename columns
    mdat1 = mdat1.rename({'CLASS': 'MCLASS', 'P': 'MP'})
    mdat2 = mdat2.rename({'CLASS': 'MCLASS', 'P': 'MP'})
    mdat3 = mdat3.rename({'CLASS': 'MCLASS', 'P': 'MP'})

    vdat1 = vdat1.rename({'CLASS': 'VCLASS', 'P': 'VP'})
    vdat2 = vdat2.rename({'CLASS': 'VCLASS', 'P': 'VP'})
    vdat3 = vdat3.rename({'CLASS': 'VCLASS', 'P': 'VP'})

    # Set crs
    mdat1 = mdat1.rio.write_crs("EPSG:4326")
    mdat2 = mdat2.rio.write_crs("EPSG:4326")
    mdat3 = mdat2.rio.write_crs("EPSG:4326")

    vdat1 = vdat1.rio.write_crs("EPSG:4326")
    vdat2 = vdat2.rio.write_crs("EPSG:4326")
    vdat3 = vdat2.rio.write_crs("EPSG:4326")

    # Check crs
    # print(mdat.rio.crs)
    # print(vdat.rio.crs)

    # Merge netcdf files and then get difference
    gdf = xr.merge([mdat1, vdat1, mdat2, vdat2, mdat3, vdat3])

    # Mask CA Coast
    geom = ca_coast
    outdat = gdf.rio.clip(geom.geometry, all_touched=False, drop=True, invert=False, from_disk=True) 
    outdat.to_netcdf(f"data/processed/oceans/ca_coast_8day_{year}.nc")
    print(f"Saving: data/processed/oceans/ca_coast_8day_{year}.nc")
    outdat2 = outdat.to_dataframe().reset_index()
    outdat2.to_csv(f"data/processed/oceans/ca_coast_8day_{year}.csv", index=False)
    print(f"Saving: data/processed/oceans/ca_coast_8day_{year}.csv")
    
    # Mask PNW Coast
    geom = pnw_coast
    outdat = gdf.rio.clip(geom.geometry, all_touched=False, drop=True, invert=False, from_disk=True) 
    outdat.to_netcdf(f"data/processed/oceans/pnw_coast_8day_{year}.nc")
    print(f"Saving: data/processed/oceans/pnw_coast_8day_{year}.nc")
    outdat2 = outdat.to_dataframe().reset_index()
    outdat2.to_csv(f"data/processed/oceans/pnw_coast_8day_{year}.csv", index=False)
    print(f"Saving: data/processed/oceans/pnw_coast_8day_{year}.csv")

    # Mask East Coast
    geom = east_coast
    outdat = gdf.rio.clip(geom.geometry, all_touched=False, drop=True, invert=False, from_disk=True) 
    outdat.to_netcdf(f"data/processed/oceans/east_coast_8day_{year}.nc")
    print(f"Saving: data/processed/oceans/east_coast_8day_{year}.nc")
    outdat2 = outdat.to_dataframe().reset_index()
    outdat2.to_csv(f"data/processed/oceans/east_coast_8day_{year}.csv", index=False)
    print(f"Saving: data/processed/oceans/east_coast_8day_{year}.csv")

    # Mask gulf of mexico
    geom = gulf_of_mexico
    outdat = gdf.rio.clip(geom.geometry, all_touched=False, drop=True, invert=False, from_disk=True) 
    outdat.to_netcdf(f"data/processed/oceans/gulf_of_mexico_8day_{year}.nc")
    print(f"Saving: data/processed/oceans/gulf_of_mexico_8day_{year}.nc")
    outdat2 = outdat.to_dataframe().reset_index()
    outdat2.to_csv(f"data/processed/oceans/gulf_of_mexico_8day_{year}.csv", index=False)
    print(f"Saving: data/processed/oceans/gulf_of_mexico_8day_{year}.csv")

    # Mask North Pacific
    geom = north_pacific
    outdat = gdf.rio.clip(geom.geometry, all_touched=False, drop=True, invert=False, from_disk=True) 
    outdat.to_netcdf(f"data/processed/oceans/north_pacific_8day_{year}.nc")
    print(f"Saving: data/processed/oceans/north_pacific_8day_{year}.nc")
    outdat2 = outdat.to_dataframe().reset_index()
    outdat2.to_csv(f"data/processed/oceans/north_pacific_8day_{year}.csv", index=False)
    print(f"Saving: data/processed/oceans/north_pacific_8day_{year}.csv")

    # Mask South Pacific
    geom = south_pacific
    outdat = gdf.rio.clip(geom.geometry, all_touched=False, drop=True, invert=False, from_disk=True) 
    outdat.to_netcdf(f"data/processed/oceans/south_pacific_8day_{year}.nc")
    print(f"Saving: data/processed/oceans/south_pacific_8day_{year}.nc")
    outdat2 = outdat.to_dataframe().reset_index()
    outdat2.to_csv(f"data/processed/oceans/south_pacific_8day_{year}.csv", index=False)
    print(f"Saving: data/processed/oceans/south_pacific_8day_{year}.csv")

    # Mask North Atlantic
    geom = north_atlantic
    outdat = gdf.rio.clip(geom.geometry, all_touched=False, drop=True, invert=False, from_disk=True) 
    outdat.to_netcdf(f"data/processed/oceans/north_atlantic_8day_{year}.nc")
    print(f"Saving: data/processed/oceans/north_atlantic_8day_{year}.nc")
    outdat2 = outdat.to_dataframe().reset_index()
    outdat2.to_csv(f"data/processed/oceans/north_atlantic_8day_{year}.csv", index=False)
    print(f"Saving: data/processed/oceans/north_atlantic_8day_{year}.csv")

    # Mask South Atlantic
    geom = south_atlantic
    outdat = gdf.rio.clip(geom.geometry, all_touched=False, drop=True, invert=False, from_disk=True) 
    outdat.to_netcdf(f"data/processed/oceans/south_atlantic_8day_{year}.nc")
    print(f"Saving: data/processed/oceans/south_atlantic_8day_{year}.nc")
    outdat2 = outdat.to_dataframe().reset_index()
    outdat2.to_csv(f"data/processed/oceans/south_atlantic_8day_{year}.csv", index=False)
    print(f"Saving: data/processed/oceans/south_atlantic_8day_{year}.csv")

    # Mask Indian Ocean
    geom = indian_ocean
    outdat = gdf.rio.clip(geom.geometry, all_touched=False, drop=True, invert=False, from_disk=True) 
    outdat.to_netcdf(f"data/processed/oceans/indian_ocean_8day_{year}.nc")
    print(f"Saving: data/processed/oceans/indian_ocean_8day_{year}.nc")
    outdat2 = outdat.to_dataframe().reset_index()
    outdat2.to_csv(f"data/processed/oceans/indian_ocean_8day_{year}.csv", index=False)
    print(f"Saving: data/processed/oceans/indian_ocean_8day_{year}.csv")

    # Mask Chukchi Beaufort sea Ocean
    geom = chukchi_sea
    outdat = gdf.rio.clip(geom.geometry, all_touched=False, drop=True, invert=False, from_disk=True) 
    outdat.to_netcdf(f"data/processed/oceans/chukchi_sea_8day_{year}.nc")
    print(f"Saving: data/processed/oceans/chukchi_sea_8day_{year}.nc")
    outdat2 = outdat.to_dataframe().reset_index()
    outdat2.to_csv(f"data/processed/oceans/chukchi_sea_8day_{year}.csv", index=False)
    print(f"Saving: data/processed/oceans/chukchi_sea_8day_{year}.csv")

    geom = beaufort_sea
    outdat = gdf.rio.clip(geom.geometry, all_touched=False, drop=True, invert=False, from_disk=True) 
    outdat.to_netcdf(f"data/processed/oceans/beaufort_sea_8day_{year}.nc")
    print(f"Saving: data/processed/oceans/beaufort_sea_8day_{year}.nc")
    outdat2 = outdat.to_dataframe().reset_index()
    outdat2.to_csv(f"data/processed/oceans/beaufort_sea_8day_{year}.csv", index=False)
    print(f"Saving: data/processed/oceans/beaufort_sea_8day_{year}.csv")





file_ = "data/processed/oceans/ca_coast_2022.nc"
file_ = "data/processed/oceans/ca_coast_8day_2017.nc"

def proc_monthly_conf_mats(file_):
    
    # Extract metadata from the filename
    filename = file_.split("/")[-1].split(".")[0]
    region = "_".join(file_.split("/")[-1].split(".")[0].split("_")[0:-1])
    year = file_.split("/")[-1].split(".")[0].split("_")[-1]

    # Load the data and convert to a Pandas DataFrame
    dat = rxr.open_rasterio(file_, decode_times=True)
    indat = dat.to_dataframe().reset_index()

    indat = indat.assign(time = indat['time'].astype(str))
    indat = indat.assign(month = pd.to_datetime(indat['time'], format="%Y-%m-%d %H:%M:%S").dt.month)

    # Extract true and predicted labels
    true_labels = indat['MCLASS'].astype(int)
    pred_labels = indat['VCLASS'].astype(int)

    # Filter out rows where true_labels are zero
    mask = true_labels != 0
    true_labels_filtered = true_labels[mask]
    pred_labels_filtered = pred_labels[mask]

    # create a list of all possible labels
    all_labels = np.unique(np.concatenate((true_labels_filtered, pred_labels_filtered)))

    # create the confusion matrix using the true and predicted labels
    cm = confusion_matrix(true_labels_filtered, pred_labels_filtered, labels=all_labels)

    # Get accuracy metrics
    full_accuracy = accuracy_score(true_labels, pred_labels)
    precision = precision_score(true_labels, pred_labels, average=None, labels=all_labels, zero_division=0)
    recall = recall_score(true_labels, pred_labels, average=None, labels=all_labels, zero_division=0)
    f1 = f1_score(true_labels, pred_labels, average=None, labels=all_labels, zero_division=0)

    # calculate the accuracy for each class
    accuracies = dict()
    for i, label in enumerate(all_labels):
        total = sum(cm[i, :])
        correct = cm[i, i]
        accuracy = correct / total if total > 0 else 0  # Avoid division by zero
        accuracies[label] = accuracy

    df_cm = pd.DataFrame(cm, index=all_labels, columns=all_labels).reset_index()
    df_cm = df_cm.rename(columns={'index': 'seascape'})
    
    df_cm.insert(1, 'region', region)
    df_cm.insert(2, 'year', year)

    df_cm.insert(3, 'acc', accuracies.values())
    df_cm.insert(4, 'full_accuracy', full_accuracy)
    df_cm.insert(5, 'precision', precision)
    df_cm.insert(6, 'recall', recall)
    df_cm.insert(7, 'f1', f1)
    
    df_cm.to_csv(f"data/processed/conf_mats/{filename}.csv", index=False)

    # Print some information about the processing
    print("------------------------------------")
    print(f"Confusion matrix saved to data/processed/conf_mats/{filename}.csv")
    print(f"Region: {region}")
    print(f"Year: {year}")
    print(f"Number of unique labels: {len(all_labels)}")
    print(f"Full accuracy: {full_accuracy:.2f}")
    print("------------------------------------")
    print("Precision by label:")
    print(precision)
    print("Recall by label:")
    print(recall)
    print("F1 score by label:")
    print(f1)


file_ = files[0]

def proc_8day_conf_mats(file_):
    
    # Extract metadata from the filename
    filename = file_.split("/")[-1].split(".")[0]
    region = "_".join(file_.split("/")[-1].split(".")[0].split("_")[0:-1])
    region = region.replace("_8day", "")
    year = file_.split("/")[-1].split(".")[0].split("_")[-1]

    # Load the data and convert to a Pandas DataFrame
    dat = rxr.open_rasterio(file_, decode_times=True)
    indat = dat.to_dataframe().reset_index()

    indat = indat.assign(time = indat['time'].astype(str))
    indat = indat.assign(month = pd.to_datetime(indat['time'], format="%Y-%m-%d %H:%M:%S").dt.month)

    # Extract true and predicted labels
    true_labels = indat['MCLASS'].astype(int)
    pred_labels = indat['VCLASS'].astype(int)

    # Filter out rows where true_labels are zero
    mask = true_labels != 0
    true_labels_filtered = true_labels[mask]
    pred_labels_filtered = pred_labels[mask]

    # create a list of all possible labels
    all_labels = np.unique(np.concatenate((true_labels_filtered, pred_labels_filtered)))

    # create the confusion matrix using the true and predicted labels
    cm = confusion_matrix(true_labels_filtered, pred_labels_filtered, labels=all_labels)

    # Get accuracy metrics
    full_accuracy = accuracy_score(true_labels, pred_labels)
    precision = precision_score(true_labels, pred_labels, average=None, labels=all_labels, zero_division=0)
    recall = recall_score(true_labels, pred_labels, average=None, labels=all_labels, zero_division=0)
    f1 = f1_score(true_labels, pred_labels, average=None, labels=all_labels, zero_division=0)

    # calculate the accuracy for each class
    accuracies = dict()
    for i, label in enumerate(all_labels):
        total = sum(cm[i, :])
        correct = cm[i, i]
        accuracy = correct / total if total > 0 else 0  # Avoid division by zero
        accuracies[label] = accuracy

    df_cm = pd.DataFrame(cm, index=all_labels, columns=all_labels).reset_index()
    df_cm = df_cm.rename(columns={'index': 'seascape'})
    
    df_cm.insert(1, 'region', region)
    df_cm.insert(2, 'year', year)

    df_cm.insert(3, 'acc', accuracies.values())
    df_cm.insert(4, 'full_accuracy', full_accuracy)
    df_cm.insert(5, 'precision', precision)
    df_cm.insert(6, 'recall', recall)
    df_cm.insert(7, 'f1', f1)
    
    df_cm.to_csv(f"data/processed/conf_mats/{filename}.csv", index=False)

    # Print some information about the processing
    print("------------------------------------")
    print(f"Confusion matrix saved to data/processed/conf_mats/{filename}.csv")
    print(f"Region: {region}")
    print(f"Year: {year}")
    print(f"Number of unique labels: {len(all_labels)}")
    print(f"Full accuracy: {full_accuracy:.2f}")
    print("------------------------------------")
    print("Precision by label:")
    print(precision)
    print("Recall by label:")
    print(recall)
    print("F1 score by label:")
    print(f1)





def check_conf_mats(file_):
    try:    
        # Load the data and convert to a Pandas DataFrame
        dat = rxr.open_rasterio(file_, decode_times=True)
    except e as Exeption:
        print(f"Error loading {file_}")




if __name__ == "__main__":

    # download_seascapes()

    years_ = ["2017", "2018", "2019", "2020", "2021", "2022"]

    # proc_seascape_data("2017")

    results = [proc_seascape_data(years) for years in years_]

    results = [proc_monthly_conf_mats(x) for x in glob.glob("data/processed/oceans/*.nc") if x i]

    results = [proc_8day_seascape_data(years) for years in years_]   

    day8_files = glob.glob("data/processed/oceans/*_8day*.nc")

    [check_conf_mats(x) for x in day8_files]

    results = [proc_8day_conf_mats(x) for x in day8_files]









