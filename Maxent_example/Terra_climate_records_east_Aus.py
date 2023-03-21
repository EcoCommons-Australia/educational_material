#EY data challenge weather data

import warnings
warnings.filterwarnings('ignore')
# Import common GIS tools
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import rasterio.features
import folium
import math
# Import Planetary Computer tools
import pystac_client
import planetary_computer

min_lon, min_lat = (148.18, -38.41) # Lower-left corner (longitude, latitude)
max_lon, max_lat = (154.44, -16.38) # Upper-right corner (longitude, latitude)

bbox = (min_lon, min_lat, max_lon, max_lat)
latitude = (min_lat, max_lat)
longitude = (min_lon, max_lon)

def _degree_to_zoom_level(l1, l2, margin = 0.0):
    degree = abs(l1 - l2) * (1 + margin)
    zoom_level_int = 0
    if degree != 0:
        zoom_level_float = math.log(360/degree)/math.log(2)
        zoom_level_int = int(zoom_level_float)
    else:
        zoom_level_int = 18
    return zoom_level_int

def display_map(latitude = None, longitude = None):
    margin = -0.5
    zoom_bias = 0
    lat_zoom_level = _degree_to_zoom_level(margin = margin, *latitude ) + zoom_bias
    lon_zoom_level = _degree_to_zoom_level(margin = margin, *longitude) + zoom_bias
    zoom_level = min(lat_zoom_level, lon_zoom_level)
    center = [np.mean(latitude), np.mean(longitude)]

    map_hybrid = folium.Map(location=center,zoom_start=zoom_level,
        tiles=" http://mt1.google.com/vt/lyrs=y&z={z}&x={x}&y={y}",attr="Google")
    line_segments = [(latitude[0],longitude[0]),(latitude[0],longitude[1]),
                     (latitude[1],longitude[1]),(latitude[1],longitude[0]),
                     (latitude[0],longitude[0])]
    map_hybrid.add_child(folium.features.PolyLine(locations=line_segments,color='red',opacity=0.8))
    map_hybrid.add_child(folium.features.LatLngPopup())

    return map_hybrid

# Plot bounding box on a map
f = folium.Figure(width=600, height=600)
m = display_map(latitude,longitude)
f.add_child(m)

import pystac
collection = pystac.read_file("https://planetarycomputer.microsoft.com/api/stac/v1/collections/terraclimate")
asset = collection.assets["zarr-https"]

import fsspec
import xarray as xr
store = fsspec.get_mapper(asset.href)
data = xr.open_zarr(store, **asset.extra_fields["xarray:open_kwargs"])

 # View the dimensions, coordinates and variables of the dataset
data

clipped_data = data.sel(lon=slice(min_lon,max_lon),lat=slice(max_lat,min_lat),time=slice('2015-01-01','2019-12-31'))

parsed_data = clipped_data[['tmax', 'tmin', 'ppt', 'soil']]

 # View the dimensions, coordinates and variables of the dataset
parsed_data

import bottleneck
import netCDF4
import dask

dataDIR = '../Terraclim_EY_E_Aus.nc'
parsed_data.to_netcdf(dataDIR)


