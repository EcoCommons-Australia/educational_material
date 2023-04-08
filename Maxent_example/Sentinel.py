# Supress Warnings
import warnings
warnings.filterwarnings('ignore')
# Import common GIS tools
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import rasterio.features
# Import Planetary Computer tools
import stackstac
import pystac_client
import planetary_computer
import xrspatial.multispectral as ms

min_lon, min_lat = (138.61,-39.34) # Lower-left corner
max_lon, max_lat = (153.78,-20.60) # Upper-right corner
bbox = (min_lon, min_lat, max_lon, max_lat)

stac = pystac_client.Client.open("https://planetarycomputer.microsoft.com/api/stac/v1")
search = stac.search(
    bbox=bbox,
    datetime="2020-01-01/2020-12-31",
    collections=["sentinel-2-l2a"],
    limit=500, # fetch items in batches of 500
    query={"eo:cloud_cover": {"lt": 20}},
)
items = list(search.get_items())
print('This is the number of scenes that touch our region:',len(items))

signed_items = [planetary_computer.sign(item).to_dict() for item in items]

# Define the pixel resolution for the final product
# Define the scale according to our selected crs, so we will use degrees
resolution = 10 # meters per pixel
scale = resolution / 111320.0 # degrees per pixel for crs=4326

data = (
    stackstac.stack(
        signed_items,
        epsg=4326, # Use common Lat-Lon coordinates
        resolution=scale, # Use degrees for crs=4326
        bounds_latlon = bbox,
#       resampling=rasterio.enums.Resampling.average, # Average resampling method (only required when resolution >10)
        assets=["B04", "B03", "B02", "B08"], # Red, Green, Blue, NIR
        chunksize=4096,
    )
    .where(lambda x: x > 0, other=np.nan) # sentinel-2 uses 0 as nodata
    .assign_coords(band=lambda x: x.common_name.rename("band")) # use common names
)

data = data.persist()

median = data.median(dim="time").compute()

image = ms.true_color(median.sel(band="red"), median.sel(band="green"), median.sel(band="blue"))

image.plot.imshow(figsize=(8,8))
plt.title("Median RGB Mosaic")
plt.axis('off')
plt.show()

ndvi_median = ms.ndvi(median.sel(band="nir"), median.sel(band="red"))

ndvi_median.plot.imshow(cmap="Greens", vmin=0.0, vmax=1.0, figsize=(10,8))
plt.title("Median NDVI")
plt.axis('off')
plt.show()

filename = "S2_mosaic_sample.tiff"


# This will report the pixel dimensions of our mosaic file. Recall that pixel resolution will impact the dimensions.
median.sel(band="red").shape

height = median.sel(band="red").shape[0]
width = median.sel(band="red").shape[1]

with rasterio.open(filename,'w',driver='GTiff',width=width,height=height,count=4,compress='lzw',dtype='float64') as dst:
    dst.write(median.sel(band="red"),1)
    dst.write(median.sel(band="green"),2)
    dst.write(median.sel(band="blue"),3)
    dst.write(median.sel(band="nir"),4)
    dst.close()

# Show the location and size of the new output file#
#!ls *.tiff -lah


# This is an example for a specific Lon-Lat location randomly selected within our sample region.
values = median.sel(x=150.71, y=-33.51, method="nearest").values
print("These are the band values (R,G,B,NIR) for the closest pixel: ", values)


 
