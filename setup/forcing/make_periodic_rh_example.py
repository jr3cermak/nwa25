#!/usr/bin/env python

# Code from James Simkins via Slack#arctic-mom6/20220418

import xarray as xr

era5 = xr.open_dataset("/Users/james/Downloads/ERA5_2m_specific_humidity_1996.nc")
f = 'huss'

times = era5['time']
latitudes = era5['latitude']
longitudes = era5['longitude']

# extend the longitude coordinate by 1
datasets = []
datasets.append(longitudes)
# add 0.25 degrees to the final longitude value because the resolution of the ERA5 is 0.25 degrees
datasets.append(longitudes[-1] + 0.25)
longitudes = xr.concat(datasets, dim='longitude')

era5_ext = xr.Dataset({
    f : xr.DataArray(
        data   = np.zeros((len(range(0,era5.dims['time'])), len(range(0,era5.dims['latitude'])), len(range(0,era5.dims['longitude']+1)))),   # enter data here
        dims   = ['time', 'latitude', 'longitude'],
        coords = {'time': times, 'latitude' : latitudes, 'longitude' : longitudes},
        attrs  = {
            'units'     : era5[f].attrs['units']
        }
       )
})

era5_ext[f].values[:,:,0:len(era5['longitude'])] = era5[f].values
era5_ext[f].values[:,:,len(era5['longitude']):len(era5['longitude']) + 1] = era5[f][:,:,0:1].values
