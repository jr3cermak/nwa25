#!/usr/bin/env python
# coding: utf-8

# # Pad ERA5 Forcing File using Xarray

# In[ ]:


import xarray as xr
import os
import cftime
import numpy as np

datadir='/glade/scratch/jsimkins/ERA5/original/'
outdir='/glade/scratch/jsimkins/ERA5/padded/'
year=1996
era5_dict = {#'ERA5_sea_ice_cover':'siconc',
            #'ERA5_10m_u_component_of_wind':'u10',
            #'ERA5_sea_surface_temperature':'sst',
            #'ERA5_10m_v_component_of_wind':'v10',
            #'ERA5_2m_temperature':'t2m',
            #'ERA5_surface_solar_radiation_downwards':'ssrd',
            #'ERA5_surface_thermal_radiation_downwards':'strd',
            #'ERA5_total_rain_rate':'trr',
            #'ERA5_mean_sea_level_pressure':'msl',
            'ERA5_2m_specific_humidity':'huss'}

for f in era5_dict.keys():
    print(f)
    # open the file for current year
    current = xr.open_dataset(f"{datadir}/{f}_{year}.nc")
    previous = xr.open_dataset(f"{datadir}/{f}_{year-1}.nc").isel(time=-1)
    next_data = xr.open_dataset(f"{datadir}/{f}_{year+1}.nc").isel(time=0)
    out = xr.concat([previous, current, next_data], dim="time")
    current.close()
    previous.close()
    next_data.close()
    all_vars = list(out.data_vars.keys()) + list(out.coords.keys())
    encodings = {v: {'_FillValue': 1.0e20} for v in all_vars}
    encodings['time'].update({'dtype':'float64', 'calendar': 'gregorian'})
    out['time'].attrs['long_name'] = 'time'
    out['time'].attrs['standard_name'] = 'time'
    out['time'].attrs['axis'] = 'T'
    out['latitude'].attrs['long_name'] = 'Latitude'
    out['latitude'].attrs['units'] = 'degrees_north'
    out['latitude'].attrs['axis'] = 'Y'
    out['longitude'].attrs['long_name'] = 'Longitude'
    out['longitude'].attrs['units'] = 'degrees_east'
    out['longitude'].attrs['axis'] = 'X'
    out=out.transpose("time", "latitude", "longitude")
    # latitude needs to be reindexed for some reason
    out=out.reindex(latitude=list(reversed(out.latitude)))
    out.to_netcdf(f'{outdir}{f}_{year}.nc', format="NETCDF4_CLASSIC", encoding=encodings, unlimited_dims='time')
    out.close()


# # Pad ERA5 Forcing file using Climate Data Operator (CDO) Mergetime

# In[ ]:


# using cdo
import xarray as xr
import os
import cftime
import numpy as np

datadir='/glade/scratch/jsimkins/ERA5/original/'
outdir='/glade/scratch/jsimkins/ERA5/padded/'
year=1996
era5_dict = {'ERA5_sea_ice_cover':'siconc',
            'ERA5_10m_u_component_of_wind':'u10',
            'ERA5_sea_surface_temperature':'sst',
            'ERA5_10m_v_component_of_wind':'v10',
            'ERA5_2m_temperature':'t2m',
            'ERA5_surface_solar_radiation_downwards':'ssrd',
            'ERA5_surface_thermal_radiation_downwards':'strd',
            'ERA5_total_rain_rate':'trr',
            'ERA5_mean_sea_level_pressure':'msl',
            'ERA5_2m_specific_humidity':'huss'}

for f in era5_dict.keys():
    # grab previous 2 timestamps and next 2 timestamps. We need 2 timestamps so time is stored as a dimension for each variable. This enables CDO mergetime to work properly.
    previous = xr.open_dataset(f"{datadir}/{f}_{year-1}.nc", decode_times=False).isel(time=slice(-3,-1))
    next_data = xr.open_dataset(f"{datadir}/{f}_{year+1}.nc", decode_times=False).isel(time=slice(0,2))

    previous.to_netcdf(f'{outdir}/prev.nc', format="NETCDF4",unlimited_dims='time')
    next_data.to_netcdf(f'{outdir}/next.nc', format="NETCDF4",unlimited_dims='time')
    previous.close()
    next_data.close()
    
    # make a system call to use CDO mergetime to concatenate our data
    os.system(str("cdo -b F32 -f nc4c mergetime " + outdir + "prev.nc " + datadir + f + "_" + str(year) + ".nc " + outdir + "next.nc " + outdir + f))  


# In[ ]:




