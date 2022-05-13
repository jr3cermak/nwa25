#!/usr/bin/env python
# coding: utf-8

# # Subset ERA5 from Global grid to NWA25 domain and calculate Specific Huimdity & Total Rain Rate


# slice down the data
import xarray as xr
import os
import cftime
import numpy as np
from glob import glob
import os

# Functions for humidity borrowed and adapted from MetPy.calc: https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.html
def mixing_ratio(partial_press, total_press, molecular_weight_ratio=0.622):
    return (molecular_weight_ratio * partial_press
                / (total_press - partial_press))


def specific_humidity_from_mixing_ratio(mr):
    return mr / (1 + mr)


def saturation_vapor_pressure(temperature):
    sat_pressure_0c = 6.112e2 # Pa
    return sat_pressure_0c * np.exp(17.67 * (temperature - 273.15) # K -> C
                                        / (temperature - 29.65))   # K -> C

def saturation_mixing_ratio(total_press, temperature):
    return mixing_ratio(saturation_vapor_pressure(temperature), total_press)

def make_periodic_append_longitude(ds, dvar, delta_long):
    # This appends a longitude point for the given variable

    times = ds['time']
    latitudes = ds['latitude']
    longitudes = ds['longitude']

    # extend the longitude coordinate by 1
    datasets = []
    datasets.append(longitudes)
    # add 0.25 degrees to the final longitude value because the resolution of the ERA5 is 0.25 degrees
    datasets.append(longitudes[-1] + delta_long)
    longitudes = xr.concat(datasets, dim='longitude')

    # Form the new data structure with a periodic longitude
    ds_ext = xr.Dataset({
        dvar : xr.DataArray(
            data   = np.zeros((len(range(0,ds.dims['time'])), len(range(0,ds.dims['latitude'])), len(range(0,ds.dims['longitude']+1)))),
            dims   = ['time', 'latitude', 'longitude'],
            coords = {'time': times, 'latitude' : latitudes, 'longitude' : longitudes},
            attrs  = {
                'units'     : ds[dvar].attrs['units']
            })
        })

    ds_ext[dvar].values[:,:,0:len(ds['longitude'])] = ds[dvar].values
    ds_ext[dvar].values[:,:,len(ds['longitude']):len(ds['longitude']) + 1] = ds[dvar][:,:,0:1].values

    return ds_ext

def save_attrs(ds):

    attr_array = {}
    for dvar in list(ds.variables):
        attr_array[dvar] = ds[dvar].attrs

    return attr_array

def fix_encoding_attrs(ds, oldattrs):

    for dvar in list(ds.variables):
        ds[dvar].encoding.update({'dtype':'float64', '_FillValue': None})
        if 'missing_value' in ds[dvar].encoding: ds[dvar].encoding.pop('missing_value')
        if 'scale_factor' in ds[dvar].encoding: ds[dvar].encoding.pop('scale_factor')
        if 'add_offset' in ds[dvar].encoding: ds[dvar].encoding.pop('add_offset')
        ds[dvar].attrs = oldattrs[dvar]

    return ds


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

# subset
years = range(1996,1998)
latsub = slice(90,39)
lonsub = slice(0,360)

# storage
# NOT USED
#subdir2 = "/Volumes/P8/workdir/james/ERA5/nwa25/subset/"
era5dir = "/import/AKWATERS/kshedstrom/ERA5/"
subdir = "/import/AKWATERS/jrcermakiii/data/ERA5_periodic_subset/"

# Tasks
#  1. Subset first
#  2. Make periodic

for f in era5_dict.keys():
    print(f)
    for y in years:
        print("-> %d" % (y))

        if f == 'ERA5_total_rain_rate':
            crr = xr.open_dataset(str(era5dir + 'ERA5_convective_rain_rate_' + str(y) + '.nc'))
            breakpoint()
            crr = xr.open_dataset(str(era5dir + 'ERA5_convective_rain_rate_' + str(y) + '.nc')).sel(latitude=latsub, longitude=lonsub)
            lsrr = xr.open_dataset(str(era5dir + 'ERA5_large_scale_rain_rate_' + str(y) + '.nc')).sel(latitude=latsub, longitude=lonsub)
            trr = crr.drop('crr')
            trr['trr'] = crr['crr'] + lsrr['lsrr']
            trr['trr'].attrs = {'units': 'kg m-2 s-1','long_name': 'Total Rainfall Rate'}
            trr['trr'].encoding = {k: v for k, v in crr.crr.encoding.items() if k in {'_FillValue', 'missing_value', 'dtype'}}
            trr.to_netcdf(str(subdir + f + '_' + str(y) + ".nc"), mode='w', format='NETCDF4_CLASSIC')
            crr.close()
            lsrr.close()
            trr.close()

        if f == 'ERA5_2m_specific_humidity':
            pair = xr.open_dataset(str(era5dir + 'ERA5_surface_pressure_' + str(y) + '.nc'))['sp'].sel(latitude=latsub, longitude=lonsub) # Pa
            tdew = xr.open_dataset(str(era5dir + 'ERA5_2m_dewpoint_temperature_' + str(y) + '.nc'))['d2m'].sel(latitude=latsub, longitude=lonsub) # K

            smr = saturation_mixing_ratio(pair, tdew)
            sphum = specific_humidity_from_mixing_ratio(smr)

            sphum.name = 'huss'
            sphum = sphum.to_dataset()

            # Remove all _FillValue
            all_vars = list(sphum.data_vars.keys()) + list(sphum.coords.keys())
            encodings = {v: {'_FillValue': None} for v in all_vars}

            # Also fix the time encoding
            encodings['time'].update({'dtype':'float64', 'calendar': 'gregorian', 'units': 'hours since 1900-01-01 00:00:00'})
            
            fout=str(subdir + f + '_' + str(y) + ".nc")
            sphum.to_netcdf(
                fout,
                format='NETCDF4_CLASSIC',
                engine='netcdf4',
                encoding=encodings,
                unlimited_dims=['time']
            )
            sphum.close()
            
        if 'total_rain_rate' not in f and 'specific_humidity' not in f:
            ds = xr.open_dataset(str(era5dir + f + '_' + str(y) + ".nc")).sel(latitude=latsub, longitude=lonsub)
            ds_attrs = save_attrs(ds)
            ds = make_periodic_append_longitude(ds, era5_dict[f], 0.25)
            ds = fix_encoding_attrs(ds, ds_attrs)
            ds.to_netcdf(str(subdir + f + '_' + str(y) + ".nc"),format="NETCDF4_CLASSIC")
            ds.close()
            breakpoint()

