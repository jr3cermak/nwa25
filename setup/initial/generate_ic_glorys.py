#!/usr/bin/env python

import numpy as np
import xarray
import xesmf

# https://github.com/raphaeldussin/HCtFlood
import os, sys, glob
from HCtFlood import kara as flood
from boundary import rotate_uv
import numpy as np
import xarray


def vgrid_to_interfaces(vgrid, max_depth=6500.0):
    if isinstance(vgrid, xarray.DataArray):
        vgrid = vgrid.data
    zi = np.concatenate([[0], np.cumsum(vgrid)])
    zi[-1] = max_depth
    return zi


def vgrid_to_layers(vgrid, max_depth=6500.0):
    if isinstance(vgrid, xarray.DataArray):
        vgrid = vgrid.data
    ints = vgrid_to_interfaces(vgrid, max_depth=max_depth)
    z = (ints + np.roll(ints, shift=1)) / 2
    layers = z[1:]
    return layers



def write_initial(glorys_dirs, glorys_file, start_date, vgrid_file, grid_file, output_file):
    vgrid = xarray.open_dataarray(vgrid_file)
    z = vgrid_to_layers(vgrid)
    ztarget = xarray.DataArray(
        z,
        name='zl',
        dims=['zl'],
        coords={'zl': z},
    )

    glorys_files = []
    for glorys_dir in glorys_dirs:
        fileGlob = sorted(glob.glob(os.path.join(glorys_dir, glorys_file)))
        if len(fileGlob) > 0:
            glorys_files = glorys_files + fileGlob
    glorys = (
        xarray.open_mfdataset(glorys_files)
        [['thetao', 'so', 'zos', 'uo', 'vo']]
        .rename({'longitude': 'lon', 'latitude': 'lat'})
    )
    #glorys = (
    #    xarray.open_dataset(glorys_file)
    #    [['thetao', 'so', 'zos', 'uo', 'vo']]
    #    .rename({'longitude': 'lon', 'latitude': 'lat'})
    #)

    # Round time down to midnight
    glorys['time'] = (('time', ), glorys['time'].dt.floor('1d'))

    # Interpolate GLORYS vertically onto target grid.
    # Depths below bottom of GLORYS are filled by extrapolating the deepest available value.
    revert = glorys.interp(depth=ztarget, kwargs={'fill_value': 'extrapolate'}).ffill('zl', limit=None)

    # Flood temperature and salinity over land.
    flooded = xarray.merge((
        flood.flood_kara(revert[v], zdim='zl') for v in ['thetao', 'so', 'uo', 'vo']
    ))

    # flood zos separately to avoid the extra z=0 added by flood_kara.
    flooded['zos'] = flood.flood_kara(revert['zos']).isel(z=0).drop('z')

    # Horizontally interpolate the vertically interpolated and flooded data onto the MOM grid.
    target_grid = xarray.open_dataset(grid_file)
    target_grid['x'] -= 360.0
    target_t = (
        target_grid
        [['x', 'y']]
        .isel(nxp=slice(1, None, 2), nyp=slice(1, None, 2))
        .rename({'y': 'lat', 'x': 'lon', 'nxp': 'xh', 'nyp': 'yh'})
    )
    # Interpolate u and v onto supergrid to make rotation possible
    target_uv = (
        target_grid
        [['x', 'y']]
        .rename({'y': 'lat', 'x': 'lon'})
    )

    regrid_kws = dict(method='bilinear', reuse_weights=False, periodic=False, ignore_degenerate=True)

    glorys_to_t = xesmf.Regridder(glorys, target_t, filename='regrid_glorys_tracers.nc', **regrid_kws)
    glorys_to_uv = xesmf.Regridder(glorys, target_uv, filename='regrid_glorys_uv.nc', **regrid_kws)

    interped_t = glorys_to_t(flooded[['thetao', 'so', 'zos']]).drop(['lon', 'lat'])

    # Interpolate u and v, rotate, then extract individual u and v points
    interped_uv = glorys_to_uv(flooded[['uo', 'vo']]).drop(['lon', 'lat'])
    urot, vrot = rotate_uv(interped_uv['uo'], interped_uv['vo'], target_grid['angle_dx'])
    uo = urot.isel(nxp=slice(0, None, 2), nyp=slice(1, None, 2)).rename({'nxp': 'xq', 'nyp': 'yh'})
    uo.name = 'uo'
    vo = vrot.isel(nxp=slice(1, None, 2), nyp=slice(0, None, 2)).rename({'nxp': 'xh', 'nyp': 'yq'})
    vo.name = 'vo'

    interped = (
        xarray.merge((interped_t, uo, vo))
        .transpose('time', 'zl', 'yh', 'yq', 'xh', 'xq')
    )

    # Rename to match MOM expectations.
    interped = interped.rename({
        'thetao': 'temp',
        'so': 'salt',
        'zos': 'ssh',
        'uo': 'u',
        'vo': 'v'
    })
    #replace time with the start date of the model
    interped['time'] = (('time', ), [start_date])
    # Fix output metadata, including removing all _FillValues.
    all_vars = list(interped.data_vars.keys()) + list(interped.coords.keys())
    # below does the same as NCO: ncatted -a _FillValue,v10,o,f,100000000000000000000 ERA5_10m_v_component_of_wind_1996.nc ERA5_10m_v_component_of_wind_1996.nc
    encodings = {v: {'_FillValue': 1.0e20} for v in all_vars}
    encodings['time'].update({'dtype':'float64', 'calendar': 'gregorian'})
    interped['zl'].attrs = {
        'long_name': 'Layer pseudo-depth, -z*',
         'units': 'meter',
         'cartesian_axis': 'Z',
         'positive': 'down'
    }

    interped.to_netcdf(
        output_file,
        format='NETCDF3_64BIT',
        engine='netcdf4',
        encoding=encodings,
        unlimited_dims='time'
    )


def main():
    #glorys_dir = '/glade/scratch/jsimkins/'
    glorys_dirs = ['/import/AKWATERS/kshedstrom/glorys', '/import/AKWATERS/kshedstrom/glorys2']
    glorys_year = 1993
    glorys_file = f'GLORYS_REANALYSIS_{glorys_year}-01-01.nc'
    start_date = np.datetime64(f'{glorys_year}-01-01T03:00:00')
    vgrid_file = '/import/AKWATERS/jrcermakiii/configs/Arctic12/INPUT2/vgrid_75_2m.nc'
    #vgrid_file = '/glade/work/jsimkins/gridInfo/nwa25/vgrid_75_2m.nc'
    grid_file = '/import/AKWATERS/jrcermakiii/configs/Arctic12/INPUT2/ocean_hgrid.nc'
    #grid_file = '/glade/work/jsimkins/gridInfo/nwa25/ocean_hgrid.nc'
    output_file = f'glorys_ic_75z_{glorys_year}.nc'
    write_initial(glorys_dirs, glorys_file, start_date, vgrid_file, grid_file, output_file)


if __name__ == '__main__':
    main()
