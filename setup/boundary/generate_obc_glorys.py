#!/usr/bin/env python

# this code has been adapted from code created by Andrew Ross - https://github.com/andrew-c-ross/nwa-shared
from glob import glob
#from subprocess import run
import os, sys
import xarray as xr
from boundary import Segment

# xarray gives a lot of unnecessary warnings
import warnings
warnings.filterwarnings('ignore')

def write_year(year, glorys_dirs, segments, is_first_year=False, glorys_file_structure='*.nc'):
    # Open all the glorys files at once?
    glorys = []
    
    for glorys_dir in glorys_dirs:
        glorys.append(
            xr.open_mfdataset(sorted(glob(os.path.join(glorys_dir, glorys_file_structure))))
            .rename({'latitude': 'lat', 'longitude': 'lon', 'depth': 'z'})
        )

    # Temporary assignment to pull two hemispheric halves together
    glorysAll = xr.Dataset()
    # Auto-detect periodicy
    #glorysAll['time'] = glorys[0]['time']
    gloryVariables = ['vo', 'uo', 'so', 'zos', 'thetao']
    for gVar in gloryVariables:
        if gVar == 'zos':
           glorysAll[gVar] = xr.concat([glorys[0][gVar][:,:,:], glorys[1][gVar][:,:,1:-1]], dim='time')
        else:
           glorysAll[gVar] = xr.concat([glorys[0][gVar][:,:,:,:], glorys[1][gVar][:,:,:,1:-1]], dim='time')

    #breakpoint()
    # Floor first time down to midnight so that it matches initial conditions
    if is_first_year:
       tnew = xr.concat((glorysAll['time'][0].dt.floor('1d'), glorysAll['time'][1:]), dim='time')
       glorysAll['time'] = ('time', tnew)

    segCt = 0
    for seg in segments:
        segCt = segCt + 1
        glorysSub = glorysAll
        #breakpoint()
        print("  SEG:", segCt)
        print("  -> regrid_velocity")
        seg.regrid_velocity(glorysSub['uo'], glorysSub['vo'], suffix=year, flood=False)
        for tr in ['thetao', 'so']:
            print("  -> regrid_tracer:",tr)
            seg.regrid_tracer(glorysSub[tr], suffix=year, flood=False)
        print("  -> regrid_tracer: zos")
        seg.regrid_tracer(glorysSub['zos'], suffix=year, flood=False)

# this is an xarray based way to concatenate the obc yearly files into one file (per variable of output)
# the former way to do this was based on ncrcat from NCO tools
def ncrcat_rename(nsegments, ncrcat_outdir, output_dir, delete_old_files=False):
    rename_dict={'thetao':'temp', 'so':'salt', 'zos':'zeta', 'uv':'uv'}
    for var in ['thetao', 'so', 'zos', 'uv']:
        for seg in range(1, nsegments+1):
            comb = xr.open_mfdataset(f'{output_dir}{var}_00{seg}_*')
            if var!='uv':
                comb=comb.rename({f'{var}_segment_00{seg}':f'{rename_dict[var]}_segment_00{seg}'})
                if var!='zos':
                    comb=comb.rename({f'dz_{var}_segment_00{seg}':f'dz_{rename_dict[var]}_segment_00{seg}'})
            # Fix output metadata, including removing all _FillValues.
            all_vars = list(comb.data_vars.keys()) + list(comb.coords.keys())
            encodings = {v: {'_FillValue': None} for v in all_vars}
            encodings['time'].update({'dtype':'float64', 'calendar': 'gregorian'})
            comb.to_netcdf(f'{ncrcat_outdir}{rename_dict[var]}_00{seg}.nc',
                           encoding=encodings,
                           unlimited_dims='time',
                           format='NETCDF3_64BIT')
            print(f'concatenated and renamed {rename_dict[var]}_00{seg}.nc')
            del(comb)
            if delete_old_files==True:
                os.remove(f'{output_dir}{var}_00{seg}_*')

def main():
    first_year = 2010
    glorys_dirs = ['/import/AKWATERS/kshedstrom/glorys', '/import/AKWATERS/kshedstrom/glorys2']
    glorys_file_structure = 'GLORYS_REANALYSIS_{year}-*.nc'
    #output_dir = '/Users/james/Documents/glorys_obc_gen/'
    #ncrcat_outdir = '/Users/james/Documents/glorys_obc_gen/'
    output_dir = '/import/AKWATERS/jrcermakiii/configs/Arctic12/OBC/glorys_obc_gen/'
    ncrcat_outdir = '/import/AKWATERS/jrcermakiii/configs/Arctic12/OBC/glorys_obc_gen/'
    #hgrid = xr.open_dataset('/Users/james/Documents/nwa25/ocean_hgrid.nc')
    hgrid = xr.open_dataset('/import/AKWATERS/jrcermakiii/configs/Arctic12/INPUT2/ocean_hgrid.nc')
    #segments = [
    #    Segment(1, 'south', hgrid, output_dir=output_dir),
    #    Segment(2, 'north', hgrid, output_dir=output_dir),
    #    Segment(3, 'east', hgrid, output_dir=output_dir)
    #    ]

    # Open boundary segments
    segments = [
        Segment(1, 'north', hgrid, output_dir=output_dir),
        Segment(2, 'west', hgrid, output_dir=output_dir),
        Segment(3, 'south', hgrid, output_dir=output_dir),
        Segment(4, 'east', hgrid, output_dir=output_dir)
    ]

    for y in range(2010,2011):
        print(y)
        # Just do one DAY!
        glorys_file_structure = f'GLORYS_REANALYSIS_{y}-01-01.nc'
        write_year(y, glorys_dirs, segments, is_first_year= (y==first_year), glorys_file_structure=glorys_file_structure)
    
    ncrcat_rename(len(segments), ncrcat_outdir, output_dir)


if __name__ == '__main__':
    main()

sys.exit()
# Post processing of time
import xarray as xr
filelist = ['uv_001.nc', 'uv_003.nc' , 'zeta_002.nc']
for f in filelist:
    print(f)
    comb = xr.open_dataset(f)
    # Fix output metadata, including removing all _FillValues.
    all_vars = list(comb.data_vars.keys()) + list(comb.coords.keys())
    encodings = {v: {'_FillValue': None} for v in all_vars}
    encodings['time'].update({'dtype':'float64', 'calendar': 'gregorian'})
    comb.to_netcdf(f,
                   encoding=encodings,
                   unlimited_dims='time',
                   format='NETCDF3_64BIT')

