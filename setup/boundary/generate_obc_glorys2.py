#!/usr/bin/env python

# this code has been adapted from code created by Andrew Ross - https://github.com/andrew-c-ross/nwa-shared
from glob import glob
import os, sys
import xarray as xr
from boundary import Segment
import datetime
import argparse

# xarray gives a lot of unnecessary warnings
import warnings
warnings.filterwarnings('ignore')

def write_day(year, glorys_dirs, segments, is_first_record=False, glorys_file_structure='*.nc', timeRef=None, flooding=False):
    glorysFiles = []

    for glorys_dir in glorys_dirs:
        '''
        glorys.append(
            xr.open_mfdataset(sorted(glob(os.path.join(glorys_dir, glorys_file_structure))))
            .rename({'latitude': 'lat', 'longitude': 'lon', 'depth': 'z'})
        )
        '''
        [glorysFiles.append(fn) for fn in sorted(glob(os.path.join(glorys_dir, glorys_file_structure)))]

    glorysSub = xr.Dataset()

    for fn in glorysFiles:
        glorysSub = xr.merge([glorysSub, xr.open_mfdataset(fn)])

    # Rename
    glorysSub = glorysSub.rename({'latitude': 'lat', 'longitude': 'lon', 'depth': 'z'})

    # Debug
    #breakpoint()

    # Floor first time down to midnight so that it matches initial conditions
    if is_first_record:
       tnew = xr.concat((glorysSub['time'][0].dt.floor('1d'), glorysSub['time'][1:]), dim='time')
       glorysSub['time'] = ('time', tnew.data)
       print("First record!")
    else:
       print("Not first record...")

    # Debug
    #breakpoint()

    segCt = 0
    for seg in segments:
        segCt = segCt + 1
        #glorysSub = glorysAll
        #breakpoint()
        print("  SEG:", segCt)
        print("  -> regrid_velocity")
        seg.regrid_velocity(glorysSub['uo'], glorysSub['vo'], suffix=year, flood=flooding)
        for tr in ['thetao', 'so']:
            print("  -> regrid_tracer:",tr)
            seg.regrid_tracer(glorysSub[tr], suffix=year, flood=flooding)
        print("  -> regrid_tracer: zos")
        seg.regrid_tracer(glorysSub['zos'], suffix=year, flood=flooding)

# this is an xarray based way to concatenate the obc yearly files into one file (per variable of output)
# the former way to do this was based on ncrcat from NCO tools
def ncrcat_rename(nsegments, ncrcat_outdir, output_dir, delete_old_files=False, fdate=None):
    rename_dict={'thetao':'temp', 'so':'salt', 'zos':'zeta', 'uv':'uv'}
    for var in ['thetao', 'so', 'zos', 'uv']:
        for seg in range(1, nsegments+1):
            fileList = glob(f'{output_dir}/{var}_00{seg}_????.nc')
            #print(var,seg,fileList)
            comb = xr.open_mfdataset(fileList)
            if var!='uv':
                comb=comb.rename({f'{var}_segment_00{seg}':f'{rename_dict[var]}_segment_00{seg}'})
                comb=comb.rename({f'dz_{var}_segment_00{seg}':f'dz_{rename_dict[var]}_segment_00{seg}'})

            # Fix output metadata, including removing all _FillValues.
            all_vars = list(comb.data_vars.keys()) + list(comb.coords.keys())
            encodings = {v: {'_FillValue': None} for v in all_vars}
            encodings['time'].update({'dtype':'float64', 'calendar': 'gregorian'})
            comb.to_netcdf(f'{ncrcat_outdir}/{rename_dict[var]}_00{seg}_{fdate}.nc',
                           encoding=encodings,
                           unlimited_dims='time',
                           format='NETCDF3_64BIT')
            print(f'concatenated and renamed {rename_dict[var]}_00{seg}.nc')
            del(comb)
            if delete_old_files==True:
                os.remove(f'{output_dir}/{var}_00{seg}_*')

def getArguments():
    '''
    Process command line arguments.

    --date="YYYY-MM-DD"
    --firstRec=True
    --doYear=YYYY
    --procDays=n
    --procDir=<processing subdirectory>
    '''

    parser = argparse.ArgumentParser(description='Process daily Glorys OBCs')
    parser.add_argument("--flooding", help='Flag for flooding missing grid points over land', type=bool, default=False)
    parser.add_argument("--firstRec", help='First record of the timeseries', type=bool, default=False)
    parser.add_argument("--date", help='Date of record in YYYY-MM-DD format; or start of period', type=str, default=None)
    parser.add_argument("--doYear", help='Process an entire year of Glorys', type=int, default=None)
    parser.add_argument("--procDays", help='Only process n days of Glorys', type=int, default=None)
    parser.add_argument("--timeRef", help='Set time reference for processing', type=str, default=None)
    parser.add_argument("--startDate", help='Set starting date for processing', type=str, default=None)
    parser.add_argument("--procDir", help='Subdir for processing', type=str, default=None)

    args = parser.parse_args()
    return args

def collectFilesByYear(doYear, glorys_dirs):

    dataDir = glorys_dirs[0]
    yearOfFiles = sorted(glob(os.path.join(dataDir, f'GLORYS_REANALYSIS_{doYear}-*.nc')))

    return yearOfFiles

def main():
    args = getArguments()

    #glorys_dirs = ['/home/cermak/workdir/datasets/global/glorys', '/home/cermak/workdir/datasets/global/glorys2']
    glorys_dirs = ['/import/AKWATERS/kshedstrom/glorys', '/import/AKWATERS/kshedstrom/glorys2']
    glorys_file_structure = 'GLORYS_REANALYSIS_{year}-*.nc'
    output_dir = os.path.join('/home/jrcermakiii/workdir/configs/Arctic12/OBC/obc_gen', args.procDir)
    ncrcat_outdir = os.path.join('/home/jrcermakiii/workdir/configs/Arctic12/OBC/obc_gen', args.procDir)
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

    fdate = args.date
    try:
        trec = datetime.datetime.strptime(fdate,"%Y-%m-%d")
    except:
        sys.exit("Date format for --date is incorrect, use YYYY-MM-DD.")

    is_first_record = args.firstRec
    doYear = args.doYear
    procDays = args.procDays
    timeRef = args.timeRef
    startDate = args.startDate
    floodFlag = args.flooding
    skipDates = False
    if startDate:
        skipDates = True

    # Determine year from fdate
    y = trec.year

    # Collect files by specified year
    if doYear:
        daysToProcess = collectFilesByYear(doYear, glorys_dirs)
    else:
        daysToProcess = [f'GLORYS_REANALYSIS_{fdate}.nc']

    nProc = 0

    for glorys_file_structure in daysToProcess:

        if nProc > 0:
            is_first_record = False

        #glorys_file_structure = f'GLORYS_REANALYSIS_{fdate}.nc'
        glorys_file_structure = os.path.basename(glorys_file_structure)
        fdate = os.path.basename(glorys_file_structure)[-13:-3]
        if skipDates:
            if fdate != startDate:
                continue
            skipDates = False
            is_first_record = False
            print("Skipped to",fdate)

        print("*****>", fdate)
        write_day(y, glorys_dirs, segments, is_first_record=is_first_record,
            glorys_file_structure=glorys_file_structure, timeRef=timeRef, flooding=floodFlag)

        ncrcat_rename(len(segments), ncrcat_outdir, output_dir, fdate=fdate)

        nProc = nProc + 1

        if procDays:
            if procDays == nProc:
                break

if __name__ == '__main__':
    main()
