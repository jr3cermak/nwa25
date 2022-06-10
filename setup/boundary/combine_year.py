#!/usr/bin/env python

import xarray as xr
import numpy as np
import datetime
import os, sys, glob
import argparse

def getArguments():
    '''
    Process command line arguments.

    --srcDir=Source directory
    --dstDir=Destination directory
    --doYear=YYYY
    '''

    parser = argparse.ArgumentParser(description='Process daily Glorys OBCs')
    parser.add_argument("--srcDir", help='Source directory', type=str, default=None)
    parser.add_argument("--dstDir", help='Destination directory', type=str, default=None)
    parser.add_argument("--doYear", help='Year of data to combine', type=int, default=None)
    parser.add_argument("--fixTime", help='Ensure hour of first record is midnight', type=bool, default=False)
    parser.add_argument("--pad", help='Copy first record to prior records', type=bool, default=False)
    parser.add_argument("--padDir", help='Directory of files requiring record padding', type=str, default=None)

    args = parser.parse_args()
    return args

def main():

    args = getArguments()
    dstDir = args.dstDir
    padDir = args.padDir
    yr = args.doYear

    # Collect files to merge using xr.open_mfdataset()
    varNames = ['salt', 'temp', 'uv', 'zeta']

    for seg in range(1,5):
        print("->",seg)
        for vName in varNames:
            destFile = f'{dstDir}/{vName}_00{seg}.nc'
            padFile = None
            doMerge = False
            
            if not(os.path.isfile(destFile)):
                doMerge = True

            if padDir:
                padFile = f'{padDir}/{vName}_00{seg}.nc'
            if not(os.path.isfile(padFile)):
                print("ERROR: Prior year file to pad was not found:", padFile)
                sys.exit()

            if doMerge:
                print("  ->",vName,"merge")

                fileNames = sorted(glob.glob(os.path.join(args.srcDir,f'{vName}_00{seg}_{yr}*')))
                nfiles = len(fileNames)
                if nfiles < 365 or nfiles > 366:
                    print("Inconsistent number of files found:",nfiles)
                    sys.exit()
                mergedData = xr.open_mfdataset(fileNames)
                #breakpoint()
                # Temporary one off
                if vName == 'zeta':
                    vrbs = list(mergedData.keys())
                    for vrb in vrbs:
                        if vrb.find("_zos_") >= 0:
                            vrbOrig = vrb
                            vrbNew = vrb.replace("_zos_","_zeta_")
                            mergedData = mergedData.rename({vrbOrig: vrbNew})
                encodings = {}
                for vrbName in list(mergedData.variables):
                    encodings[vrbName] = {'_FillValue': None}
                encodings['time'] = {'_FillValue': None}
                encodings['time'].update({'dtype':'float64', 'calendar': 'gregorian'})

                # If fixTime is True, make sure the time coordinate of the first record is
                # midnight.  Use the special .dt accessor to get to date and time components.
                # Updating values used as coordinates is harder than plain variables.
                if args.fixTime:
                    begOfYear = xr.DataArray(datetime.datetime(yr, 1, 1))
                    firstRec = mergedData['time'][0]
                    currentHour = float(firstRec.dt.hour)
                    if currentHour != 0.0:
                        diffSec = int((firstRec - begOfYear).dt.seconds)
                        if diffSec > (60 * 60 * 24):
                            print("ERROR: Adjustment period too large?", begOfYear.data, firstRec.data)
                            sys.exit()
                        print("    - Adjusting first record time coordinate by", diffSec, "seconds.")
                        md = mergedData['time'].data
                        md[0] = (firstRec - np.timedelta64(diffSec, 's')).data
                        mdUpd = mergedData.assign_coords({'time':('time',md,mergedData['time'].attrs)})
                        mergedData = mdUpd

                # Reset time valid_min and valid_max attributes
                #breakpoint()
                mergedData['time'].attrs['valid_min'] = mergedData['time'].min().dt.strftime('%Y-%m-%d %H:%M:%S').to_dict()['data']
                mergedData['time'].attrs['valid_max'] = mergedData['time'].max().dt.strftime('%Y-%m-%d %H:%M:%S').to_dict()['data']

                mergedData.to_netcdf(destFile, encoding=encodings,
                    unlimited_dims='time', format='NETCDF3_64BIT')
                #breakpoint()
            else:
                mergedData = xr.open_dataset(destFile)

            # if padding is True, the first time record of the recently written record is
            # appended to the end of the prior year.
            if args.pad:
                doPadding = False
                # Check the prior year to see if the record has already been appended
                priorYear = xr.open_dataset(padFile)

                # If time of these records do not match, set doPadding True
                doPadding = bool(priorYear['time'][-1] != mergedData['time'][0])

                if doPadding:
                    print("  ->",vName,"padding")

                    # Create empty padded dataset
                    padPriorYear = xr.Dataset()
                    
                    # Collect variable names and coordinate names
                    varNamesPad = list(set(list(priorYear.variables)) - set(list(priorYear.coords)))
                    crdNamesPad = list(priorYear.coords)

                    # Loop through variables and merge items through concatenation
                    # priorYear + mergedData[0]
                    for varNamePad in varNamesPad:
                        firstRec = mergedData[varNamePad][0]
                        priorRec = priorYear[varNamePad]
                        padPriorYear[varNamePad] = xr.concat([priorRec, firstRec], "time")

                    padPriorYear['time'].attrs['valid_max'] = padPriorYear['time'].max().dt.strftime('%Y-%m-%d %H:%M:%S').to_dict()['data']
                    encodings = {}
                    for vrbName in list(padPriorYear.variables):
                        encodings[vrbName] = {'_FillValue': None}

                    # Release resources
                    priorYear.close()

                    # Overwrite the original file
                    padPriorYear.to_netcdf(padFile, encoding=encodings,
                        unlimited_dims='time', format='NETCDF3_64BIT')
                else:
                    priorYear.close()

            # Release resources
            mergedData.close()

if __name__ == '__main__':
    main()
