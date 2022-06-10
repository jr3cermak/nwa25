#!/usr/bin/env python

import os, sys, datetime, glob
import argparse

# GLOBALS

parser = None

def main():
  args = getArguments()

  processOBCS(args)


def getArguments():
    '''
    Check and/or create chinook jobs for a given year of OBCs
    --doYear=YYYY
    --maxDays=9
    --procDir=<processing subdirectory>
    --jobDir=<directory for new slurm scripts>
    --templateScript=<location of template script for slurm jobs>
    --showJobs=<bool> Show jobs to be created
    --createJobs=<bool> Actually create the jobs
    --report=<bool> Show data report

    If procDir is `obc_gen`, then a sequence of directories will be created under procDir
    by year.  IE: obc_gen/pYYYY_01 obc_gen/pYYYY_02 ....

    If neither showJobs or createJobs is True, this script will only check for missing days.
    '''

    global parser

    parser = argparse.ArgumentParser(description='Check and/or create Glorys OBCs processing jobs')
    parser.add_argument("--showJobs", help='Show the jobs to be created', type=bool, default=False)
    parser.add_argument("--createJobs", help='Actually create the job scripts', type=bool, default=False)
    parser.add_argument("--report", help='Show data report', type=bool, default=False)
    parser.add_argument("--procDir", help='Processing subdirectory (obc_gen)', type=str, default='obc_gen')
    parser.add_argument("--jobDir", help='New job script directory (jobs)', type=str, default='jobs')
    parser.add_argument("--templateScript", help='Template script for new slurm jobs (template.slurm)', type=str, default='template.slurm')
    parser.add_argument("--doYear", help='The selected year of Glorys to process', type=int, default=None)
    parser.add_argument("--maxDays", help='Maximum number of days to allow job to run', type=int, default=9)

    args = parser.parse_args()
    return args


def showHelp(msg):

    global parser

    if msg:
        print(msg)


def getYear(yr):

    if yr is None:
        showHelp("Please specify a year to examine or process.")
        sys.exit()

    mo = 1
    dy = 1
    dtbeg = datetime.date(yr, mo, dy)
    dtend = datetime.date(yr+1, mo, dy)

    return (dtbeg, dtend)


def writeJob(args, startDate, numDays, jobSeq):

    #print(jobSeq, startDate, numDays)

    jobName = "proc%d_%02d.slurm" % (args.doYear, jobSeq)
    destFilename = os.path.join(args.jobDir, jobName)

    fpi = open(args.templateScript, 'r')
    fpo = open(destFilename, 'w')

    fpo.write(fpi.read())
    fpi.close()

    # Write environment variables
    fpo.write("\n# Environment variables\n")
    fpo.write("YR=%d\n" % (args.doYear))
    fpo.write("JDIR=p${YR}_%02d\n" % (jobSeq))
    fpo.write("PDIR=%s\n" % (args.procDir))
    fpo.write("\n")

    # Directory setup
    fpo.write("cd /home/jrcermakiii/workdir/configs/Arctic12/OBC/${PDIR}\n")
    # If job directory exist, do some clean up
    fpo.write("if [ -d \"${JDIR}\" ]; then\n")
    fpo.write("    rm ${JDIR}/*_${YR}.nc\n")
    fpo.write("    rm ${JDIR}/regrid_segment_00*.nc\n")
    fpo.write("else\n")
    # Just create the directory
    fpo.write("    mkdir ${JDIR}\n")
    fpo.write("fi\n")
    fpo.write("\n")

    # Actually start the job
    fpo.write("cd /home/jrcermakiii/workdir/configs/Arctic12/OBC\n")
    fpo.write("python generate_obc_glorys2.py --flooding=True --date=${YR}-01-01 --firstRec=True --doYear=${YR} --procDir=${JDIR} --startDate=%s --procDays %d\n" % (startDate.strftime("%Y-%m-%d"), numDays))

    fpo.close()


def createJobs(args, missingDates, maxDays=9):

    jobSeq = 0
    ct = 0
    ctMax = maxDays
    prevDate = None
    startDate = None

    for dt in missingDates:
      if ct == 0:
        startDate = dt
      else:
        deltaDt = dt - prevDate

        if deltaDt.days > 1:
          jobSeq = jobSeq + 1
          writeJob(args, startDate, ct, jobSeq)
          startDate = dt
          ct = 0

      if ct == ctMax:
        jobSeq = jobSeq + 1
        writeJob(args, startDate, ct, jobSeq)
        startDate = dt
        ct = 0

      ct = ct + 1
      
      prevDate = dt

    if ct > 0:
      jobSeq = jobSeq + 1
      writeJob(args, startDate, ct, jobSeq)

def processOBCS(args):

    yr = None
    missingDates = []

    if args.doYear:
        yr = args.doYear

    (dtbeg, dtend) = getYear(yr)
    dt = dtbeg

    while dt < dtend:
      dts = dt.strftime("%Y-%m-%d")

      # There should be 16 files with 4 bounaries * 4 variables
      fglobname = os.path.join(
          'obc_gen',
          'p%d' % (yr),
          '*_%s.nc' % (dts)
      )

      flist = glob.glob(fglobname)
      dateok = False
      if flist:
          if len(flist) == 16:
              dateok = True

      if not(dateok):
          #print("Missing", dts)
          missingDates.append(dt)
  
      dt = dt + datetime.timedelta(days=1)

    if len(missingDates) == 0:
      print("Entire year exists.")
      sys.exit()

    if args.report:
      print("Missing dates", missingDates)

    if args.createJobs:
      createJobs(args, missingDates, args.maxDays)

if __name__ == '__main__':
    main()
