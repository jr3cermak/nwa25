# Stock command

python chkYear.py --procDir obc_gen --jobDir jobs --templateScript template.slurm\
  --doYear 1994\
  --showJobs True --createJobs True --maxDays 19

# Post processing

## First year

python ../combine_year.py --doYear 1993 --srcDir p1993 --dstDir 1993 --fixTime True

## Subsequent years

python ../combine_year.py --doYear 1994 --srcDir p1994 --dstDir 1994 --pad True --padDir 1993
