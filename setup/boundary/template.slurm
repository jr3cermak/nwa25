#!/bin/bash
#SBATCH --partition=t1small
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=48:00:00
#SBATCH --constraint=[Broadwell]

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/import/AKWATERS/jrcermakiii/local/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/import/AKWATERS/jrcermakiii/local/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/import/AKWATERS/jrcermakiii/local/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/import/AKWATERS/jrcermakiii/local/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

ulimit -u 25000
ulimit -n 2048

conda activate gridTools

# Job entry should look like
#YR=2003
#JDIR=p${YR}_000
#PDIR=obc_gen

#cd /home/jrcermakiii/workdir/configs/Arctic12/OBC/${PDIR}
#mkdir ${JDIR}
#cd ${JDIR}

#python generate_obc_glorys2.py --flooding=True --date=${YR}-01-01 --firstRec=True --doYear=${YR} --procDir=${JDIR} --startDate=${YR}-09-04 --procDays 7
