#!/bin/ksh --login

# Set the queueing options
#SBATCH --ntasks=4   ## refer to the run_fcst task section in the workflow xml
#SBATCH -t 6:40:00
#SBATCH -A nrtrr
# #SBATCH -qos debug
#SBATCH -J blending 
##SBATCH --partition=sjet  ## check which jet partition was used from the run_fcst log file in the cycle that needs to be repeated
set -x

dtg0=$1
nDOMAIN=$2
wlen=$3

dtg0=2021051300
nDOMAIN=1 
wlen=1200

#echo "usage example: ./blending.ksh 2021052700 1 1200"

dtgS=$4
MPI=$5
NCU=$6

#WRFDTGCDF=/mnt/lfs4/BMC/wrfruc/Chunhua.Zhou/HRRR/hrrr_databasedir_rap_gfs.$wlen/run/$dtg0/wrfprd_bc
#WRFDTGCDF2=/mnt/lfs4/BMC/wrfruc/Chunhua.Zhou/HRRR/hrrr_databasedir_gfs/run/$dtg0/wrfprd_bc
## WRFDTGCDF3=/mnt/lfs4/BMC/wrfruc/Chunhua.Zhou/HRRR/hrrr_databasedir_gfs/run/$dtg0/wrfprd_bc

# regional 
WRFDTGCDF=/mnt/lfs4/BMC/wrfruc/RAPv5_may2021/cycle_test1/2021051223/wrfprd/wrfout_d01_2021-05-13_00_00_00
# global
WRFDTGCDF2=/mnt/lfs4/BMC/wrfruc/RAPv5_may2021/cycle/boundary/wrfinput_d01_2021-05-13_00_00_00

curdir=`pwd`
BLEND_WRK=$curdir/$dtg0/d${nDOMAIN}/$wlen
WRFLOG=${BLEND_WRK}
mkdir -p ${BLEND_WRK}
WRFBIN=`pwd`

case ${nDOMAIN} in 
1)
  DOMAIN="1" ;;
2)
  DOMAIN="1 2" ;;
3)
  DOMAIN="1 2 3" ;;
esac

cd ${BLEND_WRK}

 for nd in ${DOMAIN}
  do
   if [ ${nd} = 1 ] || [ ${nd} = 2 ]; then

## first step -- blending

#    ln -sf ${WRFDTGCDF}/wrfinput_d0${nd} fg  ## CWB deterministic analysis
#    ln -sf ${WRFDTGCDF}/${dtg0}/cold_wrfinput_d0${nd}_${dtg0}.relo bg ## NCEP GFS analysis 
#    ln -sf ${WRFDTGCDF2}/wrfinput_d0${nd} bg ## NCEP GFS analysis 

    ln -sf ${WRFDTGCDF} fg  ## regional deterministic analysis
    ln -sf ${WRFDTGCDF2} bg ## NCEP GFS analysis 
srun    ${WRFBIN}/da_blending.exe -lx ${wlen} > ${WRFLOG}/blend_main.log.${dtg0} -debug 1
    if [[ `grep "Blending completed successfully" ${WRFLOG}/blend_main.log.${dtg0} | wc -l` -eq 0 ]]; then
       echo "blending fail: please check ${WRFLOG}/blend_main.log.${dtg0}"
       exit 1
    fi
#    cp fg_blend ${WRFDTGCDF}/wrfinput_d0${nd}_${dtg0}.blend
    cp fg_blend wrfinput_d0${nd}_${dtg0}.blend

## second step -- diagnose
srun    ${WRFBIN}/da_diagnose.exe > ${WRFLOG}/blend_diag.log.${dtg0}
    if [[ `grep "Diagnose completed successfully" ${WRFLOG}/blend_diag.log.${dtg0} | wc -l` -eq 0 ]]; then
       echo "blending fail: please check ${WRFLOG}/blend_diag.log.${dtg0}"
       exit 2
    fi
    cp fg_blend_diag wrfinput_d0${nd}_${dtg0}.blend_diag
 #   cp fg_blend_diag ${WRFDTGCDF}/wrfinput_d0${nd}_${dtg0}.blend_diag
 #   cp ${WRFDTGCDF}/wrfinput_d0${nd}_${dtg0}.blend_diag \
 #      ${WRFDTGCDF}/wrfinput_d0${nd}_${dtg0}

   fi
  done
exit 0
