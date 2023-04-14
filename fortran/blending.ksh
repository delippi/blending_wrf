#!/bin/ksh
set -x

dtg0=$1
nDOMAIN=$2
Lx=$3


echo "usage example: ./blending.ksh 2021052700 1 960"

# regional
FG=/mnt/lfs4/BMC/wrfruc/RAPv5_may2021/cycle_test1/2021051223/wrfprd/wrfout_d01_2021-05-13_00_00_00
# global
BG=/mnt/lfs4/BMC/wrfruc/RAPv5_may2021/cycle/boundary/wrfinput_d01_2021-05-13_00_00_00

curdir=`pwd`
BLEND_WRK=$curdir/$dtg0/d${nDOMAIN}/$Lx
LOG=${BLEND_WRK}
mkdir -p ${BLEND_WRK}
BIN=`pwd`

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

#    ln -sf ${FG}/wrfinput_d0${nd} fg  ## CWB deterministic analysis
#    ln -sf ${FG}/${dtg0}/cold_wrfinput_d0${nd}_${dtg0}.relo bg ## NCEP GFS analysis
#    ln -sf ${BG}/wrfinput_d0${nd} bg ## NCEP GFS analysis

    ln -sf ${FG} fg  ## regional deterministic analysis
    ln -sf ${BG} bg ## NCEP GFS analysis
    ${BIN}/da_blending.exe -lx ${Lx} > ${LOG}/blend_main.log.${dtg0} -debug 0
    if [[ `grep "Blending completed successfully" ${LOG}/blend_main.log.${dtg0} | wc -l` -eq 0 ]]; then
       echo "blending fail: please check ${LOG}/blend_main.log.${dtg0}"
       exit 1
    fi
#    cp fg_blend ${FG}/wrfinput_d0${nd}_${dtg0}.blend
    cp fg_blend wrfinput_d0${nd}_${dtg0}.blend

## second step -- diagnose
    ${BIN}/da_diagnose.exe > ${LOG}/blend_diag.log.${dtg0}
    if [[ `grep "Diagnose completed successfully" ${LOG}/blend_diag.log.${dtg0} | wc -l` -eq 0 ]]; then
       echo "blending fail: please check ${LOG}/blend_diag.log.${dtg0}"
       exit 2
    fi
    cp fg_blend_diag wrfinput_d0${nd}_${dtg0}.blend_diag
 #   cp fg_blend_diag ${FG}/wrfinput_d0${nd}_${dtg0}.blend_diag
 #   cp ${FG}/wrfinput_d0${nd}_${dtg0}.blend_diag \
 #      ${FG}/wrfinput_d0${nd}_${dtg0}

   fi
  done
exit 0
