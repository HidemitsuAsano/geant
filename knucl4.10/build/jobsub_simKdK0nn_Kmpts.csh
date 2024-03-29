#!/bin/tcsh -f
set Version="1"
#set DATADIR="/group/had/knucl/e15/data/Run78/"
set OUTDIR="/group/had/knucl/e15/asano/sim/"
#set KWSKDIR="/group/had/knucl/e15/shinngo/Run78/evtracking/"

set logbasedir="/home/had/hiasano/logs/"
set date=`date +%Y%m%d_%H%M`
set logdir="${logbasedir}${date}_simK0nn_Kmpts"
echo "log files  ${logdir}"

if( ! -d $logdir) then 
  mkdir -p  $logdir
endif 

set OUTDIRSUB="${OUTDIR}simK0nn_Kmpts${Version}"
if( ! -d $OUTDIRSUB) then 
 mkdir -p $OUTDIRSUB
endif


@ i = 0
while ($i < 400)   

  set EXEC___="./knucl"
  set CONF___="conf/Run78/analyzer_kwsk_sim.conf"
  set CARD___="KnuclSetting_K0nn_Kmpts.card"
  set MAC___="run.mac"
  set jobnum=`printf  "%03d"  $i`

  set OUTFILE=${OUTDIRSUB}"/sim_K0nn_Kmpts_0${jobnum}.root"

  #echo ${INPFILE}
  echo ${OUTFILE}
  #echo ${CDSFILE}

  set logname = "${logdir}/run$i.log"
  bsub -o $logname -q l ${EXEC___} ${CARD___} ${MAC___} ${OUTFILE} ${jobnum}
    @ i ++
end
