#!/bin/tcsh -f
set Version="1"
set OUTDIR="/group/had/knucl/e15/asano/sim/"

set logbasedir="/home/had/hiasano/logs/"
set date=`date +%Y%m%d_%H%M`
set logdir="${logbasedir}${date}_simSmpip_ns"
echo "log files  ${logdir}"

if( ! -d $logdir) then 
  mkdir -p  $logdir
endif 

set OUTDIRSUB="${OUTDIR}simSmpip_ns${Version}"
if( ! -d $OUTDIRSUB) then 
 mkdir -p $OUTDIRSUB
endif


@ i = 0
while ($i < 400)   

  set EXEC___="./knucl"
  set CONF___="conf/Run78/analyzer_kwsk_sim_DoraAir.conf"
  set CARD___="KnuclSetting_Smpip_ns_DoraAir.card"
  set MAC___="run.mac"
  set jobnum=`printf  "%03d"  $i`

  #set INPFILE=${DATADIR}"run78_0${jobnum}.dat"
  set OUTFILE=${OUTDIRSUB}"/sim_Smpip_ns_0${jobnum}.root"
  #set CDSFILE=${KWSKDIR}"run78_0${jobnum}_evtracking.root"

  #echo ${INPFILE}
  echo ${OUTFILE}
  #echo ${CDSFILE}
  set logname = "${logdir}/run$i.log"
  bsub -o $logname -q l ${EXEC___} ${CARD___} ${MAC___} ${OUTFILE} ${jobnum}
    @ i ++
end

