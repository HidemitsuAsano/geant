#!/bin/tcsh -f
set Version="17"
set OUTDIR="/group/had/knucl/e15/asano/sim/"

set logbasedir="/home/had/hiasano/logs/"
set date=`date +%Y%m%d_%H%M`
set logdir="${logbasedir}${date}_sim"
echo "log files  ${logdir}"

if( ! -d $logdir) then 
  mkdir -p  $logdir
endif 

set OUTDIRSUB="${OUTDIR}sim${Version}"
if( ! -d $OUTDIRSUB) then 
 mkdir -p $OUTDIRSUB
endif


@ i = 0
while ($i < 400)   

  set EXEC___="./knucl"
  set CONF___="conf/Run78/analyzer_kwsk_sim_DoraAir.conf"
  set CARD___="KnuclSetting_nSmpip_DoraAir.card"
  set MAC___="run.mac"
  set jobnum=`printf  "%03d"  $i`

  #set INPFILE=${DATADIR}"run78_0${jobnum}.dat"
  set OUTFILE=${OUTDIRSUB}"/sim_nSmpip_0${jobnum}.root"
  #set CDSFILE=${KWSKDIR}"run78_0${jobnum}_evtracking.root"

  #echo ${INPFILE}
  echo ${OUTFILE}
  #echo ${CDSFILE}
  @ j = $i + 400 
  #j =+400
  set logname = "${logdir}/run$i.log"
  bsub -o $logname -q l ${EXEC___} ${CARD___} ${MAC___} ${OUTFILE} $j
    @ i ++
end

