#!/bin/tcsh -f
set Version="4"
#set DATADIR="/group/had/knucl/e15/data/Run78/"
set OUTDIR="/group/had/knucl/e15/asano/sim/"
#set KWSKDIR="/group/had/knucl/e15/shinngo/Run78/evtracking/"

set logbasedir="/home/had/hiasano/logs/"
set date=`date +%Y%m%d_%H%M`
set logdir="${logbasedir}${date}_simnpipiL"
echo "log files  ${logdir}"

if( ! -d $logdir) then 
  mkdir -p  $logdir
endif 

set OUTDIRSUB="${OUTDIR}simnpipiL${Version}"
if( ! -d $OUTDIRSUB) then 
 mkdir -p $OUTDIRSUB
endif


@ i = 0
while ($i < 400)   

  set EXEC___="./knucl"
#  set CONF___="conf/Run78/analyzer_kwsk_sim_DoraAir.conf"
  set CONF___="conf/Run78/analyzer_kwsk_sim.conf"
#  set CARD___="KnuclSetting_npipiL_DoraAir.card"
  set CARD___="KnuclSetting_npipiL.card"
  set MAC___="run.mac"
  set jobnum=`printf  "%03d"  $i`

  #set INPFILE=${DATADIR}"run78_0${jobnum}.dat"
  set OUTFILE=${OUTDIRSUB}"/sim_npipiL_0${jobnum}.root"
  #set CDSFILE=${KWSKDIR}"run78_0${jobnum}_evtracking.root"

  #echo ${INPFILE}
  echo ${OUTFILE}
  #echo ${CDSFILE}
  @ j = $i + 0
  set logname = "${logdir}/run$i.log"
  bsub -o $logname -q l ${EXEC___} ${CARD___} ${MAC___} ${OUTFILE} $j
    @ i ++
end

