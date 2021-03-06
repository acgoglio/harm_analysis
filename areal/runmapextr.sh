#!/bin/bash
#
# by AC Goglio (CMCC)
# annachiara.goglio@cmcc.it
#
# Written: 01/02/2021
#
#set -u
set -e
#set -x 
########################
echo "*********** Data extraction for areal hanalysis *********"
MEXTR_INIFILE='map_extr.ini'
JOB_TEMPLATE='mextrjob.temp'
JOB_TORUN='mextr.job'
SRC_DIR=$(pwd)

# Check and load ini file
if [[ -e ./${MEXTR_INIFILE} ]]; then
   echo "Reading ini file ${MEXTR_INIFILE}.."
   source ./${MEXTR_INIFILE}
   echo "..Done"
else
   echo "${MEXTR_INIFILE} file NOT FOUND in $(pwd)! Why?"
   exit 
fi

# Check job template file
if [[ -e  ./${JOB_TEMPLATE} ]]; then
   echo "Found the job template ${JOB_TEMPLATE}"
else
   echo "${JOB_TEMPLATE} file NOT FOUND in $(pwd)! Why?"
   exit 
fi

# Check work directory
if [[ -d ${ANA_WORKDIR} ]]; then
   echo "Work dir: ${ANA_WORKDIR}"
   cp ${MEXTR_INIFILE} ${ANA_WORKDIR}/
else
   echo "Work dir: ${ANA_WORKDIR} NOT FOUND! Why?"
   exit
fi

# Built the job from the template
echo "I am building the job.."
# Sed file creation and sobstitution of parematers in the templates  
SED_FILE=sed_file.txt
cat << EOF > ${ANA_WORKDIR}/${SED_FILE}
   s/%J_NAME%/${J_NAME//\//\\/}/g
   s/%J_OUT%/${J_OUT//\//\\/}/g
   s/%J_ERR%/${J_ERR//\//\\/}/g
   s/%J_QUEUE%/${J_QUEUE//\//\\/}/g
   s/%J_CWD%/${J_CWD//\//\\/}/g
   s/%J_CPUS%/${J_CPUS//\//\\/}/g
   s/%J_PROJ%/${J_PROJ//\//\\/}/g
   #
   s/%SRC_DIR%/${SRC_DIR//\//\\/}/g
EOF

      sed -f ${ANA_WORKDIR}/${SED_FILE} ${JOB_TEMPLATE} > ${ANA_WORKDIR}/${JOB_TORUN}
      rm ${ANA_WORKDIR}/${SED_FILE}
echo ".. Done"
echo "Job path/name: ${ANA_WORKDIR}/${JOB_TORUN}"

# Run the job
echo "Submitting job ${J_NAME} to queue ${J_QUEUE} (Good luck!).."
bsub<${ANA_WORKDIR}/${JOB_TORUN}
echo "Check the output in ${ANA_WORKDIR} and/or the errors in ${J_ERR}!"





