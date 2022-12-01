#!/bin/bash
#BSUB -J %J_NAME%          # Name of the job.
#BSUB -o %J_OUT%  # Appends std output to file %J.out.
#BSUB -e %J_ERR%  # Appends std error to file %J.err.
#BSUB -q %J_QUEUE%
#BSUB -n %J_CPUS%
#BSUB -P %J_PROJ% # Project number
#
#
# by AC Goglio (CMCC)
# annachiara.goglio@cmcc.it
#
# Written: 11/02/2021
#
# Script for map TS extraction 
# Ini file: map_extr.ini 
#
set -u
set -e
#set -x 
#################################################

# Source ini file
  SRC_DIR=%SRC_DIR%
  source ${SRC_DIR}/map_extr.ini

# Read and check infos (work dir, file names, archive dir, etc.)

  # Workdir check
  if [[ -d $ANA_WORKDIR ]]; then
     cd $ANA_WORKDIR
     echo "WORKDIR: $ANA_WORKDIR"
     
     # Clean workdir
     echo "WARNING: I am going to remove all files in $ANA_WORKDIR ..."
     sleep 1
     for TO_BE_RM in $( ls $ANA_WORKDIR ); do
         rm $ANA_WORKDIR/$TO_BE_RM
         echo $TO_BE_RM
     done
     # Cp the ini file to the workdir 
     cp ${SRC_DIR}/map_extr.ini ${ANA_WORKDIR}

  else
     echo "ERROR: WORKDIR $ANA_WORKDIR NOT FOUND!!"
     exit
  fi

  # Input file check
  
  # Num check
  if [[ ${#ANA_INPATHS[@]} -ne ${#ANA_INFILES_TPL[@]} ]]; then
     echo "ERROR: Check Inputs array, something is missing!"
     exit
  else
     INSET_NUM=${#ANA_INPATHS[@]}
     echo "The num of input file sets is $INSET_NUM"
  fi

  # Input File check and link

  IDX_IN=0
  while [[ $IDX_IN -lt $INSET_NUM ]]; do

    if [[ -d ${ANA_INPATHS[${IDX_IN}]} ]]; then

      IDX_DATE=$ANA_STARTDATE
      
      while [[ $IDX_DATE -le $ANA_ENDDATE ]]; do

        echo "Date: $IDX_DATE"
        ANA_INFILE=$( echo ${ANA_INFILES_TPL[${IDX_IN}]} | sed -e s/%YYYYMMDD%/$IDX_DATE/g  )

        FOUND_NUM=0
        if [[ $GRID_TO_EXTRACT != "uv2t" ]]; then
           #if [[ -e ${ANA_INPATHS[$IDX_IN]}/${IDX_DATE:0:6}/$ANA_INFILE ]]; then
              FOUND_NUM=$(( $FOUND_NUM + 1 ))
              echo "Found infile: $ANA_INFILE"
              ln -sf ${ANA_INPATHS[$IDX_IN]}/${IDX_DATE:0:6}/$ANA_INFILE .
           #else
           #   echo "NOT Found infile: $ANA_INFILE in path: ${ANA_INPATHS[$IDX_IN]}/${IDX_DATE:0:6}/"
           #fi
        else
           if [[ -e ${ANA_INPATHS[$IDX_IN]}/$ANA_INFILE ]]; then
              FOUND_NUM=$(( $FOUND_NUM + 1 ))
              echo "Found infile: $ANA_INFILE"
              ln -sf ${ANA_INPATHS[$IDX_IN]}/$ANA_INFILE .
           fi
        fi

      IDX_DATE=$( date -d "$IDX_DATE 1 day" +%Y%m%d ) 
      done    
 
    else 
      echo "ERROR: Input dir ${ANA_INPATHS[${IDX_IN}]} NOT FOUND!!"
      exit
    fi

    IDX_IN=$(( $IDX_IN + 1 ))
    done

# Read ana type from ini file and set the environment
  
  if [[ $MAP_FLAG == 1 ]]; then
    echo "Map analisys required!"
    module load $MAP_MODULE
    echo "Loading the Enviroment.."
  fi


####################### EXTRACTION ##############################

# Extract values from netCDF

 if [[ $MAP_FLAG == 1 ]]; then

  echo "------ Maps Extraction from model outputs ------"

  # Loop on model datasets
  IDX_IN=0
  while [[ $IDX_IN -lt $INSET_NUM ]]; do

    echo "I am working on ${ANA_INTAG[$IDX_IN]} dataset..."

    ############# 3D VARS #########################

    echo "I am working on 3D vars..."

    # Loop on 3D FIELDS
    for VAR in ${VAR3D_NAME[@]}; do
     echo "Extracting field: $VAR .. "

     # Input file names (linked in work dir)
     EXT_INS=$(  echo "${ANA_INFILES_TPL[${IDX_IN}]}" | sed -e "s/%YYYYMMDD%/"*"/g" )

       VLEV="allv"

       MAP3D_OUTFILE=$( echo "$MAP3D_OUTFILE_TPL" | sed -e "s/%LEV%/$VLEV/g" -e "s/%FIELD%/${VAR}/g" -e "s/%ANA_STARTDATE%/${ANA_STARTDATE}/g" -e "s/%ANA_ENDDATE%/${ANA_ENDDATE}/g" -e "s/%INDATASET%/${ANA_INTAG[$IDX_IN]}/g" )

      echo "# MAP extraction from ${ANA_INTAG[$IDX_IN]}" 
      echo "# VAR: ${VAR}" 
      echo "# LEV: ${VLEV}" 
      echo "# Period: ${ANA_STARTDATE}-${ANA_ENDDATE}"

      # Extraction

      # Loop on input daily files
      IDX_NC=0
      for MAPEXT in $EXT_INS; do
          echo "Infile: $MAPEXT"

          # Select or compute var
          if [[ $GRID_TO_EXTRACT == "uv2t" ]]; then
             echo "Computing i=sqrt(uo*uo+vo*vo)..."
             #cdo expr,'i=sqrt(uo*uo+vo*vo);ut=uo;vt=vo' $MAPEXT name_tmp_${IDX_NC}.nc
             cdo expr,'ut=uo;vt=vo' $MAPEXT name_tmp_${IDX_NC}.nc
          else
             echo "Selecting VAR:$VAR..." 
             cdo selname,$VAR $MAPEXT name_tmp_${IDX_NC}.nc
          fi

          # Select vertical level
          #cdo intlevel,${VLEV} name_tmp.nc lev_tmp_${IDX_NC}.nc  
          #rm name_tmp.nc       


      IDX_NC=$(( $IDX_NC + 1 ))
      done

      # Cat extracted files
      echo "Cat extracted files.."
      cdo mergetime name_tmp_*.nc $MAP3D_OUTFILE #Instead of cdo cat
      rm name_tmp_*.nc
      # ADD time order 
    done

    ############# 2D VARS #########################

    echo "I am working on 2D vars..."

    # Loop on 2D FIELDS
    for VAR in ${VAR2D_NAME[@]}; do
     echo "Extracting field: $VAR .. "

     # Input file names (linked in work dir)
     EXT_INS=$(  echo "${ANA_INFILES_TPL[${IDX_IN}]}" | sed -e "s/%YYYYMMDD%/"*"/g" )

       VLEV=0
       echo "Depth: $VLEV .. "

       MAP2D_OUTFILE=$( echo "$MAP2D_OUTFILE_TPL" | sed -e "s/%LEV%/$VLEV/g" -e "s/%FIELD%/${VAR}/g" -e "s/%ANA_STARTDATE%/${ANA_STARTDATE}/g" -e "s/%ANA_ENDDATE%/${ANA_ENDDATE}/g" -e "s/%INDATASET%/${ANA_INTAG[$IDX_IN]}/g" )

       AMPPHA_NEWFILE=$( echo "$AMPPHA_NEWFILE_TPL" | sed -e "s/%LEV%/$VLEV/g" -e "s/%FIELD%/${VAR}/g" -e "s/%ANA_STARTDATE%/${ANA_STARTDATE}/g" -e "s/%ANA_ENDDATE%/${ANA_ENDDATE}/g" -e "s/%INDATASET%/${ANA_INTAG[$IDX_IN]}/g" )

      echo "# MAP extraction from ${ANA_INTAG[$IDX_IN]}" 
      echo "# VAR: ${VAR}" 
      echo "# LEV: ${VLEV}" 
      echo "# Period: ${ANA_STARTDATE}-${ANA_ENDDATE}"

      # Extraction

      # Loop on input daily files
      IDX_NC=0
      for MAPEXT in $EXT_INS; do
          echo "Infile: $MAPEXT"

          # Select var
          cdo selname,$VAR $MAPEXT name_tmp_${IDX_NC}.nc

      IDX_NC=$(( $IDX_NC + 1 ))
      done

      # Cat extracted files and create the new file with the same grid and time vars
      echo "Cat extracted files.."
      cdo mergetime name_tmp_*.nc $MAP2D_OUTFILE 
      cp name_tmp_0.nc ${AMPPHA_NEWFILE}_ini
      rm name_tmp_*.nc
  
    done

   IDX_IN=$(( $IDX_IN + 1 ))
   done
   
   # Clean links to the mod
   for TOBERM in $(ls ${EXT_INS}); do
       rm -v $TOBERM
   done

 fi

