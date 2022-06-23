#!/bin/bash

skip_vmd=$4

RESET='\e[0m'
BL_IT='\033[3;36m'
BL_IT_WBG='\033[3;47;36m'
BL_IT_B='\033[1;31;36m'

echo "${BL_IT_B}----------BEGINNING CV ANALYSIS----------${RESET}"
echo " "
export VMD=/opt/vmd/vmd3
export script=$1 #this is the path of the tcl script for VMD
export scriptname="generate_cvs.tcl"
export main_path=$2 #THIS IS THE MAIN PATH, WHERE THE ORIGINAL PYTHON 'main.py' was submitted from
export script_path=$3

#-------------------------------------------------------------------------------
#replica 00
export dir_path=$main_path'/rep_0/*'

if [ "$skip_vmd" = "yes" ]; then :
else
    nohup $script_path/submit_vmd_sub.sh >/dev/null 2>&1 &&
    : #need the colon becaue && cant directly evalute "fi", so i say && do nothing (:)
fi
echo "${BL_IT}replica 0 complete, moving to replica 1${RESET}"
sleep 1s

#-------------------------------------------------------------------------------
#replica 01
export dir_path=$main_path'/rep_1/*'
if [ "$skip_vmd" = "yes" ]; then :
else
    nohup $script_path/submit_vmd_sub.sh >/dev/null 2>&1 &&
    : #need the colon becaue && cant directly evalute "fi", so i say && do nothing (:)
fi
echo "${BL_IT}replica 1 complete, moving to replica 2${RESET}"
sleep 1s

#-------------------------------------------------------------------------------
#replica 02
export dir_path=$main_path'/rep_2/*'
if [ "$skip_vmd" = "yes" ]; then :
else
    nohup $script_path/submit_vmd_sub.sh >/dev/null 2>&1 &&
    : #need the colon becaue && cant directly evalute "fi", so i say && do nothing (:)
fi
echo "${BL_IT}replica 2 complete, moving to replica 3${RESET}"
sleep 1s

#-------------------------------------------------------------------------------
#replica 03
export dir_path=$main_path'/rep_3/*'
if [ "$skip_vmd" = "yes" ]; then :
else
    nohup $script_path/submit_vmd_sub.sh >/dev/null 2>&1 &&
    : #need the colon becaue && cant directly evalute "fi", so i say && do nothing (:)
fi
echo "${BL_IT}replica 3 complete, moving to replica 4${RESET}"
sleep 1s

#-------------------------------------------------------------------------------
#replica 04
export dir_path=$main_path'/rep_4/*'
if [ "$skip_vmd" = "yes" ]; then :
else
    nohup $script_path/submit_vmd_sub.sh >/dev/null 2>&1 &&
    : #need the colon becaue && cant directly evalute "fi", so i say && do nothing (:)
fi
echo "${BL_IT}replica 4 complete, creating indicator file${RESET}"
sleep 1s

#-------------------------------------------------------------------------------
#now, need to make a file that indicates we are all done :)
#python will need to find this file before it can move forward
if [ "$skip_vmd" = "yes" ]; then
    :
else
    touch .VMD_COMPLETE.INDC
fi
echo " "
echo "${BL_IT_B}----------CV generation complete----------${RESET}"
exit

#one & will submit to background and continue, two && will make sure the prior command finishes before moving forward
