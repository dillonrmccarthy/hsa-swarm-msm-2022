#!/usr/bin/bash
#need to launch the dispdev vmd command from a SEPARATE bash file because you need the first one for the /dev/null output to supress vmd, because the nohub command doesnt work from python, so python needs to launch a shell script, which in turn, launches another shell script.

export VMD=/opt/vmd/vmd3
cd $dir_path
$VMD/vmd3 -dispdev text -e $script -args $scriptname
sleep 1s
exit
