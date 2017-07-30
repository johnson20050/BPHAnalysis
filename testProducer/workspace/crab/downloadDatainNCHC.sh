#!/usr/bin/env sh


if [ "$1" == "" ]; then
    echo "need to input link in NCHC"
    exit
fi

# use date to be the folder name
if [ "`date | cut -d" " -f3`" != "" ]; then
    MON=` date | cut -d" " -f2`
    DAY=` date | cut -d" " -f3`
    YEAR=`date | cut -d" " -f6`
else
    MON=` date | cut -d" " -f2`
    DAY=` date | cut -d" " -f4`
    YEAR=`date | cut -d" " -f7`
fi
folderName=CRABdata_${DAY}_${MON}_${YEAR}

mkdir -p ~/Data/CRABdata/$folderName
for name in `cat $1 | grep '.root'`
do
    xrdcp root://se01.grid.nchc.org.tw/${name} ~/Data/CRABdata/${folderName}/
done
