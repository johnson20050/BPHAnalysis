#!/usr/bin/env sh
# load the file named "path" in the current directory and read the path in NCHC.
# if there is no "path", show home directory in NCHC.
# You can iterate this script to enter the further folder.


checkPathFile=`ls path`
if [ "$checkPathFile" != "" ]; then
    NCHCpath=`cat path`
    rm .path
    echo -e "\n" 
    #NCHCpath=/cms/store/user/ltsai/testFile_161216/Charmonium
    for path in `xrdfs se01.grid.nchc.org.tw ls  $NCHCpath`
    do 
        # store the path in "crabPath"
        if [ "`echo $path | grep -v '.root' `" != "" ]; then
            echo $path >> .path
            echo $path
        else
            file=$path
            echo $file >> .path
            echo $file
        fi
    done
    echo -e "\n"
    if [ "`cat .path | grep -v '.root'`" != "" ]; then
        echo 'Modify path again'
    else
        echo 'Here is your file'
    fi
    mv .path path
else
    echo -e "\n"
    
    for _path in `xrdfs se01.grid.nchc.org.tw ls /cms/store/user/ltsai `
    do
        echo $_path >> path
        echo $_path
    done
    echo -e "\n"
    echo 'List of home directory in NCHC is stored in path'
    echo 'Please modify the file and run this script again'
    echo $path
        
    
fi

