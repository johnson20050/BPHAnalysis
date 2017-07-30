
if [ "$1" == "rm" ]; then
    for name in `xrdfs se01.grid.nchc.org.tw ls  $NCHCpath | grep ".root" `
    do
        xrdfs se01.grid.nchc.org.tw rm $name
    done
elif [ "$1" == "rmdir" ]; then
    for name in `xrdfs se01.grid.nchc.org.tw ls  $NCHCpath | grep ".root" `
    do
        xrdfs se01.grid.nchc.org.tw rmdir $name
    done
fi
