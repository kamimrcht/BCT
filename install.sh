#!/bin/bash



#THIS INSTALLATION REQUIRES GCC, GIT, MAKE, CMAKE



function help {
echo "Bct installation script"
echo "This installation requires GCC>=4.9, GIT, MAKE and CMAKE3"
echo "-f absolute path of folder to put the binaries"
echo "-t to use multiple thread for compilation (default 8)"
}


#~ mkdir src;
#~ cd src;

threadNumber=8
folder=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
folder+="/bin"


while getopts "hf:t:" opt; do
case $opt in
h)
help
exit
;;
f)
echo "use folder: $OPTARG" >&2
folder=$OPTARG
;;
t)
echo "use  $OPTARG threads" >&2
threadNumber=$OPTARG
;;
\?)
echo "Invalid option: -$OPTARG" >&2
exit 1
;;
:)
echo "Option -$OPTARG requires an argument." >&2
exit 1
;;
esac
done

if [ -z "$folder"  ]; then
	help
	exit 0
	fi



mkdir $folder;
echo "I put binaries in $folder";


cat src/bct_header.py>Bct.py
echo "Bct_INSTDIR = (\"$folder\")" >> Bct.py
cat src/bct_broken.py>>Bct.py
chmod +x Bct.py





git clone --recursive https://github.com/GATB/bcalm >>logCompile 2>>logCompile;
cd bcalm;
mkdir build; cd build;
cmake -DKSIZE_LIST="32 64 128" ..  >>logCompile 2>>logCompile;
make -j $threadNumber >>logCompile 2>>logCompile;
cp bcalm $folder;
cd ../..;
echo PHASE TWO, GRAPH CONSTRUCTION: BCALM;



git clone https://github.com/Malfoy/BGREAT2 >>logCompile 2>>logCompile;
cd BGREAT2;
make -j $threadNumber >>logCompile 2>>logCompile;
cp bgreat $folder;
cd ..;
echo PHASE THREE, READ MAPPING ON THE DBG: BGREAT;


git clone https://github.com/Malfoy/BTT >>logCompile 2>>logCompile;
cd BTT;
make -j $threadNumber >>logCompile 2>>logCompile;
if [ $? -ne 0 ]
       then
              echo "there was a problem with btrim compilation, check logs"
              exit 1
       fi
cp btt $folder;
cd ..;
echo PHASE FOUR GRAPH CLEANING: BTT;




echo The end !;
