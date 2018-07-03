#!/bin/bash
#commands to do a diva run

#This assumes the param.par, coast.cont,gvsampling.dat are all present in input dir already

usage()
{
cat << eof

USAGE:

   Script to perform the DIVA processing for the SOCAT reanalysis

   -d <dir>   | directory where the outputs should be copied to
   -f <fname> | filename that is to be interpolated
   -c <col>   | column number of input file to interpolate
   -p <file>  | default diva parameter file
   -h         | display this help 

eof
}

##########################################################
#These parameters are input from the python scripts via getopt
#where do we want our input/output files - this is dictated by our python processing
#Filename to interpolate
#filename=01_joined.nc.txt
#path to location of input data files
#pathname=~/process_scratch/oceanflux/diva/diva-4.6.10/DIVA3D/marktest
#Where we want our final output file to go
#destinationdir=~/process_scratch/oceanflux/diva/diva-4.6.10/DIVA3D/marktest/output/
##########################################################

while getopts "hd:f:c:p:t:" flag
do
   #echo "$flag" $OPTIND $OPTARG
   case "$flag" in
      d) destinationdir=$OPTARG ;;
      f) tmpname=$OPTARG;;
      c) column=$OPTARG;;
      p) defaultparameters=$OPTARG;;
      t) tempworkspaces=$OPTARG;;

      h) usage;exit 1;;
      ?) usage;exit 1;;
   esac
done

#Should do some testing on inputs here
if [ -z "$destinationdir" ]
then
   echo "Destination directory has not been given. Cannot interpolate file. Exiting."
   exit 1
fi

if [ -z "$tmpname" ]
then
   echo "File to interpolate has not been given. Cannot interpolate file. Exiting."
   exit 1
fi

if [ -z "$defaultparameters" ]
then
   echo "Default parameter file has not been given. Cannot interpolate file. Exiting."
   exit 1
fi

if [ ! -d "$destinationdir" ]
then
   echo "$destinationdir is not a directory. Cannot interpolate file. Exiting."
   exit 1
fi

if [ ! -e "$tmpname" ]
then
   echo "$tmpname does not exist. Cannot interpolate file. Exiting."
   exit 1
fi

if [ ! -e "$defaultparameters" ]
then
   echo "$defaultparameters does not exist. Cannot interpolate file. Exiting."
   exit 1
fi

if [ -z "$tempworkspaces" ] 
then
   echo "No temporary workspace directory given - will default to /tmp"
   tempworkspaces="/tmp"
fi

#convert to the absolute paths due to 'cd-ing' around
destinationdir=`readlink -f $destinationdir`
tmpname=`readlink -f $tmpname`
pathname=`dirname $tmpname`
filename=`basename $tmpname`


##########################################################
#These parameters are dependent on the DIVA install location
#where does the diva processing take place - this is dictated by the DIVA scripts
divainstallstripped=/home/cercache/project/oceanflux-shared/workspace/mwarren/diva-install/diva-4.6.11/DIVA3D/divastripped/
divainstallbin=/home/cercache/project/oceanflux-shared/workspace/mwarren/diva-install/diva-4.6.11/DIVA3D/bin/
##########################################################


##########################################################
#We shall copy the divainstall directory for each instance that we run
#so that it can be run on a grid / cloud based computer system
newtempinstall=`mktemp -d $tempworkspaces/diva-temporary-workspace-XXXXXXXX`
newtempinstallstripped=$newtempinstall/divastripped
#copy the install to this temporary directory
cp -r $divainstallbin $newtempinstall/
cp -r $divainstallstripped $newtempinstall/
#run path to use for executing (ie full path)
runpath=`readlink -f $newtempinstallstripped`
#input and output dirs
divainput=`readlink -f $newtempinstallstripped/input/`
divaoutput=`readlink -f $newtempinstallstripped/output/`
##########################################################

##########################################################
#These parameter names are derived from the script inputs
#full path to filename to interpolate
inputfile=$pathname/${filename}
#name to give the interpolated file
interpolatedfile=${filename}.anl
#name to give the interpolated error file
errorfile=${filename}-error.anl
##########################################################

#This next section will be looped through for each file in the input directory
#It does the actual interpolation and moves files around

#copy file to be interpolated into the diva workspace
cp ${inputfile} ${divainput}/data.dat.tmp

#we need to remove the top row and extract the columns we want (lat,lon,data)
cut -f 1,2,$column ${divainput}/data.dat.tmp | awk 'NR > 1 { print }' > ${divainput}/data.dat

#We should copy a parameter file with some defaults to the input directory
#else each consequtive run will use the previous run's settings (I think)
#in any case - we need one for the first run so may as well copy each run
cp ${defaultparameters} ${divainput}/param.par

#We also want a generalised cross validation values file
#make it here and replace the version that exists
echo "0.1" > ${divainput}/gcvsampling.dat
echo "1" >> ${divainput}/gcvsampling.dat
echo "10" >> ${divainput}/gcvsampling.dat
echo "100" >> ${divainput}/gcvsampling.dat
echo "500" >> ${divainput}/gcvsampling.dat

#Need to be in the install directory due to hard coded paths in the DIVA scripts
cd $newtempinstallstripped
#clean the workspace
$runpath/divaclean || exit 1
#fit the parameters
#$runpath/divafit -r || exit 1
#Need to check L is > 6 else it doesn't work?
#$runpath/divadatacoverage -n -r || exit
#create the mesh
$runpath/divamesh || exit 1
#estimate the signal to noise? (cross validation)
$runpath/divagcv -r || exit 1
#do the interpolation
$runpath/divacalc || exit 1

#copy the result to our requested destination
cp ${divaoutput}/fieldascii.anl ${destinationdir}/${interpolatedfile} || exit 1
cp ${divaoutput}/errorfieldascii.anl ${destinationdir}/${errorfile} || exit 1
#copy the param files also?
cp ${divainput}/param.par ${destinationdir}/${interpolatedfile}_parameters || exit 1

echo "Removing temporary workspace"
#Now clean up and remove the workspace
rm -rf $newtempinstall