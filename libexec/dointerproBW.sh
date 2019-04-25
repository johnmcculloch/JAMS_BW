#! /bin/bash
#Wrapper for running interpro on Biowulf
##Define defaults and parameters
message_use="

dointerproBW.sh
By John McCulloch
Runs Interproscan on Biowulf (module by W. Resch) using as input a fasta file.

Use: $(basename "$0") [options] -i aminoacids.fasta

-h Help
-v Version
-i File aminoacid sequences in FASTA format
-c Chunk size for splitting queries (200-1000). Default=500.
-p Prefix. (Default is to use file basename preceding last . as prefix)
-s Simulate, do not actually launch interpro sbatch job on Biowulf.

"

version="0.33"

message_ver="

dointerproBW.sh ver $version
by John McCulloch APR-2019

"
#Set defaults
chunk="500"
outdir=""
simflag="false"

#Get options
usage () { echo "$message_use"; }
version () { echo "$message_ver"; }
function die {
    echo "$@" >&2
    exit 1
}

options=':i:p:c:svh'
while getopts $options option
do
    case $option in
        i  ) fasta="$OPTARG" ;;
        p  ) prefix="$OPTARG" ;;
        c  ) chunk="$OPTARG" ;;
        s  ) simflag="true" ;;
        v  ) version; exit;;
        h  ) usage; exit;;
        \? ) echo "Unknown option: -$OPTARG" >&2; exit 1;;
        :  ) echo "Missing option argument for -$OPTARG" >&2; exit 1;;
        *  ) echo "Unimplemented option: -$OPTARG" >&2; exit 1;;
    esac
done

shift $(($OPTIND - 1))

#Find executables or die.
echo "You are on Biowulf."
module load java/1.8.0_211 || die "Could not load java module."
module load interproscan/5.29-68.0 || die "Could not load interproscan module."

#Generate sbatch script
interproscan --goterms --pathways -f tsv "$fasta" interproscan "$chunk" 1> iprocommand.txt

if [ "$simflag" == "false" ]
then
    cd interproscan
    sbatch interproscan.batch 1> ipro.job
    sbjob=`cat ipro.job | head -1`
    echo "Interpro sbatch job is $sjob"
    cd ../
else
    echo "This is only a simulation. You can run the sbatch script later."
fi
