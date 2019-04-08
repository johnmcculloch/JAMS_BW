#!/bin/bash
###Counts fastqs in parallel
#Set Defaults and version
message_use="

countfastqJAMS.sh
By John McCulloch
Counts fastqs in parallel

Use: $(basename "$0") [options] -d /path/to/folder/containing/only/fastqs

-h Help
-v Version
-d /path/to/folder/containing/only/fastqs
-t Number of threads to use.

"

version="1.0"

message_ver="

fastannotate_JAMS.sh ver $version
by John McCulloch OCT-2018

"

#Define functions
usage () { echo "$message_use"; }
version () { echo "$message_ver"; }
function die {
    echo "$@" >&2
exit 1
}

#Set defaults
threads=0
nflag="false"

#Get options
options=':d:t:vh'
while getopts $options option
do
    case $option in
        d  ) fastqfolder="$OPTARG" ;;
        t  ) threads="$OPTARG" ;;        
        v  ) version; exit;;
        h  ) usage; exit;;
        \? ) echo "Unknown option: -$OPTARG" >&2; exit 1;;
        :  ) echo "Missing option argument for -$OPTARG" >&2; exit 1;;
        *  ) echo "Unimplemented option: -$OPTARG" >&2; exit 1;;
    esac
done

shift $(($OPTIND - 1))

countfastq(){
    stats=`cat $1 | awk '((NR-2)%4==0){read=$1;total++;count[read]++;sum+=length($1)}END{print total"§"sum}'`
    fn="$1"
    echo -e $fn'§'$stats | sed $'s/§/\t/g' 
}
export -f countfastq

countfasta(){
    cat $1 | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' 
}

cd "$fastqfolder"
fqlist=`ls *.fastq`
parallel -j "$threads" -k countfastq ::: "$fqlist" >> reads.stats
