#!/bin/bash
###Fast genome annotation for use with the JAMS pipeline
##Set environment
#Set Defaults and version
message_use="

fastannotate_JAMS.sh
By John McCulloch
Fast annotates fasta files using PROKKA.

Use: $(basename "$0") [options] -i contigs.fasta

-h Help
-v Version
-i File containing contigs or sequence to be annotated in FASTA format
-n Do not try and rename contigs if they are larger than 20 characters long. (Default is to rename them.)
-j Skip tbl2asn when using prokka. This step takes too long when there are over 10000 contigs.
-p Prefix
-t Number of threads to use. Default is to use all available by setting --cpu prokka option to 0.

"

version="1.3"

message_ver="

fastannotate_JAMS.sh ver $version
by John McCulloch NOV-2018

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
skiptbl="false"

#Get options
options=':i:t:p:jnvh'
while getopts $options option
do
    case $option in
        i  ) fasta="$OPTARG" ;;
        p  ) prefix="$OPTARG" ;;
        t  ) threads="$OPTARG" ;;
        n  ) nflag="true" ;;
        j  ) skiptbl="true" ;;
        v  ) version; exit;;
        h  ) usage; exit;;
        \? ) echo "Unknown option: -$OPTARG" >&2; exit 1;;
        :  ) echo "Missing option argument for -$OPTARG" >&2; exit 1;;
        *  ) echo "Unimplemented option: -$OPTARG" >&2; exit 1;;
    esac
done

shift $(($OPTIND - 1))

echo "Input file containing sequences to be annotated is $fasta"

#Find executables or die.
#Find gff2bed
deptest=$( { convert2bed -w | grep -c "convert2bed"; } 2>&1 )
if [ "$deptest" -lt "1" ]
then
    die "convert2bed was not found. Exiting now."
else
    echo "convert2bed is installed."
fi

#Find out size of headers.
if [ ! "$prefix" ]
then
    prefix=`echo $fasta | rev | cut -f 2- -d "." | rev`
fi

echo "Input file is $fasta"
echo "Prefix is $prefix"

namelen=`grep ">" "$fasta" | sort -z | tail -2 | head -1 | wc -c`
if [ "$namelen" -ge 20 ]
then
    long="true"
else
    long="false"
fi

#Rename contigs if they are long and -n flag is absent
if [ "$nflag" == "false" ] && [ "$long" == "true" ]
then
    echo "Will rename contig headers because they are too long."
    cat "$fasta" | sed s/"$prefix"/"d0e64cW3"/g > prktmp.fasta
elif [ "$nflag" == "true" ] && [ "$long" == "true" ]
then
    die "Contig headers are too long. Aborting now."
else
    echo "Contig headers do not need renaming."
fi


#Fix prokka if skipping tbl2asn
if [ "$skiptbl" == "true" ]
then
    #Long winded way for it to work in OSX and Linux. Sigh.
    prokkapointer=`which prokka`
    prokkadir=`readlink "$prokkapointer" | sed s/"prokka$"/""/g`
    linkdir=`echo "$prokkapointer" | sed s/"\/prokka$"/""/g`
    targetdir="$linkdir"/"$prokkadir"
    cat "$prokkapointer" | sed s/"tbl2asn -V b -a r10k"/"#tbl2asn -V b -a r10k"/g > "$targetdir"/prokkaJAMS
    chmod 775 "$targetdir"/prokkaJAMS
    prokkaexe="$targetdir"/prokkaJAMS
else
    prokkaexe="prokka"
fi


#Annotate
if [ "$long" == "true" ]
then
    "$prokkaexe" --cpus "$threads" --addgenes --outdir "$prefix"_PROKKA --prefix "$prefix" --locustag "d0e64cW3" --gcode 11 prktmp.fasta
else
    "$prokkaexe" --cpus "$threads" --addgenes --outdir "$prefix"_PROKKA --prefix "$prefix" --locustag "$prefix" --gcode 11 "$fasta"
fi

#Repair if headers were too long
if [ "$long" == "true" ]
then
    echo "Fixing prefixes in output files."
    cd "$prefix"_PROKKA
    cat "$prefix".gff | sed s/"d0e64cW3"/"$prefix"/g > tmp.gff
    mv tmp.gff "$prefix".gff
    cat "$prefix".gb[a-z] | sed s/"d0e64cW3"/"$prefix"/g > tmp.gbk
    mv tmp.gbk "$prefix".gbk
    cat "$prefix".faa | sed s/"d0e64cW3"/"$prefix"/g > tmp.faa
    mv tmp.faa "$prefix".faa
    cat "$prefix".ffn | sed s/"d0e64cW3"/"$prefix"/g > tmp.ffn
    mv tmp.ffn "$prefix".ffn
    cat "$prefix".fna | sed s/"d0e64cW3"/"$prefix"/g > tmp.fna
    mv tmp.fna "$prefix".fna
    cat "$prefix".fsa | sed s/"d0e64cW3"/"$prefix"/g > tmp.fsa
    mv tmp.fsa "$prefix".fsa
    cat "$prefix".tbl | sed s/"d0e64cW3"/"$prefix"/g > tmp.tbl
    mv tmp.tbl "$prefix".tbl
    cat "$prefix".tsv | sed s/"d0e64cW3"/"$prefix"/g > tmp.tsv
    mv tmp.tsv "$prefix".tsv
    cd ../
    rm prktmp.fasta
fi

#Generate BED file and make maps
cd "$prefix"_PROKKA/
convert2bed -i gff < "$prefix".gff | sort -k1,1 -k2,2n | tr -d "\'" > "$prefix".bed
#cat "$prefix".bed | tr -d "\'" | awk -F '[\t]' '{print $4"\t"($3-$2)"\t"$8"\t"$1}' > feature2contig.map
#cat "$prefix".bed | tr -d "\'" | awk -F'[\t;]' '/product/{for(i=1;i<=NF;++i)if($i~/product/)print $4"\t"$i}' | sed s/"product="/""/g | sed s/" "/"_"/g | sed $'s/[^- \tA-Za-z0-9]/_/g' | sed s/"__"/"_"/g | sed s/"__"/"_"/g | sed s/"_-"/"-"/g | sed s/"-_"/"-"/g | sed s/"_$"/""/g > feature2product.map
#cat "$prefix".bed | tr -d "\'" | awk -F'[\t;]' '/eC_number/{for(i=1;i<=NF;++i)if($i~/eC_number/)print $4"\t"$i}' | sed s/"eC_number="/"EC_"/g > CDS2EC.map

cd ../

#Final message
echo "Annotation finished."
