#!/bin/bash
##Wrapper for building a custom Kraken2 database wil all bacteria, archaea, fungi, viruses, protozoa and the human and mouse genomes, for use with JAMS.
########################
## Define parameters  ##
########################
#Set particular parameters
start=$(date "+%Y-%m-%d %H:%M:%S")
startclock=$(date +%s)

##Get options
message_use="
buildJAMSk2db.sh
By John McCulloch
Wrapper for building a custom Kraken2 database wil all bacteria, archaea, fungi, viruses, protozoa and the human and mouse genomes, for use with JAMS.

Use: $(basename "$0") -p <path/to/outdir> -n <dbname> -t <num_threads>

-h Help;
-v Version;
-p Complete path to out directory (will be created if does not exist);
-n Database name [Default=JAMSk2db]
-t Number of threads to use [Default=32]
-m Maximum size of database, in Gb. [Default=32]
-a Test dependencies only and exit.

"
vernum="0.91"
verdate="MAR-2019"
message_ver="buildJAMSk2db.sh ver $vernum by John McCulloch $verdate"

#Define defaults
testandleave="false"
dbname="JAMSk2db"
maxsizeGB=32

#Get options
usage () { echo "$message_use"; }
version () { echo "$message_ver"; }
function die {
    echo "$@" >&2
    exit 1
}

options=':p:n:m:t:avh'
while getopts $options option
do
    case $option in
        p  ) dbpath="$OPTARG" ;;
        n  ) dbname="$OPTARG" ;;
        t  ) threads="$OPTARG" ;;
        m  ) maxsizeGB="$OPTARG" ;;
        a  ) testandleave="true" ;;
        v  ) version; exit ;;
        h  ) usage; exit ;;
        \? ) echo "Unknown option: -$OPTARG" >&2; exit 1;;
        :  ) echo "Missing option argument for -$OPTARG" >&2; exit 1;;
        *  ) echo "Unimplemented option: -$OPTARG" >&2; exit 1;;
    esac
done

shift $(($OPTIND - 1))

#Define functions
function report2log {
    justnow=$(date "+%Y-%m-%d %H:%M:%S") 
    echo "${justnow} -> $@" >> "$dbpath"/"$dbname"_Build_log.txt
    echo "$@"
}

function taxid2gensize {
        gs=`cat "$1" | grep -v "^>" | tr -d '\n' | wc -m`
        taxid=`head -1 "$1" | rev | cut -f 1 -d "|" | rev`
        echo -e $taxid'\t'$gs
}

renameheader(){
    taxid=`cat file2taxid.map | grep -w $1 | cut -f 2`
    mv "$1" "$1".tmp
    cat "$1".tmp | awk -F'[ ]' '/^>/{acc=$1}{print $1}' | sed -e '/^>/ s/$/\|kraken:taxid\|brObdiNGnag/' | sed s/"brObdiNGnag"/$taxid/g > "$1"
    rm "$1".tmp
}

#Set paths
thisdir=`pwd`
db="$dbpath"/"$dbname"

#Quit at this point if only testing dependencies.
if [ "$testandleave" == "true" ]
then
    die "Only testing dependencies for $message_ver ; Aborting now."
fi

#Start the clock
start=$(date +%s)
#Start the news
report2log "Starting to build $dbname db using buildJAMSk2db.sh ver $vernum"
report2log "Using $threads threads"
report2log "Database full path is $db"
maxsizebytes=`echo "$maxsizeGB" | awk '{print ($1 * 1000000000)}'`
report2log "Database maximum size is $maxsizeGB Gb, or $maxsizebytes bytes."

#Create directory if it does not exist
mkdir -p "$db"
cd "$db"

###Download taxononomy
report2log "Downloading Taxonomy."
kraken2-build --download-taxonomy --db "$db"
taxdir="$db"/taxonomy

#2. Get the Genomes
#Set directory structures for the genomes in fasta format
mkdir "$db"/genomes
echo -e Taxid'\t'Size > taxid2size.map

#2.1 Get genomes for bugs
#Download bacteria, archaea, fungi, viral and protozoa genomes
for bug in bacteria archaea fungi viral protozoa vertebrate_mammalian
do
    mkdir "$db"/genomes/"$bug"
    cd "$db"/genomes/"$bug"
    report2log "Getting $bug genomes"
    #get list of strains available
    report2log "Getting availabale links in GenBank..."
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/"$bug"/assembly_summary.txt

    if [ "$bug" == "bacteria" ]
    then
        #For bacteria get only complete genomes and scaffolds
        #Generate list of complete genomes
        cat assembly_summary.txt | awk -F "\t" '$11=="latest" && $12=="Complete Genome"{print $20}' | awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' > ftpfilepaths.list
        #Add Scaffold paths to list
        cat assembly_summary.txt | awk -F "\t" '$11=="latest" && $12=="Scaffold"{print $20}' | awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' >> ftpfilepaths.list

    elif [ "$bug" == "vertebrate_mammalian" ]
    then
        cat assembly_summary.txt | awk -F "\t" '$8=="Homo sapiens" && $12=="Chromosome"{print $20}' | awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' > ftpfilepaths.list
        cat assembly_summary.txt | awk -F "\t" '$8=="Mus musculus" && $12=="Chromosome"{print $20}' | awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' >> ftpfilepaths.list
    else
        #For other taxa, get all genomes available
        cat assembly_summary.txt | awk -F "\t" '$11=="latest" {print $20}' | awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' > ftpfilepaths.list
    fi

    howmany=`cat ftpfilepaths.list | wc -l | rev | cut -f 1 -d " " | rev`
    report2log "There are $howmany $bug genomes to download."

    #generate taxidmap
    cat assembly_summary.txt | awk -F "\t" '$11=="latest"{print $20"_genomic.fna\t"$6}' | sed -e s/"ftp:\/\/ftp.ncbi.nlm.nih.gov\/genomes\/all\/"/""/g > file2taxid.map

    cat assembly_summary.txt | awk -F "\t" '$11=="latest"{print $20"_genomic.fna\t"$7}' | sed -e s/"ftp:\/\/ftp.ncbi.nlm.nih.gov\/genomes\/all\/"/""/g > file2speciestaxid.map

    cat assembly_summary.txt | awk -F "\t" '$11=="latest"{print $20"_genomic.fna\t"$1}' | sed -e s/"ftp:\/\/ftp.ncbi.nlm.nih.gov\/genomes\/all\/"/""/g > file2accession.map

    list=`cat ftpfilepaths.list`

    #Download
    report2log "Downloading $bug genomes."
    startdown=$(date +%s)
    for genome in $list
    do
        wget $genome
    done

    report2log "$bug genomes have been downloaded."
    now=$(date +%s)
    time=`echo $((now-startdown)) | awk '{print ($1/3660)}'`
    report2log "This took $time hours"

    #Decompress
    report2log "Decompressing $bug genomes."
    startdec=$(date +%s)

    N=$threads
    (
    for file in *.gz; do 
       ((i=i%N)); ((i++==0)) && wait
       gunzip "$file" &
    done
    )

    report2log "$bug genomes have been decompressed."
    now=$(date +%s)
    time=`echo $((now-startdec)) | awk '{print ($1/3660)}'`
    report2log "This took $time hours"

    #rename sequence and add taxid for kraken db
    report2log "Adding NCBI TaxIDs to headers."
    startren=$(date +%s)

    #cycle through files and adjust headers
    if [ "$bug" != "vertebrate_mammalian" ]
    then
        N=$threads
        (
        for file in *.fna; do 
           ((i=i%N)); ((i++==0)) && wait
           renameheader "$file" &
        done
        )
    else
        #For the only two vertebrate genomes, do not rename in parallel, as so not to prematurely start the next phase.
        for file in *.fna
        do
            renameheader "$file"
        done
    fi
    
    report2log "Headers for all $bug genomes have been renamed."
    now=$(date +%s)
    time=`echo $((now-startren)) | awk '{print ($1/3660)}'`
    report2log "This took $time hours"

    #Count genome size
    report2log "Counting number of basepairs for each genome"
    startcnt=$(date +%s)

    #cycle through files and count the bases
    if [ "$bug" != "vertebrate_mammalian" ]
    then
        N=$threads
        (
        for file in *.fna; do 
           ((i=i%N)); ((i++==0)) && wait
           taxid2gensize "$file" >> "$db"/taxid2size.map &
        done
        )
    fi

    report2log "The number of basepairs for all $bug genomes have been counted."
    now=$(date +%s)
    time=`echo $((now-startcnt)) | awk '{print ($1/3660)}'`
    report2log "This took $time hours"

done

##2. Add genomes to dbase
cd "$db"
startadd=$(date +%s)
for group in bacteria archaea fungi viral protozoa vertebrate_mammalian
do
    report2log "Adding $group genomes to db."
    find genomes/"$group"/ -name '*.fna' -print0 | xargs -0 -P "$threads" -I{} -n1 kraken2-build --add-to-library {} --db $db
done

report2log "All genomes have been added to the db."
now=$(date +%s)
time=`echo $((now-startadd)) | awk '{print ($1/3660)}'`
report2log "This took $time hours"

##3. Build the database
report2log "Will now build the entire database. Please be VERY patient."
kraken2-build --build --threads "$threads" --db "$db" --max-db-size "$maxsizebytes"
report2log "The $dbname database has been built."
now=$(date +%s)
time=`echo $((now-startclock)) | awk '{print ($1/3660)}'`
report2log "This took $time hours overall."
dbtimestamp=$(date +%C%y%m%d%H%M)
#echo JAMSkdb_"$vernum"_"$dbtimestamp" > "$db"/JAMSKdb.ver
