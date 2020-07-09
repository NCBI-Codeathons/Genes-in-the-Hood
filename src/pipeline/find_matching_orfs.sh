#!/bin/bash
set -o pipefail
set -o errexit

source_dir="/usr/local/data/assemblies/packages"

extraction_dir="/usr/local/data/assemblies/extracted"

results_dir="Results"

verbose=true

usage()
{
cat <<EOF
Usage:
    $script_name [options] [TARGET_DIR]

Valid options:
    -h | --help                      Display this message

    -a | --accession ARG             Assembly accession

    -m | --model                     HMM model file to use

    -r | --results                   Top-level dir for results (default: $results_dir)

    -q | --quiet                     No output to stdout

EOF
    exit $1
}

accession=""
model=""

while [[ $# -gt 0 ]]
do
    case $1 in
    -h | --help)
        usage 0
        ;;
    -a | --accession)
        accession="$2"
        shift
        ;;
    -m | --model)
        model="$2"
        shift
        ;;
    -r | --results)
        results_dir="$2"
        shift
        ;;
    -q | --quiet)
        verbose=false
        ;;
    -*)
        echo "invalid option : '$1'"
        usage 10
        ;;
    *)
        echo "$script_name does not accept any positional arguments"
        usage 10
    esac
    shift
done

if [ -z "$accession" ]; then
    echo "Error missing accession value"
    usage 1
elif [ -z "$model" ]; then
    echo "Error missing model file"
    usage 1
elif [ ! -f "$model" ]; then
    echo "Error missing model file"
    usage 1
fi

mkdir -p $extraction_dir

if [ ! -d "$extraction_dir/$accession" ]; then
    unzip $source_dir/$accession.zip -d $extraction_dir/$accession
fi

if [ ! -d "$extraction_dir/$accession" ]; then
    echo "Error unable to locate $accession data"
    usage 1
fi
output_dir=$results_dir/$accession
input_files=`find $extraction_dir/$accession -name '*.fna' `
gff_files=`find $extraction_dir/$accession -name '*.gff' `

if $verbose; then
    echo "Process $accession using $model"
    echo "output: $output_dir"
    echo "input_files: $input_files"
fi

mkdir -p $output_dir
ORFfinder -in <(cat $input_files) -g 11 -s 0 -ml 300 -n t -outfmt 0 -out $output_dir/orf.faa > $output_dir/ORFfinder.stdout 2> $output_dir/ORFfinder.stderr

hmmsearch --tblout $output_dir/hmmsearch.out $model $output_dir/orf.faa > $output_dir/hmmsearch.stdout 2> $output_dir/hmm_search.stderr 

awk 'BEGIN {OFS="\t"} /^[^#]/ {split($1, a, "ORF[0-9]*_"); location=a[1]; split(a[2], b, ":"); print b[1], b[2], b[3], "name="$1"; query="$4"; e-value="$5"; score="$6"; bias="$7}' $output_dir/hmmsearch.out | 
    tee $output_dir/locations.meta |
    cut -f 1,2,3 > $output_dir/locations.bed

cat $output_dir/locations.bed

#if [ ! -z "$gff_files" ]; then 
#    TODO
#fi

