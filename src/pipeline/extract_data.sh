#!/bin/bash
set -o pipefail
set -o errexit

source_dir="/usr/local/data/assemblies/packages"

extraction_dir="/usr/local/data/assemblies/extracted"

results_dir="Results"

header=false

usage()
{
cat <<EOF
Usage:
    $script_name [options] [TARGET_DIR]

Valid options:
    -h | --help                      Display this message

    -a | --accession ARG             Assembly accession

    -r | --results                   Top-level dir for results (default: $results_dir)

    --header                         output header line

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
    -r | --results)
        results_dir="$2"
        shift
        ;;
   --header)
        header=true
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
fi

output_dir=$results_dir/$accession
gff_files=`find $extraction_dir/$accession -name '*.gff' `

if $header; then
    echo -e "#assembly\taccession\tstart\tstop\tname\tquery\te-value\tscore\tbias"
fi

if [ ! -z "$gff_files" ]; then
    have_genes=`awk 'BEGIN {OFS="\t"} /^[^#]/ {split($1, a, "ORF[0-9]*_"); location=a[1]; split(a[2], b, ":"); if ($5 <= 10E-20) {print assmacc, b[1], b[2], b[3], $1, $4, $5, $6, $7}}' assmacc=$accession $output_dir/hmmsearch.out | 
        python ../find_genes_by_location/find_genes_by_location.py  --packages-dir $source_dir  | tr ',' '\t'`
fi

if [ -z "$have_genes" ]; then
    awk 'BEGIN {OFS="\t"} /^[^#]/ {split($1, a, "ORF[0-9]*_"); location=a[1]; split(a[2], b, ":"); if ($5 <= 10E-20) {print assmacc, b[1], b[2], b[3], $1, $4, $5, $6, $7, "\t\t\t\t"}}' assmacc=$accession $output_dir/hmmsearch.out 
else
    echo "$have_genes"
fi

