#!/bin/bash

# Shorthand for grabbing a fasta from NCBI datasets


echo "Are you sure? Your active run params will be deleted in favor of your new params!"
echo " - press [y] to run"
echo " - press any other key to abort"
read -n 1
echo ""


usage() {
    echo ""
    echo "  -h,   displays this help message."
    echo "  -c,   [NUM_OF_CELLS], the (rough) number of SPCs to simulate."
    echo ""
}

while getopts "hc:" flag; do
 case $flag in
   h)   # Handle the -h flag
        # Display script help information
        usage
        exit 0
    ;;
   c)   # Handle the -c flag with an argument
        echo "$OPTARG cells will be set up for simulation."
        num_cells=$OPTARG
   ;;
 esac
done

   # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
    rm -rf ncbi_dataset
    rm -f active_run_params/*-run-composition.json
    Rscript setup_accessions_and_composition.R $num_cells


    accessions= cat -v active_run_params/accessions.txt

    while read accession
    do
        echo "$accession"
        # Below block is from ncbi datasets at nih
        datasets download genome accession $accession --include genome,seq-report
        #endblock

        # Structure it in the program hierarchy
        unzip -o ncbi_dataset.zip -x README.md md5sum.txt
        rm -rf ncbi_dataset.zip
        mv ncbi_dataset/data/assembly_data_report.jsonl ncbi_dataset/data/$accession/
        mv ncbi_dataset/data/dataset_catalog.json ncbi_dataset/data/$accession/

    done < active_run_params/accessions.txt

fi