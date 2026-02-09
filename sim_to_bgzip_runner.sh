#!/bin/bash

# Source - https://stackoverflow.com/a
# Posted by dchakarov
# Retrieved 2025-12-10, License - CC BY-SA 3.0

# Shorthand for grabbing a fasta from NCBI datasets

prefix=""
overwrite_old=false
usage() {
    echo ""
    echo "  -h,   displays this help message."
    echo "  -n,   [PREFIX], adds PREFIX in front of the outname."
    echo "  -o,   overwrites previous oldjob folder with produced results [default=false]".
    echo ""
}

while getopts "hn:o" flag; do
 case $flag in
   h)   # Handle the -h flag
        # Display script help information
        usage
        exit 0
    ;;
   o)   # Handles the -o flag
        echo "Oldjob will be overwritten"
        overwrite_old=true
   ;;
   n)   # Handle the -n flag with an argument
        echo "Prefix $OPTARG will be added."
        prefix="$OPTARG""_"
   ;;
 esac
done

for i in {1..1}
do
    timestamp=$(date +"%y%m%d__%H_%M")
    descript=$(basename active_run_params/*run-composition.json -run-composition.json)

    # Filenames
    outname=$prefix$descript"_"$timestamp
    #outname=oldjob
    outname_meta=meta_$outname".1.tsv"

    # Args are name of output file and name of output metafile
    cargo run /husky/douglas/$outname.1.tirp $outname_meta

    echo ""
    echo "Compressing \""$outname".1.tirp\" into \""$outname".1.tirp.gz\"..."
    bgzip -@ 40 /husky/douglas/$outname.1.tirp 
    echo ""

    #######################################################################
    #### Below is meant for beagle and running the sim ####################
    #######################################################################

    mkdir -p /husky/douglas/sim/$outname
    #mkdir -p /husky/douglas/sim/$outname/countsketch_mat.csv
    mkdir -p /husky/douglas/sim/$outname/meta_and_setup

    # Meant for the data-dir
    echo "(for beagle) Sending results to /husky/douglas/sim"
    
    mv $outname_meta /husky/douglas/sim/$outname/meta_and_setup
    mv /husky/douglas/$outname.1.tirp.gz /husky/douglas/sim/$outname
    cp active_run_params/*run-composition.json /husky/douglas/sim/$outname/meta_and_setup
    
    echo ""

    # Make countsketch and produce UMAP and feature plots
    Rscript /husky/douglas/sim/proc_data.R $outname

    cp /husky/douglas/sim/$outname/umap* /husky/douglas/sim/chromvsplasmid_UMAPs
    cp active_run_params/*.json /husky/douglas/sim/chromvsplasmid_UMAPs

    if $overwrite_old
        then
            # Clearing previous oldjob
            rm -r /husky/douglas/sim/oldjob/*
            #mkdir -p /husky/douglas/sim/oldjob/countsketch_mat.csv # Deprecated use since 260205
            mkdir -p /husky/douglas/sim/oldjob/meta_and_setup
            
            # Depositing results into oldjob
            cp /husky/douglas/sim/$outname/$outname.1.tirp.gz /husky/douglas/sim/oldjob/oldjob.1.tirp.gz
            cp active_run_params/*run-composition.json /husky/douglas/sim/oldjob/meta_and_setup/oldjob_run-composition.json
            cp active_run_params/*.txt /husky/douglas/sim/oldjob/meta_and_setup/oldjob_accession.txt
            cp /husky/douglas/sim/$outname/meta_and_setup/$outname_meta /husky/douglas/sim/oldjob/meta_and_setup/meta_oldjob.1.tsv
            cp /husky/douglas/sim/$outname/meta_and_setup/session* /husky/douglas/sim/oldjob/meta_and_setup/sessionInfo_oldjob.txt

            # Plots
            cp /husky/douglas/sim/$outname/umap* /husky/douglas/sim/oldjob/umap_oldjob.jpg
            cp /husky/douglas/sim/$outname/feature* /husky/douglas/sim/oldjob/featureplot_oldjob.jpg


    fi

    cargo run --manifest-path /home/douglas/spectra_data_countsketch_calc/Cargo.toml $outname
    Rscript plot_costhetas.R $outname


    # An effort to save space
    rm /husky/douglas/sim/$outname/$outname.1.tirp.gz

done