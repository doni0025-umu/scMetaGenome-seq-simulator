#!/bin/bash

# Source - https://stackoverflow.com/a
# Posted by dchakarov
# Retrieved 2025-12-10, License - CC BY-SA 3.0

# Shorthand for grabbing a fasta from NCBI datasets

descript=""
seqdepth=1
overwrite_old=false
simdir="/husky/douglas/sim"
forreport=false

startinfo=""

usage() {
    echo ""
    echo "  -h,   displays this help message."
    echo "  -d,   [sequencing_depth], specifies sequencing depth for this run. For example, entering 10 will give a simulation where each base is sequenced, roughly, 10 times. [default=1]."
    echo "  -n,   [DESCRIPT], adds DESCRIPT in front of the rundir."
    echo "  -o,   overwrites previous oldjob folder with produced results [default=false]".
    echo "  -r,   job will end up in the FOR_REPORT dir [default=false]".
    echo ""
}

while getopts "hd:n:or" flag; do
 case $flag in
   h)   # Handle the -h flag
        # Display script help information
        usage
        exit 0
    ;;
   d)   # Handle the -d flag with an argument
        startinfo=$startinfo" - Run will start with sequencing depth $OPTARG for every base.\n"
        seqdepth=$OPTARG
   ;;
   o)   # Handles the -o flag
        startinfo=$startinfo" - Oldjob will be overwritten with produced results.\n"
        overwrite_old=true
   ;;
   n)   # Handle the -n flag with an argument
        startinfo=$startinfo" - Description $OPTARG will be added.\n"
        descript="$OPTARG"
   ;;
   r)   # Handle the -r flag
        startinfo=$startinfo" - Job will reach FOR_REPORT dir.\n"
        forreport=true
   ;;
 esac
done


echo "The simulation run will start with the following settings:"
echo ""
printf "$startinfo"
echo ""
sleep 7s


timestamp=$(date +"%y%m%dT%H_%M")
composition_name=$(basename ncbi_dataset/*_data _data)

# Filenames
tirpname=$descript"_"$composition_name"_"$timestamp
runname="seqD"$seqdepth"X_"$tirpname
rundir=$runname
if $forreport
     then
          rundir="FOR_REPORT/jobs/"$runname
fi

#rundir=oldjob
rundir_meta=meta_$tirpname".1.tsv"

# Args are name of output file and name of output metafile
cargo run /husky/douglas/$tirpname.1.tirp $rundir_meta $seqdepth

echo ""
echo "Compressing \""$tirpname".1.tirp\" into \""$tirpname".1.tirp.gz\"..."
bgzip -@ 40 /husky/douglas/$tirpname.1.tirp 
echo ""

#######################################################################
#### Below is meant for beagle and running the sim ####################
#######################################################################

mkdir -p $simdir/$rundir
mkdir -p $simdir/$rundir/meta_and_setup

# Meant for the data-dir
echo "(for beagle) Sending results to $simdir"

mv $rundir_meta $simdir/$rundir/meta_and_setup
mv /husky/douglas/$tirpname.1.tirp.gz $simdir/$rundir
cp active_run_params/*.json $simdir/$rundir/meta_and_setup

echo ""

# Make countsketch and produce UMAP and feature plots
Rscript $simdir/proc_data.R $rundir $tirpname $runname

cp $simdir/$rundir/*umapGOLD* $simdir/UMAPs
cp $simdir/$rundir/*umapSEURAT* $simdir/UMAPs

if $overwrite_old
    then
        # Clearing previous oldjob
        rm -r $simdir/oldjob/*
        #mkdir -p $simdir/oldjob/countsketch_mat.csv # Deprecated use since 260205
        mkdir -p $simdir/oldjob/meta_and_setup
        mkdir -p $simdir/oldjob/svg_plots
        
        # Depositing results into oldjob
        cp $simdir/$rundir/$tirpname.1.tirp.gz $simdir/oldjob/oldjob.1.tirp.gz
        cp active_run_params/*.json $simdir/oldjob/meta_and_setup/oldjob_run-composition.json
        cp active_run_params/*.txt $simdir/oldjob/meta_and_setup/oldjob_accession.txt
        cp $simdir/$rundir/meta_and_setup/$rundir_meta $simdir/oldjob/meta_and_setup/meta_oldjob.1.tsv
        cp $simdir/$rundir/countsketch_mat.csv $simdir/oldjob/
        cp $simdir/$rundir/meta_and_setup/session* $simdir/oldjob/meta_and_setup/sessionInfo_oldjob.txt

        # Plots
        cp $simdir/$rundir/*umapSEURAT* $simdir/oldjob/umapSEURAT_oldjob.svg
        cp $simdir/$rundir/*umapGOLDSTANDARD* $simdir/oldjob/umapGOLDSTANDARD_oldjob.svg
        cp $simdir/$rundir/*feature* $simdir/oldjob/featureplot_oldjob.svg


fi

if $forreport
     then
          find $simdir/$rundir/*umapSEURAT* -exec cp {} $simdir/FOR_REPORT/UMAP_SEURAT \;
          find $simdir/$rundir/*umapGOLDSTANDARD* -exec cp {} $simdir/FOR_REPORT/UMAP_GOLDSTANDARD \;
          find $simdir/$rundir/*feature* -exec cp {} $simdir/FOR_REPORT/FEATUREPLOT \;
          find $simdir/$rundir/*confusion* -exec cp {} $simdir/FOR_REPORT/CONFUSIONMATRIX \;
          find $simdir/$rundir/*collected* -exec cp {} $simdir/FOR_REPORT/COLLECTEDPLTS \;
fi

#cargo run --manifest-path /home/douglas/spectra_data_countsketch_calc/Cargo.toml $rundir
#Rscript plot_costhetas.R $rundir


# An effort to save space
rm $simdir/$rundir/$tirpname.1.tirp.gz
