# Source - https://stackoverflow.com/a
# Posted by dchakarov
# Retrieved 2025-12-10, License - CC BY-SA 3.0

for i in {1..1}
do
    timestamp=$(date +"%y%m%d_%H:%M")

    # Filenames
    outname=gut_1000_cells$timestamp
    #outname=oldjob
    outname_meta=meta_$outname".1.tsv"

    # Args are name of output file, name of output metafile and lastly the number of bytes for chromosome-ref chunklength.
    cargo run $outname.1.tirp $outname_meta 10000

    echo ""
    echo "Compressing \""$outname".1.tirp\" into \""$outname".1.tirp.gz\"..."
    bgzip $outname.1.tirp 
    echo ""

    #######################################################################
    #### Below is meant for beagle and running the sim ####################
    #######################################################################

    mkdir -p ~/data/sim/$outname
    mkdir -p ~/data/sim/$outname/countsketch_mat.csv

    # Meant for the data-dir
    echo "(for beagle) Sending results to /home/douglas/data/sim"
    mv $outname_meta /home/douglas/data/sim/$outname
    mv $outname.1.tirp.gz /home/douglas/data/sim/$outname
    cp run-setup.json "/home/douglas/data/sim/"$outname"/"$outname"run-setup.json"
    
    echo ""

    # Make countsketch and produce UMAP and feature plots
    Rscript ~/data/sim/proc_data.R $outname

    cp /home/douglas/data/sim/$outname/umap* /home/douglas/data/sim/chromvsplasmid_UMAPs
    cp /home/douglas/data/sim/$outname/*.json /home/douglas/data/sim/chromvsplasmid_UMAPs

    # An effort to save space
    rm /home/douglas/data/sim/$outname/$outname.1.tirp.gz

done