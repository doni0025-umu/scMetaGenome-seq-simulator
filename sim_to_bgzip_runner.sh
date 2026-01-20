# Source - https://stackoverflow.com/a
# Posted by dchakarov
# Retrieved 2025-12-10, License - CC BY-SA 3.0
timestamp=$(date +"%y%m%d_%H:%M")

# Filenames
outname=$timestamp"out"
outname_meta=$outname".1_meta.tsv"

# Args are name of output file, name of output metafile and lastly the number of bytes for chromosome-ref chunklength.
cargo run $outname.1.tirp $outname_meta 10000

echo ""
echo "Compressing \""$outname".1.tirp\" into \""$outname".1.tirp.gz\"..."
bgzip $outname.1.tirp 
echo ""

#######################################################################
#### Below is meant for beagle and running the sim ####################
#######################################################################

mkdir ~/data/sim/$outname
mkdir ~/data/sim/$outname/countsketch_mat.csv

# Meant for the data-dir
echo "(for beagle) Sending results to /home/douglas/data/sim"
mv $outname_meta /home/douglas/data/sim/$outname
mv $outname.1.tirp.gz /home/douglas/data/sim/$outname
cp run-setup.json "/home/douglas/data/sim/"$outname"/"$outname"run-setup.json"
echo ""

## Meant for the bascetRoot dir
#echo "(for beagle) Sending results to /home/douglas/git/zorn"
#cp $outname_meta /home/douglas/git/zorn
#cp $outname.tirp /home/douglas/git/zorn
#cp $outname.tirp.gz /home/douglas/git/zorn

Rscript ~/data/sim/proc_data.R $outname