# Source - https://stackoverflow.com/a
# Posted by dchakarov
# Retrieved 2025-12-10, License - CC BY-SA 3.0
timestamp=$(date +%s)

# Filenames
outname="out1"
outname_meta=$outname"_meta.tsv"

cargo run $outname.tirp $outname_meta

echo ""
echo "Compressing \"$outname.tirp\" into \"$outname.tirp.gz\"..."
bgzip -c $outname.tirp > $outname.tirp.gz 
echo ""

# Below is meant for beagle and running the sim

# Meant for the data-dir
echo "(for beagle) Sending results to /home/douglas/data/sim1"
cp $outname_meta /home/douglas/data/sim1
cp $outname.tirp /home/douglas/data/sim1
cp $outname.tirp.gz /home/douglas/data/sim1 

# Meant for the bascetRoot dir
echo "(for beagle) Sending results to /home/douglas/git/zorn"
cp $outname_meta /home/douglas/git/zorn
cp $outname.tirp /home/douglas/git/zorn
cp $outname.tirp.gz /home/douglas/git/zorn