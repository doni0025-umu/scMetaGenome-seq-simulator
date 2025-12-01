cargo run testThur.fna hackThurREMOVEeod.tirp
bgzip -c testThur.tirp > testThur.tirp.gz

cargo run ncbi_dataset/data/GCF_000006945.2/GCF_000006945.2_ASM694v2_genomic.fna firstout.tirp
bgzip -c firstout.tirp > firstout.tirp.gz