# Source - https://stackoverflow.com/a
# Posted by dchakarov
# Retrieved 2025-12-10, License - CC BY-SA 3.0
timestamp=$(date +%s)

cargo run out.tirp out_meta.tsv
bgzip -c out.tirp > out.tirp.gz