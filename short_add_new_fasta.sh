# Shorthand for grabbing a fasta from NCBI datasets

rm -r ncbi_dataset

accessions= cat -v accessions.txt

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

done < accessions.txt
