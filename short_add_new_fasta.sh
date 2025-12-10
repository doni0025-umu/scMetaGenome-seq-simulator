# Shorthand for grabbing a fasta from NCBI datasets

for accession in "$@" 
do
    echo "$accession"
    # Below block is from ncbi datasets at nih
    datasets download genome accession $accession --include genome,seq-report
    #endblock

    # Structure it in the program hierarchy
    unzip -o ncbi_dataset.zip
    rm -rf ncbi_dataset.zip
    mv ncbi_dataset/data/assembly_data_report.jsonl ncbi_dataset/data/$accession/
    mv ncbi_dataset/data/dataset_catalog.json ncbi_dataset/data/$accession/

done

# Current list:
# GCF_000005845.2 GCF_000006945.2 GCF_000013425.1 GCF_000161615.1 GCF_000412675.1
