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
    jq . ncbi_dataset/data/assembly_data_report.jsonl > pretty_asmbly_jsons/asmbly$accession.json
    jq . ncbi_dataset/data/$accession/sequence_report.jsonl > pretty_seq_jsons/seq$accession.json
    mv ncbi_dataset/data/assembly_data_report.jsonl ncbi_dataset/data/$accession/
    mv ncbi_dataset/data/dataset_catalog.json ncbi_dataset/data/$accession/


done


# 
