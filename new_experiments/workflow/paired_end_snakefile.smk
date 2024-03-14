rule all:
    input:
        "data/bloodCulture_qiime/bloodCulture_paired_end_seqs.qza"


rule import_paired_end:
    input:
        seqs_directory = directory("data/bloodCulture_qiime/raw_seqs")
    output:
        seqs = "data/bloodCulture_qiime/bloodCulture_paired_end_seqs.qza"
    conda:
        "qiime2-2023.5"
    shell:
        """
        qiime tools import \
            --type EMPPairedEndSequences \
            --input-path {input.seqs_directory} \
            --output-path {output.seqs}
        """