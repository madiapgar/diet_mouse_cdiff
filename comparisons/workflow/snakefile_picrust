rule all:
    input:
        "data/misc/dna-sequences.fasta",
        "data/picrust/tmp_out_pipeline",
        "data/picrust/out_pipeline"


rule rep_seqs2fasta:
    input:
        "data/misc/all_comp_rep_seqs.qza"
    output:
        "data/misc/dna-sequences.fasta"
    conda:
        "qiime2-2023.5"
    shell:
        """
        qiime tools export \
            --input-path {input} \
            --output-path ./data/misc
        """

rule picrust2:
    input:
        seqs_fasta = "data/misc/dna-sequences.fasta",
        biom = "data/misc/comp_tss_table.biom"
    output:
        tmp_out = directory("data/picrust/tmp_out_pipeline"),
        out = directory("data/picrust/out_pipeline")
    conda:
        "picrust2"
    shell:
        """
        picrust2_pipeline.py \
            -s {input.seqs_fasta} \
            -i {input.biom} \
            -o {output.tmp_out} \
            --stratified \
            --per_sequence_contrib \
            -p 32
        
        mv {output.tmp_out} {output.out}
        """