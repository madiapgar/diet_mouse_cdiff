rule all:
    input:
        "data/picrust/out_pipeline",
        "data/picrust/meta_contrib.tsv",
        "plots/butyrate_kinase.pdf",
        "plots/butyryl_coa_transferase.pdf",
        "plots/baiH.pdf",
        "plots/baiI.pdf",
        "stats/buty_enzyme_lm.tsv",
        "stats/buty_enzyme_dunn.tsv",
        "stats/bile_enzyme_lm.tsv",
        "stats/bile_enzyme_dunn.tsv",
        "plots/buty_stat_vis.pdf",
        "plots/bile_stat_vis.pdf"



rule picrust2:
    input:
        seqs_fasta = "data/cecal_qiime/fasta_files/dna-sequences.fasta",
        biom = "data/cecal_qiime/total_sum_scaling.biom"
    output:
        out = directory("data/picrust/out_pipeline")
    conda:
        "picrust2"
    shell:
        """
        picrust2_pipeline.py \
            -s {input.seqs_fasta} \
            -i {input.biom} \
            -o {output.out} \
            --stratified \
            --per_sequence_contrib \
            -p 32
        """


rule ko_contrib_filter:
    input:
        ko_in = "data/picrust/out_pipeline/KO_metagenome_out/pred_metagenome_contrib.tsv.gz"
    output:
        ko_out = "data/picrust/meta_contrib.tsv"
    conda:
        "r_env"
    shell:
        """
        Rscript scripts/ko_contrib_filter.R --ko_in {input.ko_in} \
                                            --ko_out {output.ko_out}
        """


rule butyrate_plots:
    input:
        metadata = "data/misc/cecal_processed_metadata.tsv",
        taxonomy = "data/cecal_qiime2/taxonomy.qza",
        ko_contrib = "data/picrust/meta_contrib.tsv"
    output:
        buk = "plots/butyrate_kinase.pdf",
        but = "plots/butyryl_coa_transferase.pdf"
    conda:
        "r_env"
    shell:
        """
        Rscript scripts/butyrate_plots.R --metadata {input.metadata} \
                                         --taxonomy {input.taxonomy} \
                                         --ko {input.ko_contrib} \
                                         --buk_plot {output.buk} \
                                         --but_plot {output.but}
        """


rule bile_acid_plots:
    input:
       metadata = "data/misc/cecal_processed_metadata.tsv",
       taxonomy = "data/cecal_qiime2/taxonomy.qza",
       ko_contrib = "data/picrust/meta_contrib.tsv"
    output:
        baiH = "plots/baiH.pdf",
        baiI = "plots/baiI.pdf"
    conda:
        "r_env"
    shell:
        """
        Rscript scripts/bile_acid_plots.R --metadata {input.metadata} \
                                          --taxonomy {input.taxonomy} \
                                          --ko {input.ko_contrib} \
                                          --baiH_plot {output.baiH} \
                                          --baiI_plot {output.baiI}
        """


rule butyrate_bile_stats:
    input:
       metadata = "data/misc/cecal_processed_metadata.tsv",
       ko_contrib = "data/picrust/meta_contrib.tsv",
       taxonomy = "data/cecal_qiime2/taxonomy.qza"
    output:
        buty_lm = "stats/buty_enzyme_lm.tsv",
        buty_dunn = "stats/buty_enzyme_dunn.tsv",
        bile_lm = "stats/bile_enzyme_lm.tsv",
        bile_dunn = "stats/bile_enzyme_dunn.tsv",
        buty_stat_vis = "plots/buty_stat_vis.pdf",
        bile_stat_vis = "plots/bile_stat_vis.pdf"
    conda:
        "r_env"
    shell:
        """
        Rscript scripts/butyrate_bile_stats.R --metadata {input.metadata} \
                                              --ko {input.ko_contrib} \
                                              --taxonomy {input.taxonomy} \
                                              --buty_lm {output.buty_lm} \
                                              --buty_dunn {output.buty_dunn} \
                                              --bile_lm {output.bile_lm} \
                                              --bile_dunn {output.bile_dunn} \
                                              --buty_stat_vis {output.buty_stat_vis} \
                                              --bile_stat_vis {output.bile_stat_vis}
        """