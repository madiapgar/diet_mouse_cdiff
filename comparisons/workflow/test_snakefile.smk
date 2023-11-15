rule all:
    input:
        "data/test_core"



rule core_metrics_analysis:
    input:
        tree = "data/comp_qiime/tree.qza",
        meta_filt_tax = "data/comp_qiime/tax_s1_filt.qza",
        metadata = "data/misc/s1_filt_comp_metadata.tsv"
    output:
        out_dir = directory("data/test_core")
    conda:
        "qiime2-2023.5"
    params:
        sampling_depth=99631
    shell:
        """
        qiime diversity core-metrics-phylogenetic \
            --i-phylogeny {input.tree} \
            --i-table {input.meta_filt_tax} \
            --p-sampling-depth {params.sampling_depth} \
            --m-metadata-file {input.metadata} \
            --output-dir {output}
        """