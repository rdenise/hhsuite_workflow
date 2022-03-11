# Rule for plotting figures

##########################################################################
##########################################################################


rule plot_similarities:
    input:
        tsv=expand(
            os.path.join(
                OUTPUT_FOLDER,
                "hhblits",
                "tsv",
                "{name}.tsv",
            ),
            name=PROT,
        ),
    output:
        png=os.path.join(
            OUTPUT_FOLDER,
            "plots",
            "similarities_hhm_profiles.png",
        ),
        pdf=os.path.join(
            OUTPUT_FOLDER,
            "plots",
            "similarities_hhm_profiles.pdf",
        ),
    params:
        e_val = plots_e_val,
        cov = plots_cov,
        pid = plots_qid,
        border_color = config['default_values_plot']['colored_border']
    log:
        os.path.join(OUTPUT_FOLDER, "logs", "plots", "similarities_hhblits.log"),
    conda:
        "../envs/pandas_plots.yaml"
    script:
        "../scripts/similarities_hhblits.py"


##########################################################################
##########################################################################