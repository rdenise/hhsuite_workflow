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
        annotations=annotations_file,
    output:
        png=report(
            os.path.join(
            OUTPUT_FOLDER,
            "plots",
            "similarities_hhm_profiles.png",
            ),
            caption="../report/similarities_png.rst",
            category="Plots",
        ),
        pdf=report(
            os.path.join(
            OUTPUT_FOLDER,
            "plots",
            "similarities_hhm_profiles.pdf",
            ),
            caption="../report/similarities_pdf.rst",
            category="Plots",
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
