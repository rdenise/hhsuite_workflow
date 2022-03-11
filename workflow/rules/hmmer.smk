# Module containing all the hmmer related rules

##########################################################################
##########################################################################


rule hmmbuild:
    input:
        aln=lambda wildcards: os.path.join(
            config["alignments"],
            f"{wildcards.name_aln}.{aln_ext}",
        ),
    output:
        hmm=os.path.join(
            OUTPUT_FOLDER,
            "hmm",
            "{name_aln}.hmm",
        ),
    log:
        os.path.join(OUTPUT_FOLDER,
                     "logs",
                     "hmmbuild",
                     "{name_aln}.log"),
    resources:
        cpus=1,
    conda:
        "../envs/hmmer.yaml"
    threads: 1
    shell:
        """
        hmmbuild --cpu {threads} --informat afa --amino {output.hmm} {input.aln} &> "{log}"
        """


##########################################################################
##########################################################################
