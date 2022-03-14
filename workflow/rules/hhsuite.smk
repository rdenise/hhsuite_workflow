# Module containing all the hhsuite related rules

##########################################################################
##########################################################################


rule reformat_a3m:
    input:
        aln=lambda wildcards: os.path.join(
            config["alignments"],
            f"{wildcards.aln_file}.{aln_ext}",
        ),
    output:
        a3m=os.path.join(
            OUTPUT_FOLDER,
            "alignments",
            "a3m",
            "{aln_file}.a3m",
        ),
    params:
        format_aln=aln_format,
    log:
        os.path.join(OUTPUT_FOLDER,
                     "logs",
                     "reformat",
                     "{aln_file}.log"),
    resources:
        cpus=1,
    conda:
        "../envs/hhsuite.yaml"
    threads: 1
    shell:
        """
        HHSCRIPTS=$(which hhmake | sed -E 's/bin.+/scripts/')

        $HHSCRIPTS/reformat.pl "{params.format_aln}" a3m "{input.aln}" "{output.a3m}" &> "{log}"
        """


##########################################################################
##########################################################################


rule hhmake:
    input:
        aln=os.path.join(
            FOLDER_ALN,
            "{name_aln}.a3m",
        ),
    output:
        hhm=os.path.join(
            OUTPUT_FOLDER,
            "hhm",
            "{name_aln}.hhm",
        ),
    params:
        input_aln=hhmake_input_aln,
        name="{name_aln}",
        id_max=hhmake_id,
        diff=hhmake_diff,
        cov=hhmake_cov,
        qid=hhmake_qid,
        qsc=hhmake_qsc,
    log:
        os.path.join(OUTPUT_FOLDER,
                     "logs",
                     "hhmake",
                     "{name_aln}.log"),
    resources:
        cpus=1,
    conda:
        "../envs/hhsuite.yaml"
    threads: 1
    shell:
        """
        hhmake -i "{input.aln}" -o "{output.hhm}" -cov "{params.cov}"\
               -id "{params.id_max}" -diff "{params.diff}"\
               -qid "{params.qid}" -qsc "{params.qsc}" -name "{params.name}"\
               -add_cons &> "{log}"
        """


##########################################################################
##########################################################################


rule hhsuite_db:
    input:
        hhms=expand(
            os.path.join(
                OUTPUT_FOLDER,
                "hhm",
                "{name}.hhm",
            ),
            name=PROT,
        ),
    output:
        a3m_fffiles=multiext(
                os.path.join(
                    OUTPUT_FOLDER,
                    "databases",
                    "{database_name}_a3m"),
                    ".ffdata", ".ffindex"
        ),
        cs219_fffiles=multiext(
                os.path.join(
                    OUTPUT_FOLDER,
                    "databases",
                    "{database_name}_cs219"),
                    ".ffdata", ".ffindex"
        ),
        hhm_fffiles=multiext(
                os.path.join(
                    OUTPUT_FOLDER,
                    "databases",
                    "{database_name}_hhm"),
                    ".ffdata", ".ffindex"
        ),
    params:
        format_aln=aln_format,
        glob_hhms=os.path.join(
                OUTPUT_FOLDER,
                "hhm",
                "*.hhm",
        ),
        glob_a3ms=os.path.join(
                FOLDER_ALN,
                "*.a3m",
        ),
        database_name=lambda w, output: output.hhm_fffiles[0].split("_hhm")[0],
        comparison=comparison,
    log:
        os.path.join(OUTPUT_FOLDER,
                     "logs",
                     "hhsuite_db",
                     "{database_name}.log"),
    resources:
        cpus=4,
    conda:
        "../envs/hhsuite.yaml"
    threads: 4
    script:
        "../scripts/hhsuitedb.py"


##########################################################################
##########################################################################


rule hhblits:
    input:
        hhm=os.path.join(
                OUTPUT_FOLDER,
                "hhm",
                "{name_hhm}.hhm",
        ),
        a3m_fffiles=multiext(
                os.path.join(
                    OUTPUT_FOLDER,
                    "databases",
                    f"{config['project_name']}_a3m"),
                    ".ffdata", ".ffindex"
        ),
        cs219_fffiles=multiext(
                os.path.join(
                    OUTPUT_FOLDER,
                    "databases",
                    f"{config['project_name']}_cs219"),
                    ".ffdata", ".ffindex"
        ),        
        hhm_fffiles=multiext(
                os.path.join(
                    OUTPUT_FOLDER,
                    "databases",
                    f"{config['project_name']}_hhm"),
                    ".ffdata", ".ffindex"
        ),        
    output:
        hhr=os.path.join(
            OUTPUT_FOLDER,
            "hhblits",
            "hhr",
            "{name_hhm}.hhr",
        ),
        tsv=os.path.join(
            OUTPUT_FOLDER,
            "hhblits",
            "tsv",
            "{name_hhm}.tsv",
        ),    
    params:
        database_name=lambda w, input: input.hhm_fffiles[0].split("_hhm")[0],
        qid=hhblits_qid,
        e_val=hhblits_e_val,
        cov=hhblits_cov,
    log:
        os.path.join(OUTPUT_FOLDER,
                     "logs",
                     "hhblits",
                     "{name_hhm}.log"),        
    resources:
        cpus=1,
    conda:
        "../envs/hhsuite.yaml"
    threads: 1
    shell:
        """
        hhblits -i "{input.hhm}" -o "{output.hhr}" -d "{params.database_name}"\
                 -maxres 500000 -cov "{params.cov}" -e "{params.e_val}"\
                 -qid "{params.qid}" -cpu {threads} -blasttab {output.tsv} &> "{log}"
        """


##########################################################################
##########################################################################
