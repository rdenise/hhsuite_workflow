##########################################################################
##########################################################################
##
##                                Library
##
##########################################################################
##########################################################################

import os, sys
import pandas as pd
import numpy as np
from snakemake.utils import validate

##########################################################################
##########################################################################
##
##                               Functions
##
##########################################################################
##########################################################################


def get_final_output():
    """
    Generate final output name
    """
    final_output = (
        multiext(
            os.path.join(
            OUTPUT_FOLDER,
            "plots",
            "similarities_hhm_profiles",
            ),
            ".png", ".pdf"
        ),
        expand(
            os.path.join(
                OUTPUT_FOLDER,
                "hmm",
                "{name_aln}.hmm",
            ),
            name_aln=PROT
        ),
    )
    return final_output


##########################################################################


def check_color_seed(annotations_df):
    """
    Infer color if color is not set by the user in the annotations' file
    """

    if "color" not in annotations_df.columns:
        annotations_df['color'] = config["default_values_plot"]["color"]
    else:
        annotations_df.fillna(value={'color':config["default_values_plot"]["color"]}, inplace=True)

    return annotations_df


##########################################################################
##########################################################################
##
##                                Variables
##
##########################################################################
##########################################################################

# Validation of the config.yaml file
validate(config, schema="../schemas/config.schema.yaml")

# path to seeds sheet (TSV format, columns: HMM, color)
annotations_file = config["default_values_plot"]["annotation_table"]

# Validation of the seed file
annotations_dtypes = {
    "HMM": "string",
    "color": "string",
}

annotations_table = pd.read_table(annotations_file, dtype=annotations_dtypes)

# Check color of the seeds
annotations_table = check_color_seed(annotations_table)

##########################################################################
##########################################################################
##
##                        Core configuration
##
##########################################################################
##########################################################################

## Store some workflow metadata
config["__workflow_basedir__"] = workflow.basedir
config["__workflow_basedir_short__"] = os.path.basename(workflow.basedir)
config["__workflow_workdir__"] = os.getcwd()

if workflow.config_args:
    tmp_config_arg = '" '.join(workflow.config_args).replace("=", '="')
    config["__config_args__"] = f' -C {tmp_config_arg}"'
else:
    config["__config_args__"] = ""

with open(os.path.join(workflow.basedir, "../config/VERSION"), "rt") as version:
    url = "https://github.com/rdenise/hhsuite_workflow/releases/tag"
    config["__workflow_version__"] = version.readline()
    config["__workflow_version_link__"] = f"{url}/{config['__workflow_version__']}"


##########################################################################
##########################################################################
##
##                           Options
##
##########################################################################
##########################################################################

# Name your project
project_name = config["project_name"]

# Result folder
OUTPUT_FOLDER = os.path.join(config["output_folder"], project_name)

# Adding to config for report
config["__output_folder__"] = os.path.abspath(OUTPUT_FOLDER)

# hhmake default max id threshold
hhmake_id = config["default_hhmake_options"]["id"]

# hhmake default diff threshold
hhmake_diff = config["default_hhmake_options"]["diff"]

# hhmake default min coverage threshold
hhmake_cov = config["default_hhmake_options"]["cov"]

# hhmake default min id threshold
hhmake_qid = config["default_hhmake_options"]["qid"]

# hhmake default score per column threshold
hhmake_qsc = config["default_hhmake_options"]["qsc"]

# hhmake default alignment format threshold
hhmake_input_aln = config["default_hhmake_options"]["input_aln"]

# hhblits option max evalue
hhblits_e_val = config["default_hhblits_options"]["e_val"]

# hhblits option percentage identity minimum
hhblits_qid = config["default_hhblits_options"]["qid"]

# hhblits option minimum coverage
hhblits_cov = config["default_hhblits_options"]["cov"]

# Option hhsuitedb 
comparison = config["default_hhsuitedb_options"]["comparison"]

# Plots option max evalue
plots_e_val = config["default_values_plot"]["e_val"]

# Plots option percentage identity minimum
plots_qid = config["default_values_plot"]["qid"]

# Plots option minimum coverage
plots_cov = config["default_values_plot"]["cov"]

##########################################################################

# Get the alignment format
aln_format = config["format_aln"][1:] if config["format_aln"].startswith('.') else config["format_aln"]

# Get the alignment file name
PROT, = glob_wildcards(os.path.join(config["alignments"], 
                       "{prot}." + config["aln_ext"]))

# Get extension alignment
aln_ext = config["aln_ext"]

# Get the alignment folder right for hhmake
if aln_format != "a3m":
    FOLDER_ALN = os.path.join(OUTPUT_FOLDER, "alignments", "a3m")
else :
    FOLDER_ALN = config["alignments"]

##########################################################################
