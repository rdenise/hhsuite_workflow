============================
HHsuite workflow report
============================

HHsuite_workflow_ is a Snakemake workflow that produce a graph visualization of the similarities between profiles. This analysis is based on commit version {{ snakemake.config["__workflow_version__"] }}_.

The analysis can be rerun with the following command:

Local:

.. code-block:: shell

   cd {{ snakemake.config["__workflow_workdir__"] }}
{% if snakemake.config["__workflow_basedir__"] == snakemake.config["__workflow_workdir__"] %}
   snakemake -j 1 --use-conda {{ snakemake.config["__config_args__"] }}
{% else %}
   snakemake -j 1 --use-conda {{ snakemake.config["__config_args__"] }} -s {{ snakemake.config["__workflow_basedir_short__"] }}/Snakefile
{% endif %}

.. note::

   Since the workflow is still work in progress, make sure to first 
   run the commands with the `--dry-run` (`-n`) flag to make sure you 
   don't inadvertedly have to regenerate large parts of the results.
   Many workflow dependencies are complex and in particular when 
   running smaller parts of the workflow, unexpected things may 
   happen.  


Workflow summary
----------------

The workflow runs the following steps:

1. Converting the alignment in a3m format
2. Cretating the hhm profiles with hhmake
3. Create the database for hhblits
4. Compare the profiles using hhblits
5. Plotting the network of similarities 


Data organization
-----------------

.. code-block:: text

   {{ snakemake.config["project_name"] }}/                          <- top-level project folder
   │
   │
   ├── report                                <- Folder with the report.html file inside    
   │
   ├── logs                                  <- Collection of log outputs
   │
   ├── databases                             <- Generated analysis database related files
   │    
   ├── alignments                            <- If conversion needed the a3m format of the alignment will be here
   │
   ├── hhblits                               <- Folder with the results of hhblits
   │   ├── hhr                               <- Results of hhblits in hhr format
   │   └── tsv                               <- Results of hhblits in tsv format
   │
   ├── hhm                                   <- HMM profiles in HHM format made with hhmake
   ├── hmm                                   <- HMM profiles in HMM format made with hmmbuild
   │
   └── plots                                 <- Folder of the final results with the figure in pdf and png 



General results
---------------

Figure
******

.. figure:: {{ snakemake.config["__output_folder__"] }}/results/plots/similarities_hhm_profiles.png
   :width: 60%
   :align: center

Workflow graph
--------------


.. _sORTholog: https://github.com/rdenise/hhsuite_workflow
.. _{{ snakemake.config["__workflow_version__"] }}: {{ snakemake.config["__workflow_version_link__"] }}
