# Supplemental File 1

This archive contains the existing Section 7 code and configuration from the local Supplemental_files/scripts directory. Machine-specific path settings have been replaced with project-relative paths and PATH-resolved executables. Analysis logic, parameters, figures, and results were not rerun or changed.

Place scripts/ under a project directory with sibling results/, annotation/, reference/, and bam/ directories. Stage the input data described in the manuscript first. Run manuscript commands from the project root; parenthesized blocks enter subdirectories temporarily. Run Snakemake from scripts/ with snakemake_env activated; python3 and bash are resolved through PATH. Configuration paths are relative to scripts/. The archived NCBI source files are in ../../essay/resources/downloads/ from scripts/, a sibling project resource directory. Preserve this relative layout or edit both source paths in config.section7.yaml and pass matching --gene-info-source and --gene2go-source paths to build_gene2go.py. That script also searches annotation/ under the project root before the sibling resource directory. Newly written manifest paths are relative to the project root. Existing result/provenance files remain unchanged historical records. Standard /usr/bin/env shebangs are portable interpreter launch directives, not machine-specific data paths, and are retained. The code requires the software in the Table of Materials and the Python libraries imported by the scripts, including PyYAML, NumPy, SciPy, and Matplotlib.

protocol_commands.txt contains manuscript command excerpts in reading order for reference, not an automatically executable end-to-end script. Run commands only after completing the corresponding input preparation steps.

standardize_callers.py initializes strand_consistent to an empty value; annotate_candidates.py subsequently derives a gene-strand compatibility annotation. The displayed GO analysis uses strand_mode=all. The separate consistent analysis is available in the code but is not the source of manuscript Figure 5.

The plotting script outputs figure1 through figure4, corresponding to manuscript Figures 2 through 5. Figure 1 is the workflow schematic and is supplied separately. figure_data.tsv uses script-level figure numbers. Original figures and analysis results were not regenerated during packaging.

The repository is https://github.com/djdjdashabi/RNA-seq. This archive fixes the manuscript code snapshot independently of later changes to the repository.
