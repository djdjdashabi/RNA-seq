configfile: "config.section7.yaml"

import shutil

BASH = shutil.which("bash")
if BASH is None:
    raise RuntimeError("Activate an environment with bash available on PATH.")
shell.executable(BASH)

RUNS = [str(value) for value in config["runs"]]
SAMPLE_IDS = {str(key): value for key, value in config["sample_ids"].items()}
CALLERS = ["reditools", "sprint", "jacusa2"]
EXPECTED_ARGS = " ".join(
    f"--expected {SAMPLE_IDS[run]}:{caller}"
    for run in RUNS
    for caller in CALLERS
)
PYTHON = config["python_executable"]
GTF = config["references"]["gtf"]
ALU_BED = config["references"]["alu_bed"]
GENE_INFO_SOURCE = config["references"]["gene_info_source"]
GENE2GO_SOURCE = config["references"]["gene2go_source"]
GENE_INFO = config["references"]["gene_info_human"]
GENE2GO = config["references"]["gene2go_human"]
GENE2GO_MANIFEST = config["references"]["gene2go_manifest"]
GO_FDR = float(config["parameters"]["go_fdr"])

NATIVE_PATTERNS = {
    "reditools": "../results/redi_{run}.out",
    "sprint": "../results/sprint_{run}.out",
    "jacusa2": "../results/jacusa_{run}.out",
}

STANDARDIZED = expand(
    "../results/standardized/{run}/{caller}.tsv",
    run=RUNS,
    caller=CALLERS,
)
ANNOTATED = expand("../results/annotation/{run}.annotated.tsv", run=RUNS)
GENE_LINKS = expand("../results/annotation/{run}.gene_links.tsv", run=RUNS)
FIGURES = expand(
    "../results/figures/figure{number}.{extension}",
    number=range(1, 5),
    extension=["png", "pdf"],
)

wildcard_constraints:
    run="69|70",
    caller="reditools|sprint|jacusa2"


def native_input(wildcards):
    return NATIVE_PATTERNS[wildcards.caller].format(run=wildcards.run)


rule all:
    input:
        STANDARDIZED,
        ANNOTATED,
        "../results/site_gene_links.tsv",
        "../results/candidate_summary.tsv",
        "../results/tool_overlap.tsv",
        "../results/site_support.tsv",
        "../results/region_summary.tsv",
        "../results/alu_summary.tsv",
        "../results/go_enrichment.all.tsv",
        "../results/go_enrichment.consistent.tsv",
        "../results/canonical_summary.tsv",
        "../results/figures/figure_data.tsv",
        FIGURES,


rule standardize:
    input:
        STANDARDIZED


rule summaries:
    input:
        "../results/candidate_summary.tsv",
        "../results/tool_overlap.tsv",
        "../results/site_support.tsv",
        "../results/region_summary.tsv",
        "../results/alu_summary.tsv",
        "../results/go_enrichment.all.tsv",
        "../results/go_enrichment.consistent.tsv",
        "../results/canonical_summary.tsv",


rule figures:
    input:
        FIGURES,
        "../results/figures/figure_data.tsv",


rule build_gene2go:
    input:
        gene_info=GENE_INFO_SOURCE,
        gene2go=GENE2GO_SOURCE,
        code="build_gene2go.py",
    output:
        gene_info=GENE_INFO,
        gene2go=GENE2GO,
        manifest=GENE2GO_MANIFEST,
    threads: 1
    shell:
        """
        {PYTHON} build_gene2go.py \
          --gene-info-source {input.gene_info} --gene2go-source {input.gene2go} \
          --gene-info-output {output.gene_info} --gene2go-output {output.gene2go} \
          --manifest {output.manifest}
        """


rule standardize_caller:
    input:
        native=native_input,
        code="standardize_callers.py",
    output:
        "../results/standardized/{run}/{caller}.tsv"
    params:
        sample_id=lambda wildcards: SAMPLE_IDS[wildcards.run]
    threads: 1
    shell:
        """
        {PYTHON} standardize_callers.py --caller {wildcards.caller} \
          --sample {params.sample_id} --input {input.native} --output {output}
        """


rule combine_run_candidates:
    input:
        candidates=lambda wildcards: expand(
            "../results/standardized/{run}/{caller}.tsv",
            run=wildcards.run,
            caller=CALLERS,
        ),
        code="concat_tsv.py",
    output:
        "../results/standardized/{run}/all_callers.tsv"
    params:
        inputs=lambda wildcards, input: " ".join(
            f"--input {path}" for path in input.candidates
        )
    threads: 1
    shell:
        "{PYTHON} concat_tsv.py {params.inputs} --output {output}"


rule annotate_run:
    input:
        candidates="../results/standardized/{run}/all_callers.tsv",
        gtf=GTF,
        code="annotate_candidates.py",
    output:
        sites="../results/annotation/{run}.annotated.tsv",
        links="../results/annotation/{run}.gene_links.tsv",
    threads: 8
    resources:
        mem_mb=8000
    shell:
        """
        {PYTHON} annotate_candidates.py --input {input.candidates} --gtf {input.gtf} \
          --sites-output {output.sites} --links-output {output.links}
        """


rule summarize_candidates:
    input:
        candidates=STANDARDIZED,
        code="summarize_candidates.py",
    output:
        summary="../results/candidate_summary.tsv",
        overlap="../results/tool_overlap.tsv",
        support="../results/site_support.tsv",
    params:
        inputs=lambda wildcards, input: " ".join(
            f"--input {path}" for path in input.candidates
        ),
        expected=EXPECTED_ARGS,
    threads: 1
    shell:
        """
        {PYTHON} summarize_candidates.py {params.inputs} {params.expected} --summary {output.summary} \
          --overlap {output.overlap} --support {output.support}
        """


rule summarize_regions:
    input:
        candidates=ANNOTATED,
        code="summarize_regions.py",
    output:
        "../results/region_summary.tsv"
    params:
        inputs=lambda wildcards, input: " ".join(
            f"--input {path}" for path in input.candidates
        ),
        expected=EXPECTED_ARGS,
    threads: 1
    shell:
        "{PYTHON} summarize_regions.py {params.inputs} {params.expected} --output {output}"


rule alu_overlap:
    input:
        candidates=STANDARDIZED,
        alu=ALU_BED,
        code="calculate_alu_overlap.py",
    output:
        summary="../results/alu_summary.tsv",
        membership="../results/alu_membership.tsv",
    params:
        inputs=lambda wildcards, input: " ".join(f"--input {path}" for path in input.candidates),
        expected=EXPECTED_ARGS,
    threads: 1
    shell:
        """
        {PYTHON} calculate_alu_overlap.py {params.inputs} {params.expected} --alu-bed {input.alu} \
          --summary {output.summary} --membership {output.membership}
        """


rule combine_gene_links:
    input:
        links=GENE_LINKS,
        code="concat_tsv.py",
    output:
        "../results/site_gene_links.tsv"
    params:
        inputs=lambda wildcards, input: " ".join(f"--input {path}" for path in input.links)
    threads: 1
    shell:
        "{PYTHON} concat_tsv.py {params.inputs} --output {output}"


rule go_enrichment:
    input:
        links="../results/site_gene_links.tsv",
        gene_info=GENE_INFO,
        gene2go=GENE2GO,
        code="go_enrichment.py",
    output:
        all_go="../results/go_enrichment.all.tsv",
        all_map="../results/go_mapping.all.tsv",
        consistent_go="../results/go_enrichment.consistent.tsv",
        consistent_map="../results/go_mapping.consistent.tsv",
    threads: 1
    shell:
        """
        {PYTHON} go_enrichment.py --links {input.links} --gene-info {input.gene_info} \
          --gene2go {input.gene2go} --strand-mode all --fdr {GO_FDR} \
          --enrichment-output {output.all_go} --mapping-output {output.all_map}
        {PYTHON} go_enrichment.py --links {input.links} --gene-info {input.gene_info} \
          --gene2go {input.gene2go} --strand-mode consistent --fdr {GO_FDR} \
          --enrichment-output {output.consistent_go} --mapping-output {output.consistent_map}
        """


rule canonical_summary:
    input:
        candidates="../results/candidate_summary.tsv",
        overlap="../results/tool_overlap.tsv",
        support="../results/site_support.tsv",
        regions="../results/region_summary.tsv",
        alu="../results/alu_summary.tsv",
        go_all="../results/go_enrichment.all.tsv",
        go_consistent="../results/go_enrichment.consistent.tsv",
        code="assemble_canonical_summary.py",
    output:
        "../results/canonical_summary.tsv"
    threads: 1
    shell:
        """
        {PYTHON} assemble_canonical_summary.py --candidate-summary {input.candidates} \
          --overlap {input.overlap} --support {input.support} --alu {input.alu} \
          --regions {input.regions} \
          --go {input.go_all} --go {input.go_consistent} --output {output}
        """


rule plot_figures:
    input:
        summary="../results/canonical_summary.tsv",
        code="plot_results.py",
    output:
        figures=FIGURES,
        data="../results/figures/figure_data.tsv",
    threads: 1
    shell:
        """
        {PYTHON} plot_results.py --summary {input.summary} --output-dir ../results/figures \
          --data-output {output.data}
        """
