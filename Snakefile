import itertools
import os
from collections import namedtuple

import yaml

from scripts import utils


configfile: "run_config.yaml"


# load in parameter configurations for STAR, STARSolo, demultiplexing, and trimming
parameters = {}
for step, each in config["params"].items():
    with open(each, "r") as f:
        parameters[step] = yaml.load(f, yaml.SafeLoader)


RUNS = [x for x in config["runs"].keys() if x in config['to_run']]                        
LIBRARIES = []
for k, v in config['runs'].items():
    if k in RUNS:
        LIBRARIES += list(v['libraries'].values())
LANES = [f"L00{n}" for n in range(1, 5)]
READS = ["R1", "R2", "R3", "R4"]


rule all:
    input:
        expand(os.path.join(config["outdir"],
                               "{run}",
                               "plots",
                               "alignment.png"),
               run=RUNS),
        expand(os.path.join(config["outdir"],
                               "{run}",
                               "plots",
                               "summary.png"),
               run=RUNS),


rule count_reads:
    input:
        lambda wc : [os.path.join(config['runs'][wc.run]['dir'],
                        "Data",
                        "Intensities",
                        "BaseCalls",
                        "Undetermined_S0_{}".format(lane) + "_{read}_001.fastq.gz")\
                     for lane in LANES]
    output:
        "output/{run}/{read}_counts.csv"
    shell:
        "n_reads=$(zcat {input} | wc -l); "
        'echo "{wildcards.read},$n_reads" > {output}'



def get_extraction_output(wc, lane=None, read=None, expand=False):
    if lane is None:
        lane = '{lane}'
    if read is None:
        read = '{read}'
    out = os.path.join(config['runs'][wc.run]['dir'],
                        "Data",
                        "Intensities",
                        "BaseCalls",
                        "Undetermined_S0_{}_{}_001.fastq.gz".format(lane, read))
    if expand:
        return expand(out, lane=LANES, read=READS)
    return out

rule extract_fastqs:
    input:
        rundir=lambda wc: config['run'][wc.run]['dir']
    params:
        rundir=lambda wc: config['run'][wc.run]['dir'],
    output:
        expand(os.path.join(params.rundir, 'Data', 'Intensities', 'BaseCalls',
                            "Undetermined_S0_{lane}_{read}_001.fastq.gz"),
               lane=LANES, read=READS, allow_missing=True)
    shell:
        """
        cd {input.rundir};
        bcl2fastq --use-bases-mask y*,y*,y*,y* --mask-short-adapter-reads 0 --minimum-trimmed-read-length 0
        """

# if configs['run_to_combine'] is not None:
#     rule combine_runs:
#         input:
#             read=lambda wc: [os.path.join(config['runs'][each]['dir'],
#                                           "Data",
#                                           "Intensities",
#                                           "BaseCalls",
#                                           "Undetermined_S0_{lane}_{read}_001.fastq.gz")\]

# fake_wc = 
# def runs_to_combine(wc):
#     fake_wc = namedtuple('fake_wc', 'run')

rule unzip_cdna_fastqs:
    input:
        lambda wc: get_extraction_output(wc, read='R1')
    output:
        temp(os.path.join(config['outdir'],
                          "{run}",
                          "tmp",
                          "{lane}_R1_001.fastq")),
    shell:
        "gunzip -c {input} > {output}"


# interleaf library index, 1/2 barcode, 1/2 barcode + umi reads to single
# fastq file
rule interleaf_r3r2r4:
    input:
        r3=lambda wc: get_extraction_output(wc, read='R3'),
        r2=lambda wc: get_extraction_output(wc, read='R2'),
        r4=lambda wc: get_extraction_output(wc, read='R4'),
    output:
        temp(
            os.path.join(
                config["outdir"],
                "{run}",
                "fastq",
                "interleafed",
                "{lane}_lib_bc_umi.fastq",
            )
        ),
    shell:
        "paste <(zcat {input.r3}) <(zcat {input.r2}) <(zcat {input.r4}) | "
        "awk '{{if (NR%4==1 || NR%4==3) {{print $1}} else "
        "{{print $1 $2 $3}} }}' > {output}"


# combine, demultiplex, trim, pair, run
rule combine_fastqs:
    input:
        cdna=expand(os.path.join(config['outdir'],
                                 "{run}",
                                 "tmp",
                                 "{lane}_R1_001.fastq"),
                    lane=LANES,
                    allow_missing=True),
        bc_umi=expand(os.path.join(config["outdir"],
                                   "{run}",
                                   "fastq",
                                   "interleafed",
                                   "{lane}_lib_bc_umi.fastq"),
                      lane=LANES,
                      allow_missing=True),
    output:
        cdna=temp(os.path.join(config["outdir"],
                               "{run}",
                               "fastq",
                               "combined",
                               "cdna.fastq")),
        bc_umi=temp(os.path.join(config["outdir"],
                                 "{run}",
                                 "fastq",
                                 "combined",
                                 "lib_bc_umi.fastq")),
    shell:
        "cat {input.cdna} > {output.cdna}; "
        "cat {input.bc_umi} > {output.bc_umi};"


checkpoint demultiplex_libraries:
    input:
        cdna=os.path.join(
            config["outdir"],
            "{run}",
            "fastq",
            "combined",
            "cdna.fastq",
        ),
        lib_bc_umi=os.path.join(
            config["outdir"],
            "{run}",
            "fastq",
            "combined",
            "lib_bc_umi.fastq",
        ),
        config="files/configs/demultiplex.yaml",
    params:
        barcodes=lambda wc: " ".join(f"-g {v}=^{k}" for k, v in config['runs'][wc.run]["libraries"].items()),
        # barcodes='',
        error=parameters["demultiplex"]["error"],
        cores=parameters["demultiplex"]["cores"],
        bc_prefix=lambda wc: os.path.join(config["outdir"],
                               wc.run,
                               "fastq",
                               "demultiplexed",
                               "{name}_bc_umi.fastq"),
        cdna_prefix=lambda wc: os.path.join(config["outdir"],
                                 wc.run,
                                 "fastq",
                                 "demultiplexed",
                                 "{name}_cdna.fastq"),
    output:
        directory(os.path.join(config["outdir"], "{run}", "fastq",
                               "demultiplexed"))
    log:
        os.path.join("logs", "{run}", "demultiplex", "{run}.log"),
    shell:
        "mkdir {output}; "
        "(cutadapt --cores {params.cores} -e {params.error} "
        "--no-indels {params.barcodes} -o {params.bc_prefix} "
        "-p {params.cdna_prefix} {input.lib_bc_umi} {input.cdna} -q 0) > {log}" # do not actually filter reads


# quality trim reads using cutadapt b/c it allows nextseq-trim
rule trim_reads:
    input:
        cdna=os.path.join(config['outdir'],
                          "{run}",
                          "fastq",
                          "demultiplexed",
                          "{library}_cdna.fastq"),
        config="files/configs/cutadapt.yaml",
    output:
        cdna=os.path.join(
                config["outdir"],
                "{run}",
                "fastq",
                "trimmed",
                "{library}_cdna.fastq")
    log:
        os.path.join("logs", "{run}", "cutadapt", "{library}.log"),
    params:
        extra=" ".join(f"--{k}={v}" for k, v in parameters["cutadapt"].items()),
    shell:
        "(cutadapt --minimum-length 1 {params.extra} "
        "--output {output.cdna} {input.cdna}) > {log}"


# paired barcode + umi reads with filtered cdna reads
rule pair_reads:
    input:
        cdna=os.path.join(
            config["outdir"],
            "{run}",
            "fastq",
            "trimmed",
            "{library}_cdna.fastq",
        ),
        bc_umi=os.path.join(
            config["outdir"],
            "{run}",
            "fastq",
            "demultiplexed",
            "{library}_bc_umi.fastq",
        ),
    log:
        os.path.join("logs", "{run}", "fastq_pair", "{library}.log"),
    output:
        cdna_paired=temp(
            os.path.join(
                config["outdir"],
                "{run}",
                "fastq",
                "trimmed",
                "{library}_cdna.fastq.paired.fq",
            )
        ),
        cdna_single=temp(
            os.path.join(
                config["outdir"],
                "{run}",
                "fastq",
                "trimmed",
                "{library}_cdna.fastq.single.fq",
            )
        ),
        bc_umi_single=temp(
            os.path.join(
                config["outdir"],
                "{run}",
                "fastq",
                "demultiplexed",
                "{library}_bc_umi.fastq.single.fq",
            )
        ),
        bc_umi_paired=temp(os.path.join(
                config["outdir"],
                "{run}",
                "fastq",
                "demultiplexed",
                "{library}_bc_umi.fastq.paired.fq",
        ))
    shell:
        "n_reads=$(awk 'NR % 4 == 2' {input.bc_umi} | wc -l); "
        "(fastq_pair -t $n_reads {input.cdna} {input.bc_umi}) > {log};"

rule collect_pairs:
    input:
        cdna=os.path.join(
            config["outdir"],
            "{run}",
            "fastq",
            "trimmed",
            "{library}_cdna.fastq",
        ),
        bc_umi=os.path.join(
                config["outdir"],
                "{run}",
                "fastq",
                "demultiplexed",
                "{library}_bc_umi.fastq.paired.fq",
        ),
    output:
        cdna=os.path.join(config['outdir'],
                          '{run}',
                          'fastq',
                          'paired',
                          '{library}_cdna.fastq.gz'),
        bc_umi=os.path.join(config['outdir'],
                          '{run}',
                          'fastq',
                          'paired',
                          '{library}_bc_umi.fastq.gz'),
    shell:
        "gzip {input.cdna} -c > {output.cdna}; "
        "gzip {input.bc_umi} -c > {output.bc_umi}; "
        

rule build_star_index:
    input:
        fasta=config["genome"]["fasta"],
        gtf=config["genome"]["gtf"],
    params:
        index_dir=config["STAR"]["index"],
        chr_n_bits=utils.estimate_STAR_ChrBinNbits(config["genome"]["fasta"], 60),
        extra=" ".join(f"--{k} {v}" for k, v in parameters["star_index"].items()),
    log:
        "logs/star_index.log",
    output:
        protected(os.path.join(config["STAR"]["index"], "Genome")),
    shell:
        "(STAR --runMode genomeGenerate --genomeDir {params.index_dir} "
        "--genomeFastaFiles {input.fasta} --sjdbGTFfile {input.gtf} "
        "--genomeChrBinNbits {params.chr_n_bits} {params.extra}) > {log}"


rule run_star_solo:
    input:
        index=os.path.join(config["STAR"]["index"], "Genome"),
        cdna=os.path.join(
            config["outdir"],
            "{run}",
            "fastq",
            "paired",
            "{library}_cdna.fastq.gz",
        ),
        bc_umi=os.path.join(
            config["outdir"],
            "{run}",
            "fastq",
            "paired",
            "{library}_bc_umi.fastq.gz",
        ),
        whitelist="ref/gel_barcode3_list.txt",
        config="files/configs/star_solo.yaml",
    params:
        index=config["STAR"]["index"],
        tmpdir=config['STAR']['tmpdir'],
        out=(
            os.path.join(config["outdir"], "{run}", "STAR", "{library}")
            + "/"
        ),
        solo=" ".join(f"--{k} {v}" for k, v in parameters["star_solo"].items())
    output:
        logfile=os.path.join(config["outdir"],
                             "{run}",
                             "STAR",
                             "{library}",
                             "Log.final.out"),
        barcode_stats=os.path.join(config["outdir"],
                                   "{run}",
                                   "STAR",
                                   "{library}",
                                   "Solo.out",
                                   "Barcodes.stats"),
        gene_summary=os.path.join(config["outdir"],
                                  "{run}",
                                  "STAR",
                                  "{library}",
                                  "Solo.out",
                                  "Gene",
                                  "Summary.csv"),
        gene_stats=os.path.join(config["outdir"],
                                "{run}",
                                "STAR",
                                "{library}",
                                "Solo.out",
                                "Gene",
                                "Features.stats"),
        sam=temp(
            os.path.join(
                config["outdir"],
                "{run}",
                "STAR",
                "{library}",
                "Aligned.out.sam"
            )
        ),
        mtx=os.path.join(
            config["outdir"],
            "{run}",
            "STAR",
            "{library}",
            "Solo.out",
            "Gene",
            "raw",
            "matrix.mtx",
        ),
        barcodes=os.path.join(
            config["outdir"],
            "{run}",
            "STAR",
            "{library}",
            "Solo.out",
            "Gene",
            "raw",
            "barcodes.tsv",
        ),
        features=os.path.join(
            config["outdir"],
            "{run}",
            "STAR",
            "{library}",
            "Solo.out",
            "Gene",
            "raw",
            "features.tsv",
        ),
    log:
        "logs/{run}/star_solo/{library}.log",
    shell:
        "tmpdir={params.tmpdir}/`uuidgen`; "
        "(STAR --genomeDir {params.index} --readFilesCommand zcat "
        "--readFilesIn {input.cdna} {input.bc_umi} --soloType CB_UMI_Simple "
        "--soloCBwhitelist {input.whitelist} --soloFeatures Gene SJ GeneFull Velocyto "
        "--soloStrand Unstranded [Reverse/Forward] --outFileNamePrefix {params.out} "
        "--soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 6 "
        "--outTmpDir $tmpdir {params.solo}) 2> {log}"


def aggregate_starsolo_output(wildcards, out='mtx'):
    checkpoint_output = checkpoints.demultiplex_libraries\
                                   .get(**wildcards)\
                                   .output[0]
    libraries = glob_wildcards(os.path.join(checkpoint_output,
                                            "{library}_cdna.fastq")).library
    if out == 'mtx': 
        path = os.path.join(config["outdir"],
                            "{run}",
                            "STAR",
                            "{library}",
                            "Solo.out",
                            "Gene",
                            "raw",
                            "matrix.mtx")
    elif out == 'log':
        path = os.path.join(config["outdir"],
                            "{run}",
                            "STAR",
                            "{library}",
                            "Log.final.out")
    elif out == 'barcode':
        path = os.path.join(config["outdir"],
                            "{run}",
                            "STAR",
                            "{library}",
                            "Solo.out",
                            "Barcodes.stats")
    elif out == 'summary':
        path = os.path.join(config["outdir"],
                            "{run}",
                            "STAR",
                            "{library}",
                            "Solo.out",
                            "Gene",
                            "Summary.csv")
    return expand(path,
                  run=wildcards.run,
                  library=[x for x in libraries if x != 'unknown'])

rule finished_run:
    input:
        aggregate_starsolo_output
    output:
        'finished/{run}.txt'
    shell:
        "touch {output}"


rule plot_star_solo_results:
    input:
        star=lambda wc: aggregate_starsolo_output(wc, 'log'),
        summary=lambda wc: aggregate_starsolo_output(wc, 'summary'),
    params:
        plotdir=os.path.join(config['outdir'], '{run}', 'plots')
    output:
        alignment=report(os.path.join(config["outdir"],
                               "{run}",
                               "plots",
                               "alignment.png"),
                         caption="report/alignment.rst",
                         category="Alignment",
                         subcategory="{run}",
        ),
        summary=report(os.path.join(config["outdir"],
                             "{run}",
                             "plots",
                             "summary.png"),
                       caption="report/quantification.rst",
                       category="STARSolo",
                       subcategory="{run}",
        )
    script:
        "scripts/plot_starsolo.py"    


    
