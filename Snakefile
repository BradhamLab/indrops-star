import itertools
import os
from scripts import utils

configfile: "files/2019-10-24_config.yaml"
# shell.prefix("module load miniconda; source activate indrops-star; ")

LIBRARIES = config['project']['libraries'].values()
# LIBRARIES = ['ASW-18hpf']
LANES = [f'L00{n}' for n in range(1, 5)]
# print(LANES)
READS = ['R1', 'R2', 'R3', 'R4']

rule all:
    input:
        expand(os.path.join(config['project']['dir'],
                                'processed', 'STAR', '{library}', 'Solo.out',
                                'Gene', 'raw', 'matrix.mtx'),
                   library=LIBRARIES),
        # expand(os.path.join(config['project']['dir'], 
        #                     'processed', 'fastq', 'demultiplexed',
        #                     '{lane}_{name}_bc_umi.fastq'),
        #        lane=LANES, name=LIBRARIES),
        # os.path.join(config['project']['dir'],
        #                     'processed', 'fastq', 'trimmed',
        #                     'L001_lib_bc_umi.fastq'),
        # expand(os.path.join(config['project']['dir'], 
        #                   'processed', 'fastq', 'combined',
        #                   '{library}_cdna.fastq'),
        #        library=LIBRARIES)
            #    lane=LANES, allow_missing=True)
        # expand(os.path.join(config['project']['dir'],
        #                     'processed', 'STAR', '{library}',
        #                     'Solo.out/Gene/raw/matrix.mtx'),
        #        library=LIBRARIES),
        # os.path.join(config['project']['dir'],
        #              'summaries', 'reads', 'read_counts.png'),
rule extract_fastqs:
    output:
        expand(os.path.join(config['project']['dir'],
                            'Data', 'Intensities', 'BaseCalls',
                            'Undetermined_S0_{lane}_{r}_001.fastq.gz'),
               lane=LANES, r=READS)
    params:
        project_dir=config['project']['dir']
    shell:
        """
        cd {params.project_dir} 
        bcl2fastq --use-bases-mask y*,y*,y*,y* --mask-short-adapter-reads 0 --minimum-trimmed-read-length 0
        """

rule unzip_fastqs:
    input:
        os.path.join(config['project']['dir'],
                      'Data', 'Intensities', 'BaseCalls',
                      'Undetermined_S0_{lane}_R1_001.fastq.gz')
    output:
        temp(os.path.join(config['project']['dir'],
                     'tmp', '{lane}_R1_001.fastq')) #TODO make these temporary again
    shell:
        "gunzip -c {input} > {output}"

# workflow should be unzip (temp) -> combine (temp output) -> trim -> match ->
# create library fasta -> demultiplex -> combine

# interleaf library index, 1/2 barcode, 1/2 barcode + umi reads to single
# fastq file
rule interleaf_r3r2r4:
    input:
        r3=os.path.join(config['project']['dir'],
                      'Data', 'Intensities', 'BaseCalls',
                      'Undetermined_S0_{lane}_R3_001.fastq.gz'),
        r2=os.path.join(config['project']['dir'],
                      'Data', 'Intensities', 'BaseCalls',
                      'Undetermined_S0_{lane}_R2_001.fastq.gz'),
        r4=os.path.join(config['project']['dir'],
                      'Data', 'Intensities', 'BaseCalls',
                      'Undetermined_S0_{lane}_R4_001.fastq.gz'),
    output:
        temp(os.path.join(config['project']['dir'], 'processed', 'fastq', 'interleafed',
                     '{lane}_lib_bc_umi.fastq')) # TODO this should become temporary file

    shell:
        "paste <(zcat {input.r3}) <(zcat {input.r2}) <(zcat {input.r4}) | "
        "awk '{{if (NR%4==1 || NR%4==3) {{print $1}} else "
        "{{print $1 $2 $3}} }}' > {output}"


# quality trim reads using cutadapt b/c it allows nextseq-trim 
rule trim_reads:
    input:
        cdna=os.path.join(config['project']['dir'],
                        'tmp', '{lane}_R1_001.fastq')
    output:
        cdna=temp(os.path.join(config['project']['dir'], 
                               'processed', 'fastq', 'trimmed',
                               '{lane}_cdna.fastq')),
    log:
        os.path.join('logs', 'cutadapt', '{lane}.log')
    params:
        cores=config['params']['cutadapt']['cores'],
        extra=' '.join(f"--{k}={v}" for k,v in \
                       config['params']['cutadapt']['extra'].items())
    shell:
        "(cutadapt --cores {params.cores} --minimum-length 1 "
        "--output {output.cdna} {params.extra} {input.cdna}) 2> {log}"

# paired barcode + umi reads with filtered cdna reads 
rule pair_reads:
    input: 
        cdna=os.path.join(config['project']['dir'], 
                          'processed', 'fastq', 'trimmed',
                          '{lane}_cdna.fastq'),
        bc_umi=os.path.join(config['project']['dir'],
                            'processed', 'fastq', 'interleafed',
                            '{lane}_lib_bc_umi.fastq')
    log:
        os.path.join('logs', 'fastq_pair', '{lane}.log')
    params:
        bc_umi_paired=temp(os.path.join(config['project']['dir'],
                           'processed', 'fastq', 'interleafed',
                           '{lane}_lib_bc_umi.fastq.paired.fq')),
    output:
        cdna_paired=temp(os.path.join(config['project']['dir'], 
                                      'processed', 'fastq', 'trimmed',
                                      '{lane}_cdna.fastq.paired.fq')),
        cdna_single=temp(os.path.join(config['project']['dir'], 
                                      'processed', 'fastq', 'trimmed',
                                      '{lane}_cdna.fastq.single.fq')),
        bc_umi_single=temp(os.path.join(config['project']['dir'],
                                        'processed', 'fastq', 'interleafed',
                                        '{lane}_lib_bc_umi.fastq.single.fq')),
        bc_keep=temp(os.path.join(config['project']['dir'],
                             'processed', 'fastq', 'trimmed',
                             '{lane}_lib_bc_umi.fastq'))
    shell:
        "n_reads=$(awk 'NR % 4 == 2' {input.cdna} | wc -l); "
        "(fastq_pair -t $n_reads {input.cdna} {input.bc_umi}) 2> {log}; "
        "mv {params.bc_umi_paired} {output.bc_keep}"

demult_dir = os.path.join(config['project']['dir'], 'processed', 'fastq', 'demultiplexed')
# gonna need some *cores*
rule cutadapt_demultiplex:
    input:
        lib_bc_umi=os.path.join(config['project']['dir'],
                                'processed', 'fastq', 'trimmed',
                                '{lane}_lib_bc_umi.fastq'),
        cdna=os.path.join(config['project']['dir'], 
                          'processed', 'fastq', 'trimmed',
                          '{lane}_cdna.fastq')
    params:
        barcodes=" ".join(f"-g {v}=^{k}" for k,v in config['project']['libraries'].items()),
        error=config['params']['demultiplex']['error'],
        cores=config['params']['demultiplex']['cores'],
        bc_prefix=lambda wc: os.path.join(demult_dir,
                                          f"{wc.lane}" + "_{name}_bc_umi.fastq"),
        cdna_prefix=lambda wc: os.path.join(demult_dir,
                                            f"{wc.lane}" + "_{name}_cdna.fastq"), 
    output:
        bc=temp(expand(os.path.join(demult_dir, '{{lane}}_{library}_bc_umi.fastq'),
                  library=LIBRARIES)),
        cdna=temp(expand(os.path.join(demult_dir, '{{lane}}_{library}_cdna.fastq'),
                    library=LIBRARIES))
    log:
        os.path.join('logs', 'demultiplix', '{lane}.log')
    shell:
        "cutadapt --cores {params.cores} -e {params.error} --no-indels {params.barcodes} "
        "-o {params.bc_prefix} -p {params.cdna_prefix} "
        "{input.lib_bc_umi} {input.cdna}"


rule combine_fastqs:
    input:
        cdna=expand(os.path.join(config['project']['dir'],
                                 'processed', 'fastq', 'trimmed',
                                 '{lane}_{{library}}_cdna.fastq'),
                    lane=LANES, allow_missing=True),
        bc_umi=expand(os.path.join(config['project']['dir'],
                                   'processed', 'fastq', 'trimmed',
                                   '{lane}_{{library}}_bc_umi.fastq'),
                      lane=LANES, allow_missing=True)
    output:
        cdna=os.path.join(config['project']['dir'], 
                          'processed', 'fastq', 'combined',
                          '{library}_cdna.fastq'),
        bc_umi=os.path.join(config['project']['dir'],
                            'processed', 'fastq', 'combined',
                            '{library}_bc_umi.fastq')
    shell:
        "cat {input.cdna} > {output.cdna}; "
        "cat {input.bc_umi} > {output.bc_umi};"


# rule count_reads_per_barcode:
#     input:
#         bc_umi=os.path.join(config['project']['dir'],
#                             'processed', 'fastq', 'combined',
#                             '{library}_bc_umi.fastq'),
#     output:
#         os.path.join(config['project']['dir'], 'processed', 'reads',
#                      '{library}_reads_per_barcode.csv')
#     shell:
#         "awk 'NR % 4 == 2' {input.bc_umi} | cut -c 1-16 | sort | uniq -c | sort -nr > {output}"

# rule match_and_plot_barcode_reads:
#     input:
#         counts=os.path.join(config['project']['dir'], 'processed', 'reads',
#                             '{library}_reads_per_barcode.csv'),
#         whitelist="ref/gel_barcode3_list.txt"
#     output:
#         csv=os.path.join(config['project']['dir'], 'summaries',
#                          '{library}_reads_per_barcode.csv'),
#         summary=os.path.join(config['project']['dir'], 'summaries',
#                          '{library}_reads_per_barcode_summary.csv'),
#         pdf=os.path.join(config['project']['dir'], 'summaries',
#                          '{library}_read_per_barcode_pdf.png'),
#         cdf=os.path.join(config['project']['dir'], 'summaries',
#                          '{library}_read_per_barcode_cdf.png')
#     script:
#         "scripts/reads_per_barcode.py"

# ## summarize reads per library + ambig
# rule count_reads_per_library:
#     input:
#         cdna=os.path.join(config['project']['dir'], 
#                           'processed', 'fastq', 'combined',
#                           '{library}_cdna.fastq')
#     output:
#         temp(os.path.join(config['project']['dir'],
#                                  'summaries', 'reads',
#                                  '{library}_read_counts.csv'))
#     shell:
#         "echo $(awk {{s++}}END{{print\ s/4}} {input.cdna}),{wildcards.library} > {output}"


# rule count_ambig_reads:
#     input:
#         expand(os.path.join(config['project']['dir'],
#                             'summaries', 'reads',
#                             '{library}_read_counts.csv'), library=LIBRARIES)
#     params:
#         weave=os.path.join(config['project']['dir'],
#                                  'processed', 'fastq', 'weaved')
#     output:
#         temp(os.path.join(config['project']['dir'],
#                                  'summaries', 'reads',
#                                  'ambig_read_counts.csv'))
#     shell:
#         "echo $(awk {{s++}}END{{print\ s/4}} {params.weave}/*cdna.fastq),Ambiguous > {output}"

# rule summarize_library_read_counts:
#     input:
#         libraries=expand(os.path.join(config['project']['dir'],
#                                       'summaries', 'reads',
#                                       '{library}_read_counts.csv'),
#                          library=LIBRARIES),
#         ambig=os.path.join(config['project']['dir'],
#                                  'summaries', 'reads',
#                                  'ambig_read_counts.csv')
#     output:
#         os.path.join(config['project']['dir'],
#                      'summaries', 'reads', 'read_counts.csv')
#     shell:
#         "cat {input.libraries} {input.ambig} > {output}"

# rule plot_library_read_counts:
#     input:
#         csv=os.path.join(config['project']['dir'],
#                          'summaries', 'reads', 'read_counts.csv')
#     output:
#         img=os.path.join(config['project']['dir'],
#                          'summaries', 'reads', 'read_counts.png')
#     script:
#         "scripts/plot_library_read_counts.py"


rule build_star_index:
    input:
        fasta=config['genome']['fasta'],
        gtf=config['genome']['gtf'],
    params:
        index_dir=config['STAR']['index'],
        chr_n_bits=utils.estimate_STAR_ChrBinNbits(config['genome']['fasta'], 60),
        extra=config['params']['star_index']
    log:
        "logs/star_index.log"
    output:
        os.path.join(config['STAR']['index'], 'Genome')
    shell:
        "(STAR --runMode genomeGenerate --genomeDir {params.index_dir} "
        "--genomeFastaFiles {input.fasta} --sjdbGTFfile {input.gtf} "
        "--genomeChrBinNbits {params.chr_n_bits} {params.extra}) > {log}"

rule run_star_solo:
    input:
        index=os.path.join(config['STAR']['index'], 'Genome'),
        cdna=os.path.join(config['project']['dir'], 
                          'processed', 'fastq', 'combined',
                          '{library}_cdna.fastq'),
        bc_umi=os.path.join(config['project']['dir'],
                            'processed', 'fastq', 'combined',
                            '{library}_bc_umi.fastq'),
        whitelist="ref/gel_barcode3_list.txt"
    params:
        index=config['STAR']['index'],
        out=os.path.join(config['project']['dir'], 'processed', 'STAR', '{library}') + '/',
        solo=" ".join(f"--{k} {v}" for k, v in config['params']['star_solo'].items())
    output:
        mtx=os.path.join(config['project']['dir'],
                                'processed', 'STAR', '{library}', 'Solo.out',
                                'Gene', 'raw', 'matrix.mtx'),
        barcodes=os.path.join(config['project']['dir'],
                                     'processed', 'STAR', '{library}', 'Solo.out',
                                     'Gene', 'raw', 'barcodes.tsv'),
        features=os.path.join(config['project']['dir'],
                                     'processed', 'STAR', '{library}', 'Solo.out',
                                     'Gene', 'raw', 'features.tsv')
    log:
        "logs/{library}_starsolo.log"
    shell:
        "(STAR --genomeDir {params.index} "
        "--readFilesIn {input.cdna} {input.bc_umi} --soloType CB_UMI_Simple "
        "--soloCBwhitelist {input.whitelist} --soloFeatures Gene SJ GeneFull "
        "--soloStrand Unstranded [Reverse/Forward] --outFileNamePrefix {params.out} "
        "--soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 6 "
        "{params.solo}) 2> {log}"

