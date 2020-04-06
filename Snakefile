import itertools
import os
from scripts import utils

configfile: "files/config.yaml"
shell.prefix("source activate indrops-star; ")

LIBRARIES = config['project']['libraries'].values()
SPLITS = ['L{}{}'.format('0'*(3 - len(str(n))), n)\
          for n in range(1, len(LIBRARIES) + 1)]
READS = ['R1', 'R2', 'R3', 'R4']

rule all:
    input:
        expand(os.path.join(config['project']['dir'],
                            'processed', 'STAR', '{library}',
                            'Solo.out/Gene/raw/matrix.mtx'),
               library=LIBRARIES),
        os.path.join(config['project']['dir'],
                     'summaries', 'reads', 'read_counts.png')

rule extract_fastqs:
    output:
        expand(os.path.join(config['project']['dir'],
                            'Data', 'Intensities', 'BaseCalls',
                            'Undetermined_S0_{split}_{r}_001.fastq.gz'),
               split=SPLITS, r=READS)
    params:
        project_dir=config['project']['dir']
    shell:
        """
        cd {params.project_dir} 
        module load bcl2fastq
        bcl2fastq --use-bases-mask y*,y*,y*,y* --mask-short-adapter-reads 0 --minimum-trimmed-read-length 0
        """

rule unzip_fastqs:
    input:
        os.path.join(config['project']['dir'],
                      'Data', 'Intensities', 'BaseCalls',
                      'Undetermined_S0_{split}_{r}_001.fastq.gz')
    output:
        temp(os.path.join(config['project']['dir'],
                          'tmp', '{split}_{r}_001.fastq'))
    shell:
        "gunzip -c {input} > {output}"

rule weave_fastqs:
    input:
        r1=os.path.join(config['project']['dir'],
                        'tmp', '{split}_R1_001.fastq'),
        r2=os.path.join(config['project']['dir'],
                        'tmp', '{split}_R2_001.fastq'),
        r3=os.path.join(config['project']['dir'],
                        'tmp', '{split}_R3_001.fastq'),
        r4=os.path.join(config['project']['dir'],
                        'tmp', '{split}_R4_001.fastq')
    params:
        prefix= os.path.join(config['project']['dir'],
                             'processed', 'fastq', 'weaved',
                             '{split}_'),
        libraries=config['project']['libraries'],
        nsubs=config['params']['weave_fastqs']['mismatches']
    output:
        temp(expand(os.path.join(config['project']['dir'],
                                 'processed', 'fastq', 'weaved',
                                 '{{split}}_' + '{library}_cdna.fastq'),
               library=LIBRARIES, allow_missing=True)),
        temp(expand(os.path.join(config['project']['dir'],
                                 'processed', 'fastq', 'weaved',
                                 '{{split}}_' + '{library}_bc_umi.fastq'),
               library=LIBRARIES, allow_missing=True))
    script:
        "scripts/weave_fastqs.py"

rule combine_fastqs:
    input:
        cdna=expand(os.path.join(config['project']['dir'],
                                 'processed', 'fastq', 'weaved',
                                 '{split}_{{library}}_cdna.fastq'),
                    split=SPLITS, allow_missing=True),
        bc_umi=expand(os.path.join(config['project']['dir'],
                                   'processed', 'fastq', 'weaved',
                                   '{split}_{{library}}_bc_umi.fastq'),
                      split=SPLITS, allow_missing=True)
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

## summarize reads per library + ambig
rule count_reads_per_library:
    input:
        cdna=os.path.join(config['project']['dir'], 
                          'processed', 'fastq', 'combined',
                          '{library}_cdna.fastq')
    output:
        temp(os.path.join(config['project']['dir'],
                                 'summaries', 'reads',
                                 '{library}_read_counts.csv'))
    shell:
        "echo $(awk {{s++}}END{{print\ s/4}} {input.cdna}),{wildcards.library} > {output}"

rule count_ambig_reads:
    input:
        expand(os.path.join(config['project']['dir'],
                            'summaries', 'reads',
                            '{library}_read_counts.csv'), library=LIBRARIES)
    params:
        weave=os.path.join(config['project']['dir'],
                                 'processed', 'fastq', 'weaved')
    output:
        temp(os.path.join(config['project']['dir'],
                                 'summaries', 'reads',
                                 'ambig_read_counts.csv'))
    shell:
        "echo $(awk {{s++}}END{{print\ s/4}} {params.weave}/*cdna.fastq),Ambiguous > {output}"

rule summarize_library_read_counts:
    input:
        libraries=expand(os.path.join(config['project']['dir'],
                                      'summaries', 'reads',
                                      '{library}_read_counts.csv'),
                         library=LIBRARIES),
        ambig=os.path.join(config['project']['dir'],
                                 'summaries', 'reads',
                                 'ambig_read_counts.csv')
    output:
        os.path.join(config['project']['dir'],
                     'summaries', 'reads', 'read_counts.csv')
    shell:
        "cat {input.libraries} {input.ambig} > {output}"

rule plot_library_read_counts:
    input:
        csv=os.path.join(config['project']['dir'],
                         'summaries', 'reads', 'read_counts.csv')
    output:
        img=os.path.join(config['project']['dir'],
                         'summaries', 'reads', 'read_counts.png')
    script:
        "scripts/plot_library_read_counts.py"

rule build_star_index:
    input:
        fasta=config['genome']['fasta'],
        gtf=config['genome']['gtf'],
    params:
        index_dir=config['STAR']['index'],
        chr_n_bits=utils.estimate_STAR_ChrBinNbits(config['genome']['fasta'], 60),
        indexNBases=13
    log:
        "logs/star_index.log"
    output:
        os.path.join(config['STAR']['index'], 'Genome')
    shell:
        "(STAR --runMode genomeGenerate --genomeDir {params.index_dir} "
        "--genomeFastaFiles {input.fasta} --sjdbGTFfile {input.gtf} "
        "--sjdbOverhang 60 --genomeChrBinNbits {params.chr_n_bits} "
        "--sjdbGTFtagExonParentTranscript Parent "
        "--genomeSAindexNbases {params.indexNBases}) 2> {log}"

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
        out=os.path.join(config['project']['dir'],
                         'processed', 'STAR', '{library}')
    output:
        mtx=os.path.join(config['project']['dir'],
                         'processed', 'STAR', '{library}',
                         'Solo.out/Gene/raw/matrix.mtx'),
        barcodes=os.path.join(config['project']['dir'],
                              'processed', 'STAR', '{library}',
                              'Solo.out/Gene/raw/barcodes.tsv'),
        features=os.path.join(config['project']['dir'],
                              'processed', 'STAR', '{library}',
                              'Solo.out/Gene/raw/features.tsv')
    log:
        "logs/{library}_starsolo.log"
    shell:
        "(STAR --genomeDir {params.index} "
        "--readFilesIn {input.cdna} {input.bc_umi} --soloType CB_UMI_Simple "
        "--soloCBwhitelist {input.whitelist} --soloFeatures Gene SJ GeneFull "
        "--soloStrand Unstranded Forward --outFileNamePrefix {params.out} "
        "--soloCBstart 1 --soloCBlen 16 --soloUMIstart 18 --soloUMIlen 6) 2> {log}"