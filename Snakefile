import itertools
import os

configfile: "files/config.yaml"

LIBRARIES = config['project']['libraries'].values()
SPLITS = ['L{}{}'.format('0'*(3 - len(str(n))), n)\
          for n in range(1, len(LIBRARIES) + 1)]
READS = ['R1', 'R2', 'R3', 'R4']

rule all:
    input:
        # expand(temp(os.path.join(config['project']['dir'],
        #                          'procesed', 'fastq', 'weaved',
        #                          '{{split}}_' + '{library}_bc_umi.fastq')),
        #        library=LIBRARIES)  
        expand(os.path.join(config['project']['dir'], 
                            'processed', 'fastq', 'combined',
                            '{library}_cdna.fastq'),
               library=LIBRARIES)

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
        "gzcat {input} > {output}"

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
        expand(temp(os.path.join(config['project']['dir'],
                                 'processed', 'fastq', 'weaved',
                                 '{{split}}_' + '{library}_cdna.fastq')),
               library=LIBRARIES, allow_missing=True),
        expand(temp(os.path.join(config['project']['dir'],
                                 'processed', 'fastq', 'weaved',
                                 '{{split}}_' + '{library}_bc_umi.fastq')),
               library=LIBRARIES, allow_missing=True)      
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