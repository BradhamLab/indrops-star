import json
import sys
sys.path.append('scripts')
import seq_utils

if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        sequences = (line.strip() for line in open(snakemake.input['whitelist'], 'r'))
        neighborhood = seq_utils.build_sequence_neighborhoods(sequences)
        with open(snakemake.output['neighbors'], 'w') as f:
            json.dump(neighborhood, f)