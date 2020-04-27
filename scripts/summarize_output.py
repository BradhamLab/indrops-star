import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

def plot_umis(umi_file):
    umis = []
    with open(umi_file, 'r') as f:
        for line in f:
            umis.append(int(line.strip()))