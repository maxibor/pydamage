import numpy as np


def check_extension(filename):
    extension = filename.split(".")[-1]
    modes = {'bam': 'rb', 'sam': 'r', 'cram': 'rc'}
    try:
        return(modes[extension])
    except KeyError:
        raise Exception(f"{extension} file extension not supported")


def get_x_y(data):
    xdata, counts = np.unique(data, return_counts=True)
    ydata = counts/counts.sum()
    return(xdata, ydata)
