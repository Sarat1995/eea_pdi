"""util.py"""

import numpy as np

# Function for formatting scientific notation.
def latex_float(f):
    float_str = "{0:.2g}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"{0} \times 10^{{{1}}}".format(base, int(exponent))
    else:
        return float_str

def time_segments(t, segments):
    # Define time segments for different annihilation constants.    
    time_split = np.array_split(np.arange(len(t)), segments)

    t_segments = []

    for i in range(segments):

        t_segments.append(t[time_split[i][0] : time_split[i][-1] + 2])

    return t_segments

def segments_like(a, b):

    c = b.transpose()

    output = []

    start = 0

    for i in range(len(a)):

        output.append(c[start:start + len(a[i])].transpose())

        start += len(a[i]) - 1

    return output
