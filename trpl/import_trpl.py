"""import_trpl.py"""

import numpy as np

def import_trpl(csv_file, norm_average):
    
    # Import data.
    with open(csv_file) as f:    
        lines = f.readlines()

    # Parse CSV format.
    data = []
    
    for line in lines:
        data.append([float(item) for item in line.split(',')])

    data = np.array(data)

    # Shift data to t = 0. Time in units of ns.
    t = data[:, 0]
    t -= t[0]

    trpl = data[:, 1:].transpose()
    
    #trpl=np.delete(trpl,1,0)
    #print(trpl)

    # Align different powers to maximum TRPL value.
    maxclip = 0
    for row in trpl:
        maxind = np.argmax(row)
        row[:len(row)-maxind] = row[maxind:]
        maxclip = max(maxclip, maxind)
    
    t = t[:len(t)-maxclip]
    trpl = trpl[:, :trpl.shape[1]-maxclip]

    if norm_average:
        # Normalize data to average of first four time points.
        for row in trpl:
            row /= np.mean(row[:4])

    else:
        # Normalize data to first time point.
        trpl /= trpl[:, [0]]

    return t, trpl
