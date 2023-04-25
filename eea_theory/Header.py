def InitFigure():

    import matplotlib.pyplot as plt

    plt.rc('text.latex', preamble='\\usepackage{amsmath}')
    plt.rc('font', **{'family': 'arial', 'size': 7, 'weight': 'bold'})

    #plt.rc('text',usetex=True)

    plt.rcParams['lines.linewidth'] = 1.25

    plt.rc('xtick',labelsize=7)
    plt.rc('ytick',labelsize=7)
    plt.rc('xtick.major',width=.5)
    plt.rc('ytick.major',width=.5)
    plt.rc('axes',labelsize=7,linewidth=.5,labelweight='bold',titlesize=7,titleweight='bold') #,tick_params=dict(width=.5))
    plt.rc('legend', fontsize=7)

    Colors = dict()

    Colors['B'] = [0,0.4470,0.7410]
    Colors['DB'] = [0,0.2235,0.3705]
    #Colors['LB'] = [0.3010,0.7450,0.9330]
    Colors['LB'] = [0.75,0.9,1]
    Colors['R'] = [0.8500,0.3250,0.0980]
    Colors['LR'] = [0.9250,0.6625,0.5490]
    Colors['DR'] = [0.6350,0.0780,0.1840]
    Colors['Y'] = [0.9290,0.6940,0.1250]
    Colors['P'] = [0.4940,0.1840,0.5560]
    Colors['G'] = [0.4660,0.6740,0.1880]
    Colors['LG'] = [0.7330,0.837,0.5940]
    Colors['DG'] = [0.2330,0.3370,0.0940]

    FSize = dict()
    FSize['Label'] = 16
    FSize['Legend'] = 16

    Vars = dict()
    Vars['Width'] = 3.8
    Vars['DWidth'] = 6.69

    return Colors, FSize, Vars
