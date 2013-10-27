import pylab as pl
import string as st

def importWeightsAsArray(_fname):
    f = open(_fname, 'r')
    arr = []

    for l in f.readlines():
        arr.append( float(st.split(l)[0]) )
    f.close()

    return arr

def calculateDifferenceOfWeigths(_w1_primal, _w2):
    return ( sum( [ abs(_w1_primal[i] - _w2[i]) 
        for i in range(len(_w1_primal)) ] ) / sum(_w1_primal) )

def plotWeigthHistogram(_arr, _fi = [1], _colnum = 10, _comm = ''):
    pl.figure(_fi[0])
    pl.clf()
    pl.hist(_arr, _colnum)
    pl.xlabel('weight')
    pl.ylabel('amount of synapses')
    pl.axis([0, 120, 0, 100])
    pl.savefig('fig' + str(_fi[0]) + _comm + '.png')
    _fi[0] += 1

def getEvolutionOfWeights(_step, _max_files):
    times = [int(i) for i in pl.linspace(0, _step * (_max_files-1), _max_files) ]
    difs = []
    wtar = importWeightsAsArray('../data/weight_target.txt')

    for i in times:
        w = importWeightsAsArray('../data/we/weight_'+str(i)+'.txt')
        difs.append(calculateDifferenceOfWeigths(wtar, w))
    return [times, difs]

def pltSomePlot(_x, _y, _fi = [1], _comm = '', _xlabel = '', _ylabel = ''):
    pl.figure(_fi[0])
    pl.clf()
    pl.plot(_x, _y)
    pl.xlabel(_xlabel)
    pl.ylabel(_ylabel)
    pl.axis([0, max(_x), 0, 1])
    pl.savefig('fig' + str(_fi[0]) + _comm + '.png')

    _fi[0] += 1

def calcAPlotMainPoints(_fig_iter):
    trg = importWeightsAsArray('../data/weight_target.txt')
    stw = importWeightsAsArray('../data/weight_2sim.txt')
    wfi = importWeightsAsArray('../data/weight_final.txt')

    plotWeigthHistogram(trg, _fig_iter, 30, '_target_weights')
    plotWeigthHistogram(stw, _fig_iter, 30, '_new_gamma_weights')
    plotWeigthHistogram(wfi, _fig_iter, 30, '_final_weights')

    print 'Difference between target and start 2: ',
    print calculateDifferenceOfWeigths(trg, stw)
    print 'Difference between target and end 2: ',
    print calculateDifferenceOfWeigths(trg, wfi)

def calcAPlotAllOutput(_fig_iter):
    [ti, di] = getEvolutionOfWeights(100, 5001)
    pltSomePlot(ti, di, _fig_iter, 'change of weight through time',
     'time, ms', 'weight difference from target')
    midi = min(di)
    miti = ti[di.index(midi)]
    print 'Minimum of weight difference %f achieved at %f' %( midi, miti )
    plotWeigthHistogram(
        importWeightsAsArray(
            '../data/we/weight_'+str(miti)+'.txt'), _fig_iter, 30, '_target_weights')

def wa():
    fig_iter = [1]
    calcAPlotMainPoints(fig_iter)
    calcAPlotAllOutput(fig_iter)
    pl.show()
    print '%d plots created' %(fig_iter[0]-1)