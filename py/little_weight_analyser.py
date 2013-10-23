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
    pl.savefig('fig' + str(_fi[0]) + _comm + '.png')
    _fi[0] += 1

def wa():
    fig_iter = [1]
    trg = importWeightsAsArray('../data/weight_1.txt')
    stw = importWeightsAsArray('../data/weight_2.txt')
    cur = importWeightsAsArray('../data/weight_3.txt')

    plotWeigthHistogram(trg, fig_iter, 10, '_starting_weights')
    plotWeigthHistogram(stw, fig_iter, 10, '_new_weights')
    plotWeigthHistogram(cur, fig_iter, 10, '_final_weights')

    print calculateDifferenceOfWeigths(trg, stw)