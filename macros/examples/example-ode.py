#************************************
#**  author: Antonella Succurro    **
#**  email:a.succurro[AT]gmail.com **
#**                                **
#**  created:     2016/06/01       **
#**  last modified: 2016/10/27     **
#************************************

import sys
sys.path.append('../../code/python/')
import classSimpleODE as sode
import numpy as np
import matplotlib.pyplot as plt
import optparse
import subprocess, os
import json
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]    
colors_rgb = []
# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
for i in range(len(tableau20)):    
    r, g, b = tableau20[i]
    colors_rgb.append((r / 255., g / 255., b / 255.))

def main():
    '''
    Examples of usage of the SimpleODE class
    '''

    opts, args = options()

    if opts.lotkavolterra:
        runLotkaVolterra()

    if opts.michaelismenten:
        runMichaelisMenten()

    return

def runLotkaVolterra():
    '''
    Define a Lotka-Volterra model for Rabbits and Foxes populations
    This examples uses more human readable synthax that is then translated by the object constructor
    [ a*X[0] -   b*X[0]*X[1] ,
    -c*X[1] + d*b*X[0]*X[1] ]
    '''
    parameters = {'a': 1., 'b': 0.1, 'c': 1.5, 'd': 0.75}
    ykeys = ['RABBITS','FOXES']
    y0 = [10, 5]
    frabbit = 'self.p["a"]*RABBITS - self.p["b"]*RABBITS*FOXES'
    ffoxes = '-self.p["c"]*FOXES + self.p["d"]*self.p["b"]*RABBITS*FOXES'
    for k in ykeys:
        frabbit = frabbit.replace(k, writey(k, ykeys))
        ffoxes = ffoxes.replace(k, writey(k, ykeys))
    functions = {'RABBITS': frabbit, 'FOXES': ffoxes}
    odesys = sode.SimpleODE(functions, parameters, ykeys, y0)
    odesys.integrateODE(maxT=16)
    f1 = plt.figure()
    plt.plot(np.array(odesys.T), odesys.Y['RABBITS'], 'r-', label='Rabbits')
    plt.plot(np.array(odesys.T), odesys.Y['FOXES']  , 'b-', label='Foxes')
    plt.grid()
    plt.legend(loc='best')
    plt.xlabel('time')
    plt.ylabel('population')
    plt.title('Evolution of fox and rabbit populations')
    f1.savefig('plot-example-ode_lotka-volterra.png')

def buildODESys_mm():
    '''
    Define a ODE system for Michaelis Menten enzymatic reactions
    This example builds explicitly the ODE system in the python way (y array)
    '''
    parameters = {"k1" : 10, "k-1": 5, "k2": 4}

    ykeys = ['E', 'S', 'SE', 'P']
    y0 = [0.1, 0.1, 0.0, 0.0]
    dS = '-self.p["k1"]*y[0]*y[1] + self.p["k-1"]*y[2]'
    dE = '-self.p["k1"]*y[0]*y[1] + (self.p["k-1"]+self.p["k2"])*y[2]'
    dSE = 'self.p["k1"]*y[0]*y[1] - (self.p["k-1"]+self.p["k2"])*y[2]'
    dP = 'self.p["k2"]*y[2]'
    functions = {'S' : dS, 'E': dE, 'SE': dSE, 'P': dP}

    return functions, parameters, ykeys, y0

def runMichaelisMenten():
    f1 = plt.figure()
    functions, parameters, ykeys, y0 = buildODESys_mm()
    odesys = run_coop(functions, parameters, ykeys, y0)
    for i,y in enumerate(ykeys):
        plt.plot(np.array(odesys.T), odesys.Y[y], color=colors_rgb[i], label='[%s]' % y)
    ax = plt.gca()
    plt.legend(loc='best')
    plt.grid()
    plt.xlabel('time')
    plt.ylabel('Concentration')
    plt.title('S + E <-> SE -> P + E')
    f1.savefig('plot-example-ode_michaelis-menten.png')

def writey(keyword, keylist):
    return 'y[%d]'%keylist.index(keyword)
    
def run_coop(functions, parameters, ykeys, y0):

    odesys = sode.SimpleODE(functions, parameters, ykeys, y0)
    odesys.integrateODE(maxT=10, stepChoiceLevel=(0, 10, 100), alwayspositive=False)

    return odesys

def options():
    '''define here in-line arguments'''
    parser = optparse.OptionParser(description='Parsing options')
    parser.add_option('-v', '--verbose', dest='verbose', help='increase output verbosity', action='store_true')
    parser.add_option('-l', '--lotkavolterra', dest='lotkavolterra', help='run the lotka volterra simulation', action='store_true')
    parser.add_option('-p', '--plotmode', dest='plotmode', help='plotmode', action='store_true')
    parser.add_option('-m', '--michaelismenten', dest='michaelismenten', help='run the michaelis menten simulation', action='store_true')
    opts, args = parser.parse_args()
    if opts.verbose:
        print "verbosity turned on"
        print opts
        print args
    return opts, args


if __name__=="__main__":
    main()
