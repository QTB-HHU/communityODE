#************************************
#**  author: Antonella Succurro    **
#**  email:a.succurro[AT]gmail.com **
#**                                **
#**  created:     2016/06/01       **
#**  last modified: 2016/10/27     **
#************************************

import sys
sys.path.append('../../code/python/')
import numpy as np
import matplotlib.pyplot as plt
import optparse
import subprocess, os
import json
import os.path
import matplotlib as mpl
mpl.rcParams['legend.numpoints'] = 1

import classSimpleODE as sode
from objectivesAndRanges import *
from mycolors import *

orgkeys = ['DIATOM', 'PSEUDOALTERO', 'FLAVO', 'ALTEROMO', 'PSEUDOMO']
orgkeysshort = {'DIATOM':'D', 'PSEUDOALTERO':'PA', 'FLAVO':'F', 'ALTEROMO':'A', 'PSEUDOMO':'P'}

def main():
    '''
    Build and solve an ODE model with the options passed through the parser
    See available options under the options() function
    '''

    opts, args = options()
    jsondicts = opts.jsondicts.split(' ')
    functions, parameters, ykeys, y0 = buildODESys(jsondicts[0], opts.fullmedia, opts.minimalmedia, opts.diatom, opts.pseudoaltero, opts.altero, opts.flavo, opts.pseudo, opts.verbose)
    odesys = run_coop(functions, parameters, ykeys, y0, int(opts.ltime), opts.isteps, opts.verbose)
    ofile = opts.ofile
    if len(ofile) < 1:
        ofile = '%s/fig' % getOfileName(opts.diatom, opts.pseudoaltero, opts.altero, opts.flavo, opts.pseudo)

    global obj, dbm
    obj_fm, obj_mm = getObjective(opts.pseudoaltero, opts.flavo, opts.altero, opts.pseudo)
    dbm_fm, dbm_mm = getObjective(False, False, False, False)

    bk = [opts.pseudoaltero, opts.flavo, opts.altero, opts.pseudo]
    global bacteriakeys
    bacteriakeys = filter(lambda x: bk[orgkeys.index(x) - 1], orgkeys[1:])

    psets = []
    for j in jsondicts:
        if len(j.split('/')) < 2:
            psets.append(j.split('.')[-2])
        else:
            psets.append(j.split('/')[-2])

    if opts.fullmedia:
        obj = obj_fm
        dbm = dbm_fm
    else:
        obj = obj_mm
        dbm = dbm_mm

    if opts.paperfig:
        plotFigure5(odesys, opts.fullmedia, int(opts.ltime), ('%s-%s' % (ofile, psets[0])), opts.otitle)
        plotFigureMetabolites(odesys, opts.fullmedia, int(opts.ltime), ('%s-%s' % (ofile, psets[0])), opts.otitle)
        return

        
    plotBiomassAndAbundances(odesys, opts.fullmedia, int(opts.ltime), ('%s-%s' % (ofile, psets[0])), jsondicts[0])
    plotMetabolites(odesys, opts.fullmedia, int(opts.ltime), ('%s-%s' % (ofile, psets[0])), jsondicts[0])

    for i, jd in enumerate(jsondicts[1:]):
        if not os.path.isfile(jd):
            continue
        with open(jd) as infile:
            paramdicts = json.load(infile)
        odesys = run_coop(functions, mergeParameters(paramdicts), ykeys, y0, int(opts.ltime), opts.isteps, opts.verbose)
        plotBiomassAndAbundances(odesys, opts.fullmedia, int(opts.ltime), ('%s-%s' % (ofile, psets[i+1])), jd)
        plotMetabolites(odesys, opts.fullmedia, int(opts.ltime), ('%s-%s' % (ofile, psets[i+1])), jd)


def getOfileName(diatom, pseudoaltero, altero, flavo, pseudo):
    '''
    If the output path is not specified, creates one according to the organisms in the ODE system
    '''
    dirname = './output'
    if diatom:
        dirname = dirname+'-Diatom'
    if pseudoaltero:
        dirname = dirname+'-Pseudoaltero'
    if flavo:
        dirname = dirname+'-Flavo'
    if altero:
        dirname = dirname+'-Altero'
    if pseudo:
        dirname = dirname+'-Pseudo'
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    return dirname

def run_coop(functions, parameters, ykeys, y0, runltime, runisteps, verbose=False):
    '''
    Solve the ODE system initializing the SimpleODE object and calling the integrateODE method
    '''
    odesys = sode.SimpleODE(functions, parameters, ykeys, y0)
    odesys.integrateODE(maxT=runltime, verbose=verbose, stepChoiceLevel=(0, 10, runisteps))

    return odesys

def mergeParameters(pdict):
    '''
    Helper to translate dictionary keys
    '''
    parameters = {}
    for k1 in pdict.keys():
        for k2 in pdict[k1].keys():
            parameters[k1+'_'+k2] = pdict[k1][k2]
    if not parameters.get('k_diatom_K_fe3_b', False):
        parameters['k_diatom_K_fe3_b'] = parameters['k_diatom_K_fe3']
    if not parameters.get('k_diatom_K_vit_b', False):
        parameters['k_diatom_K_vit_b'] = parameters['k_diatom_K_vit']
    return parameters

def buildODESys(jsondicts, runfullmedia, runminimalmedia, rundiatom, runpseudoaltero, runaltero, runflavo, runpseudo, verbose=False):
    '''
    Build the system of Ordinary Differential Equation according to the options passed (media type, organisms etc.)
    Loads the parameters from a json file or from 
    '''
    if '.json' in jsondicts:
        with open(jsondicts) as infile:
            paramdicts = json.load(infile)
    else:
        paramdicts = jsondicts

    parameters = mergeParameters(paramdicts)

    ykeys = ['DIATOM', 'PSEUDOALTERO', 'FLAVO', 'ALTEROMO', 'PSEUDOMO', 'VIT', 'FE', 'DOCA', 'DOCPA', 'COP', 'EPA', 'BAC', 'DOMA']
    y0 = [0.001, 0.9, 0.06, 0.02, 0.02, 1., 1.,  0.1, 0.1, 0.1, 0.0, 0.0, 0.0]
    if runminimalmedia:
        y0[ykeys.index('VIT')] = 0.1
        y0[ykeys.index('FE')] = 0.1
        y0[ykeys.index('DIATOM')] = 0.003
    elif runfullmedia:
        y0[ykeys.index('VIT')] = 1.
        y0[ykeys.index('FE')] = 1.
        y0[ykeys.index('DIATOM')] = 0.04
    #y0[ykeys.index('FE')] = 10.
    #y0[ykeys.index('VIT')] = 10.
    gammad = '(self.p["k_diatom_growth"]*(VIT/(VIT+self.p["k_diatom_K_vit"]))*(FE/(FE+self.p["k_diatom_K_fe3"]))*(1 - DIATOM/self.p["k_diatom_cc"]))'
    deltad = 'self.p["k_diatom_death"]*1/(1 + '+gammad+')'
    fd = gammad+'*DIATOM - '+deltad+'*DIATOM'
    fpa = '0.'
    ffl = '0.'
    fa = '0.'
    fp = '0.'

    fdocacop = '((VIT**4)/((VIT**4)+self.p["k_diatom_K_vit_b"]))*((FE**4)/((FE**4)+self.p["k_diatom_K_fe3_b"]))'
    fdocpa = 'self.p["k_diatom_docpa_prod"]*'+gammad+'*DIATOM'
    fdoca = '(self.p["k_diatom_doca_prod"] + self.p["k_diatom_V_doca-cop"]*( 1 - '+fdocacop+'))*DIATOM*'+gammad
    fc = '(self.p["k_diatom_cop_prod"] + self.p["k_diatom_V_doca-cop"]*'+fdocacop+')*DIATOM*'+gammad
    fv = ' - self.p["k_diatom_vit_cons"]*'+gammad+'*DIATOM'
    ff = ' - self.p["k_diatom_fe3_cons"]*'+gammad+'*DIATOM'
    fe = 'self.p["k_diatom_epa_prod"]*DIATOM*'+gammad+' - self.p["k_degradation_epa"]*EPA'
    fb = ' - self.p["k_degradation_bac"]*BAC'
    fdoma = 'self.p["k_diatom_dom_prod"]*'+deltad+'*DIATOM'

    if runpseudoaltero:
        gammapa = '(self.p["k_pseudoaltero_growth"]*(DOCPA/(DOCPA+self.p["k_pseudoaltero_K_docpa"]))*(1 - PSEUDOALTERO/self.p["k_pseudoaltero_cc"]))'
        fpa = gammapa+'*PSEUDOALTERO - self.p["k_pseudoaltero_death"]*PSEUDOALTERO*(1 + self.p["k_pseudoaltero_Vmax_epa"]*EPA/(EPA+self.p["k_pseudoaltero_K_epa"]))*1/(1 + '+gammapa+')'
        fdocpa = fdocpa+' - self.p["k_pseudoaltero_docpa_cons"]*PSEUDOALTERO*'+gammapa
        fb = fb+' + self.p["k_pseudoaltero_bac_prod"]*PSEUDOALTERO*'+gammapa
    else:
        y0[ykeys.index('PSEUDOALTERO')] = 0.

    if runflavo:
        gammaf = '(self.p["k_flavo_growth"]*(COP/(COP+self.p["k_flavo_K_cop"]))*(1 - FLAVO/self.p["k_flavo_cc"]))'
        ffl = gammaf+'*FLAVO - self.p["k_flavo_death"]*FLAVO*(1 + self.p["k_flavo_Vmax_bac"]*BAC/(BAC+self.p["k_flavo_K_bac"]))*1/(1 + '+gammaf+')'
        fc = fc+' - self.p["k_flavo_cop_cons"]*FLAVO*'+gammaf
    else:
        y0[ykeys.index('FLAVO')] = 0.

    if runaltero:
        gammaa_doca = '(self.p["k_altero_doca_frac"]*DOCA/(DOCA+self.p["k_altero_K_doca"]))'
        gammaa_dom =  '(self.p["k_altero_dom_frac"]*DOMA/(DOMA+self.p["k_altero_K_dom"]))'
        gammaa_com = '(self.p["k_altero_growth"]*(1 - ALTEROMO/self.p["k_altero_cc"]))'
        gammaa = '('+gammaa_dom+' + '+gammaa_doca+')*'+gammaa_com
        gammaa_doca = gammaa_doca+'*'+gammaa_com
        gammaa_dom = gammaa_dom+'*'+gammaa_com
        fa = gammaa+'*ALTEROMO - self.p["k_altero_death"]*ALTEROMO*(1 + self.p["k_altero_Vmax_bac"]*BAC/(BAC+self.p["k_altero_K_bac"]))*1/(1 + '+gammaa+')'
        fdoca = fdoca+' - self.p["k_altero_doca_cons"]*ALTEROMO*'+gammaa_doca
        fdoma = fdoma+' - self.p["k_altero_dom_cons"]*ALTEROMO*'+gammaa_dom
        fv = fv+' + self.p["k_altero_vit_prod"]*ALTEROMO*'+gammaa_doca
        ff = ff+' + self.p["k_altero_fe3_prod"]*ALTEROMO*'+gammaa_doca
    else:
        y0[ykeys.index('ALTEROMO')] = 0.

    if runpseudo:
        gammap = '(self.p["k_pseudo_growth"]*(DOMA/(DOMA+self.p["k_pseudo_K_dom"]))*(1 - PSEUDOMO/self.p["k_pseudo_cc"]))'
        fp = gammap+'*PSEUDOMO - self.p["k_pseudo_death"]*PSEUDOMO*(1 + self.p["k_pseudo_Vmax_bac"]*BAC/(BAC+self.p["k_pseudo_K_bac"]))*1/(1 + '+gammap+')'
        fdoma = fdoma+' - self.p["k_pseudo_dom_cons"]*PSEUDOMO*'+gammap
    else:
        y0[ykeys.index('PSEUDOMO')] = 0.

    if not rundiatom:
        y0[ykeys.index('DIATOM')] = 0.
                

    if verbose:
        functions = {'DIATOM':fd, 'VIT':fv, 'FE':ff, 'PSEUDOALTERO': fpa, 'FLAVO': ffl, 'ALTEROMO': fa, 'PSEUDOMO': fp, 'DOCA':fdoca, 'DOCPA':fdocpa, 'COP':fc, 'EPA':fe, 'BAC':fb, 'DOMA': fdoma}
        for k,v in functions.items():
            print('d[%s]/dt = %s' % (k, v))
        
    for k in ykeys:
        fd = fd.replace(k, writey(k, ykeys))
        fa = fa.replace(k, writey(k, ykeys))
        fp = fp.replace(k, writey(k, ykeys))
        fpa = fpa.replace(k, writey(k, ykeys))
        ffl = ffl.replace(k, writey(k, ykeys))
        fv = fv.replace(k, writey(k, ykeys))
        ff = ff.replace(k, writey(k, ykeys))
        fdoca = fdoca.replace(k, writey(k, ykeys))
        fdocpa = fdocpa.replace(k, writey(k, ykeys))
        fdoma = fdoma.replace(k, writey(k, ykeys))
        fe = fe.replace(k, writey(k, ykeys))
        fc = fc.replace(k, writey(k, ykeys))
        fb = fb.replace(k, writey(k, ykeys))

    functions = {'DIATOM':fd, 'VIT':fv, 'FE':ff, 'PSEUDOALTERO': fpa, 'FLAVO': ffl, 'ALTEROMO': fa, 'PSEUDOMO': fp, 'DOCA':fdoca, 'DOCPA':fdocpa, 'COP':fc, 'EPA':fe, 'BAC':fb, 'DOMA': fdoma}

    return functions, parameters, ykeys, y0


def plotFigure5(odesys, fm, ltime, oname, otitle):
    '''
    Plot Biomasses over time (top panel) and bacteria relative abundances over time (bottom panel)
    '''
    fig5, (ax5a, ax5b) = plt.subplots(2, sharex=True)
    figtitle = 'Community dynamics in '
    if fm:
        figtitle = figtitle+'Complete Media'
    else:
        figtitle = figtitle+'Minimal Media'
    figtitle = figtitle+otitle
    ax5a.set_title(figtitle)
    ax5aa = ax5b.twinx()
    ax5b.set_xlabel('Time [a.u.]')
    ax5a.set_ylabel('Biomass [a.u.]')
    ax5aa.set_ylabel('Tot Bacteria Biomass [a.u.]')
    ax5b.set_ylabel('Relative abundances')
    ax5b.set_xlim([0, ltime])
    timex = np.array(odesys.T)

    totbacbm = np.zeros(len(odesys.Y['PSEUDOALTERO']))
    for nk in orgkeys:
        ax5a.plot(timex, odesys.Y[nk], '-', color=orgcolors[nk], label=orgkeysshort[nk])
        if 'DIATOM' not in nk:
            totbacbm += odesys.Y[nk]
    objarr = {'time': np.array(dbm.keys(), dtype=float),
              'data': np.array(dbm.values(), dtype=float)}
    ax5aa.plot(timex, totbacbm, ':k', label='Tot Bacteria')
    ax5a.plot(objarr['time'], objarr['data'], 's', color=orgcolors['DIATOM'])
    ax5a.plot(-10, 0, 's', color='k', label='Data')
    l5 = ax5a.legend(loc='best', prop={'size':10}, scatterpoints=1)
    for i, k in enumerate(bacteriakeys):
        ax5b.plot(timex, odesys.Y[k]/totbacbm, '-', color=orgcolors[k], label=k.lower()+' : tot bacteria')
        for t, v in obj.items():
            ax5b.plot(t, v[i], 's', color=orgcolors[k])
    for f1 in [ax5a, ax5b]:
        f1.axvspan(0, 40, alpha=0.08, color='green')
        f1.axvspan(40, 72, alpha=0.08, color='blue')
        endofstat = 152
        mstr = 'mm'
        if fm:
            endofstat = 128
            mstr = 'fm'
        f1.axvspan(72, endofstat, alpha=0.06, color='red')
        f1.axvspan(endofstat, ltime, alpha=0.1, color='yellow')
    l5 = ax5aa.legend(loc='best', prop={'size':10}, scatterpoints=1)
    fig5.savefig('%s-fig5-%s.png'%(oname,mstr))
    return

def plotFigureMetabolites(odesys, fm, ltime, oname, otitle):
    '''
    Plot Metabolites over Time
    '''
    fig5, (ax5a, ax5b, ax5c) = plt.subplots(3, sharex=True)
    figtitle = 'Metabolite dynamics in '
    if fm:
        figtitle = figtitle+'Complete Media'
    else:
        figtitle = figtitle+'Minimal Media'
    figtitle = figtitle+otitle
    ax5a.set_title(figtitle)
    ax5c.set_xlabel('Time [a.u.]')
    ax5c.set_xlim([0, ltime])
    ax5a.set_ylabel('Micronutrients [a.u.]')
    ax5b.set_ylabel('Bactericides [a.u.]')
    ax5c.set_ylabel('Nutrients [a.u.]')
    timex = np.array(odesys.T)

    c = 0
    for k in ['VIT', 'FE']:
        ax5a.plot(timex, odesys.Y[k], '-', color=colors_rgb[c], label=k.lower())
        c +=1
    for k in ['EPA', 'BAC']:
        ax5b.plot(timex, odesys.Y[k], '-', color=colors_rgb[c], label=k.lower())
        c +=1
    for k in ['DOCA', 'DOCPA', 'COP', 'DOMA']:
        ax5c.plot(timex, odesys.Y[k], '-', color=colors_rgb[c], label=k.lower().replace('doma', 'dom'))
        c +=1

    l5a = ax5a.legend(loc='best', prop={'size':10}, scatterpoints=1)
    l5b = ax5b.legend(loc='best', prop={'size':10}, scatterpoints=1)
    l5c = ax5c.legend(loc='best', prop={'size':10}, scatterpoints=1)
    for f1 in [ax5a, ax5b, ax5c]:
        f1.axvspan(0, 40, alpha=0.08, color='green')
        f1.axvspan(40, 72, alpha=0.08, color='blue')
        endofstat = 152
        mstr = 'mm'
        if fm:
            endofstat = 128
            mstr = 'fm'
        f1.axvspan(72, endofstat, alpha=0.06, color='red')
        f1.axvspan(endofstat, ltime, alpha=0.1, color='yellow')
    fig5.savefig('%s-metabolites-%s.png'%(oname,mstr))
    return

def plotBiomassAndAbundances(odesys, fm, ltime, oname, otitle):
    '''
    Plot Biomasses over time (top panel) and bacteria relative abundances over time (bottom panel)
    '''
    fig5, (ax5a, ax5b) = plt.subplots(2, sharex=True)
    ax5a.set_title(otitle)
    ax5b.set_xlabel('Time [a.u.]')
    ax5a.set_ylabel('Biomass [a.u.]')
    ax5b.set_ylabel('Relative abundances')
    ax5b.set_xlim([0, ltime])
    timex = np.array(odesys.T)

    totbacbm = np.zeros(len(odesys.Y['PSEUDOALTERO']))
    for nk in orgkeys:
        ax5a.plot(timex, odesys.Y[nk], '-', color=orgcolors[nk], label=nk.lower())
        if 'DIATOM' not in nk:
            totbacbm += odesys.Y[nk]
    objarr = {'time': np.array(dbm.keys(), dtype=float),
              'data': np.array(dbm.values(), dtype=float)}
    ax5a.plot(objarr['time'], objarr['data'], 's', color=orgcolors['DIATOM'], label='Data')
    l5 = ax5a.legend(loc='best', prop={'size':10}, scatterpoints=1)
    for i, k in enumerate(bacteriakeys):
        ax5b.plot(timex, odesys.Y[k]/totbacbm, '-', color=orgcolors[k], label=k.lower()+' : tot bacteria')
        for t, v in obj.items():
            ax5b.plot(t, v[i], 's', color=orgcolors[k])
    for f1 in [ax5a, ax5b]:
        f1.axvspan(0, 40, alpha=0.08, color='green')
        f1.axvspan(40, 72, alpha=0.08, color='blue')
        endofstat = 152
        mstr = 'mm'
        if fm:
            endofstat = 128
            mstr = 'fm'
        f1.axvspan(72, endofstat, alpha=0.06, color='red')
        f1.axvspan(endofstat, ltime, alpha=0.1, color='yellow')
    fig5.savefig('%s-biomasses-abundances-%s.png'%(oname,mstr))
    return

def plotMetabolites(odesys, fm, ltime, oname, otitle):
    '''
    Plot Metabolites over Time
    '''
    fig5, (ax5a, ax5b, ax5c) = plt.subplots(3, sharex=True)
    ax5a.set_title(otitle)
    ax5c.set_xlabel('Time [a.u.]')
    ax5c.set_xlim([0, ltime])
    ax5a.set_ylabel('Micronutrients [a.u.]')
    ax5b.set_ylabel('Bactericides [a.u.]')
    ax5c.set_ylabel('Nutrients [a.u.]')
    timex = np.array(odesys.T)

    c = 0
    for k in ['VIT', 'FE']:
        ax5a.plot(timex, odesys.Y[k], '-', color=colors_rgb[c], label=k.lower())
        c +=1
    for k in ['EPA', 'BAC']:
        ax5b.plot(timex, odesys.Y[k], '-', color=colors_rgb[c], label=k.lower())
        c +=1
    for k in ['DOCA', 'DOCPA', 'COP', 'DOMA']:
        ax5c.plot(timex, odesys.Y[k], '-', color=colors_rgb[c], label=k.lower().replace('doma', 'dom'))
        c +=1

    l5a = ax5a.legend(loc='best', prop={'size':10}, scatterpoints=1)
    l5b = ax5b.legend(loc='best', prop={'size':10}, scatterpoints=1)
    l5c = ax5c.legend(loc='best', prop={'size':10}, scatterpoints=1)
    for f1 in [ax5a, ax5b, ax5c]:
        f1.axvspan(0, 40, alpha=0.08, color='green')
        f1.axvspan(40, 72, alpha=0.08, color='blue')
        endofstat = 152
        mstr = 'mm'
        if fm:
            endofstat = 128
            mstr = 'fm'
        f1.axvspan(72, endofstat, alpha=0.06, color='red')
        f1.axvspan(endofstat, ltime, alpha=0.1, color='yellow')
    fig5.savefig('%s-metabolites-%s.png'%(oname,mstr))
    return

def writey(keyword, keylist):
    '''
    Helper function for translation of equation
    '''
    return 'y[%d]'%keylist.index(keyword)

def options():
    '''define here in-line arguments'''
    parser = optparse.OptionParser(description='Parsing options')
    parser.add_option('-v', '--verbose', dest='verbose', help='increase output verbosity', action='store_true')
    parser.add_option('-f', '--fullmedia', dest='fullmedia', help='run with full media conditions', action='store_true')
    parser.add_option('-m', '--minimalmedia', dest='minimalmedia', help='run with minimal media conditions', action='store_true')
    parser.add_option('-a', '--altero', dest='altero', help='add alteromonadacea organism to ODE model', action='store_true')
    parser.add_option('-b', '--pseudoaltero', dest='pseudoaltero', help='add pseudoalteromonadacea organism to ODE model', action='store_true')
    parser.add_option('-c', '--pseudo', dest='pseudo', help='add pseudomonadacea organism to ODE model', action='store_true')
    parser.add_option('-d', '--diatom', dest='diatom', help='add diatom organism to ODE model', action='store_true')
    parser.add_option('-e', '--flavo', dest='flavo', help='add flavobacteriadacea to ODE model', action='store_true')
    #parser.add_option('-x', '--fluxplots', dest='fluxplots', help='make flux plots', action='store_true')
    #parser.add_option('-s', '--scan', dest='scan', help='scan over some parameters', default='[]')
    parser.add_option('-o', '--ofile', dest='ofile', help='output file name', default='')
    parser.add_option('-t', '--otitle', dest='otitle', help='output plot title', default='')
    parser.add_option('-i', '--isteps', dest='isteps', help='set maximum number of integration steps', default='3000')
    parser.add_option('-l', '--ltime', dest='ltime', help='set maximum time for solver', default='300')
    parser.add_option('-j', '--jsondicts', dest='jsondicts', help='space separated list of json file of parameter dictionaries', default='parametersMM.json parametersCM.json')
    parser.add_option('-P', '--paperfig', dest='paperfig', help='producing paper fig', action='store_true')
    parser.add_option('-g', '--debug', dest='debug', help='debug option', action='store_true')
    opts, args = parser.parse_args()
    if opts.verbose:
        print "verbosity turned on"
        print opts
        print args
    return opts, args

if __name__=="__main__":
    main()
