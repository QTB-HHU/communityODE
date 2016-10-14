#************************************
#**  author: Antonella Succurro    **
#**  email:a.succurro[AT]gmail.com **
#**                                **
#**  created: 2016/06/01           **
#**  last modified: 2016/08/18     **
#************************************

#from __future__ import print_function
import numpy as np
from scipy.integrate import ode, odeint
import warnings

'''
Begin helper code
'''
'''
End helper code
'''

class SimpleODE(object):
    '''
    Class to build and solve Ordinary Differential Equation systems
    '''
    def __init__(self, functions, parameters, ykeys, y0):
        '''
        functions = dictionary defining the equations
        parameters = dictionary of parameters
        '''
        self.f = functions
        self.p = parameters
        self.yk = ykeys
        self.dy = np.zeros(len(self.yk))
        self.y = [np.array(y0)]
        self.T = [0]
        self.Y = {}
        
    def odesys(self, y, verbosity=False):
        '''
        Interpret the ODE system equations from human-readable to python code
        Build the self.dy array for the solver
        '''
        for i, k in enumerate(self.yk):
            self.dy[i] = eval(self.f[k])
            if verbosity:
                print self.dy[i]
        return self.dy

    def integrateODE(self, integrator='vode', maxT=200, stepChoiceLevel=(0, 10, 1000), verbose=False, alwayspositive=True):
        '''
        Routine to integrate the ODE system
        '''
        def dydt(t, y, obj, verb=False):
            return obj.odesys(y, verb)

        solver = ode(dydt).set_integrator(integrator)
        if integrator == 'dopri5':
            solver.set_integrator(integrator, nsteps=1, max_step=stepChoiceLevel[1])
        elif integrator == 'vode':
            solver.set_integrator(integrator, max_step=stepChoiceLevel[1])
        else:
            print 'This option is available only for vode and dopri5 SciPy solvers!'
            return
        #set the parameters of the differential function dydt: model and verbose
        solver.set_f_params(self, verbose)
        solver._integrator.iwork[2] = -1
        warnings.filterwarnings("ignore", category=UserWarning)
        solver.set_initial_value(self.y[0], self.T[0])
        if verbose:
            print 'Initial Conditions to start integration: '
            print '\t t0: ', self.T[0]
            print '\t y0: ', self.y[0]

        step = 0
        while (solver.successful() and self.T[-1]<maxT and step < stepChoiceLevel[2]):
            step += 1
            solver.integrate(maxT, step=True)
            handleZeroEvent, metIdx = self.checkForNegativeSolutions(solver.y)
            if handleZeroEvent and alwayspositive:
                if verbose:
                    print '****************************'
                    print '!in handleZeroEvent!'
                slope = float(self.dy[metIdx])
                zeroDT = 0.-self.y[-1][metIdx]/slope
                self.y.append(np.zeros(len(self.yk)))
                for i,dy in enumerate(self.dy):
                    new_dy = float(dy)*zeroDT
                    if i == metIdx:
                        self.y[-1][i] = 0.
                    else:
                        self.y[-1][i] = self.y[-2][i]+new_dy
                self.T.append(self.T[-1]+zeroDT)
                solver.set_initial_value(self.y[-1], self.T[-1])
                if verbose:
                    print 'new initial conditions: '
                    print '\t t0: ', self.T[-1], ', solver t: ',solver.t
                    print '\t y0: ', self.y[-1]
            else:
                self.y.append(solver.y)
                if verbose:
                    print '    dt: ', solver.t-self.T[-1]
                    print '    dydt: ', self.dy
                    print '    dy: ', self.dy*(solver.t-self.T[-1])
                    print [self.y[-2:][m] for m in range(len(self.yk))]
                self.T.append(solver.t)
        warnings.resetwarnings()

        provarr = np.array(self.y)
        for i,k in enumerate(self.yk):
            self.Y[k] = provarr[:,i]
        return

    def checkForNegativeSolutions(self, ysol, yThr=0):
        '''
        Checks if some quantity went below zero
        Return boolean and the negative quantity index (None in case all positive)
        '''
        minY = min(ysol)
        isNeg = minY < yThr
        if isNeg:
            idNeg = np.where(ysol==minY)[0][0]
            return isNeg, idNeg
        else:
            return isNeg, None
