"""Brownian motion example
This module contains functions for simulating Brownian motion
and analyzing the results (lectures 4,5 and lab 2, task 2)
"""
import numpy as np
import matplotlib.pyplot as plt
from time import time

def brown0(Nt,M,dt=1):
    """ Inefficient version of brown1
    """
    X = np.zeros((Nt+1,M)) #initialize array

    #Run M Nt-step simulations
    for j in range(M):
        for i in range(Nt):
            X[i+1,j] = X[i,j] + np.sqrt(dt)*np.random.randn(1)

    #Compute statistics
    Xm = np.mean(X,axis=1)
    Xv = np.var(X,axis=1)

    return X,Xm,Xv



def brown1(Nt,M,dt=1):
    """Run M Brownian motion simulations each consisting of Nt time steps
    with time step = dt
    Returns: X: the M trajectories; Xm: the mean across these M samples; Xv:
    the variance across these M samples
    """
    from numpy.random import randn

    N = np.sqrt(dt)*randn(Nt+1,M) #Generate all needed random numbers
    N[0,:] = 0. #Initial conditions, left out during lecture 5

    X = np.cumsum(N,axis=0) #Compute

    Xm = np.mean(X,axis=1)
    Xv = np.var(X,axis=1)
    return X,Xm,Xv


def analyze(display=False):
    """Lab2 task 2: Compute variance in M Brownian motion simulations with varying M
    Compute and return error along with M values and variances.
    Plot error if input variable display=True
    To run this function, first import brown in python terminal, then: brown.analyze(True)
    """

    Mvalues = [10,100,1000,10000] #Compute error for these values of M
    Nt=100
    Xvarray = np.zeros(len(Mvalues)) #initialize array to store variances at t=Nt+1

    #Compute variances for each M
    for i,M in enumerate(Mvalues):
        _,_,Xv = brown1(Nt,M)
        print(i,M,Xv[-1])
        Xvarray[i] = Xv[-1]

    errorv = np.abs(Xvarray-Nt) #we expect Xv(t) = t for dt = 1

    if display:
        plt.figure()
        plt.loglog(Mvalues,errorv,'x--')

        #construct least-squares fit to: error  = A M^n
        #we expect n=-1, but in practice, the observed value
        #can vary substantially
        p=np.polyfit(np.log(Mvalues),np.log(errorv),1)
        n = p[0]; A = p[1]
        plt.plot(Mvalues,np.exp(A)*(Mvalues)**n,'r--')
        plt.legend(('simulation','least-squares fit with n=%f'%(n)),loc='best')
        plt.xlabel('M')
        plt.ylabel('$\epsilon$')
        plt.title('Variation of variance with sample size')
    return Mvalues,Xvarray,errorv


if __name__ == '__main__':
    #Compare speed of brown0 and brown1, use run command in terminal to
    #execute this block
    t1 = time()
    _,_,_ = brown0(1,500,200)
    t2 = time()
    dt0 = t2-t1

    t3 = time()
    _,_,_ = brown1(500,200)
    t4 = time()
    dt1 = t4-t3

    print("dt0=",dt0)
    print("dt1=",dt1)
