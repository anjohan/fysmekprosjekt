import numpy as np
import matplotlib.pylab as plt
from NAtoms import *
import time

def test():

    #N = 50
    L = 15
    T = 10
    dt = 0.02
    nt = int(T/dt) + 1
    m = 1

    #r0 = np.array([-2, 0, 0, 2, 0, 0, 0, -2, 0, 0, 2, 0])
    #r0 = np.random.uniform(0.0, L, size=(3*N))

    r0 = fcc(4, L)

    v0 = np.random.uniform(-2, 2, size=(len(r0)))
    #v0 = np.zeros(len(r0))
    #v0 = np.array([-1, 0, 0, 1, 0, 0])

    r, v, t = simulate(r0, v0, L, m, T, dt, 'VelVerlet', wrt_file = True)

    """
    KinEng, PotEng = Energy(r, v, t, m)

    plt.plot(t, (KinEng + PotEng))
    plt.title("Total energy as function of time")
    plt.xlabel(r'$t$')
    plt.ylabel(r'$E(t)$')
    plt.grid(True)
    plt.show()
    """

def test_fcc():

    #N = 108 # 4, 32, 108, 256, 500, 864
    nBoxes = 4
    L = 15
    m = 1

    T = 10
    dt = 0.05

    r0 = fcc(nBoxes,L)
    v0 = np.zeros(len(r0))

    t = simulate_sparse(r0, v0, m, T, dt, 'VelVerlet', 'fcc.xyz')

if __name__ == "__main__":
    #test_fcc()
    test()
