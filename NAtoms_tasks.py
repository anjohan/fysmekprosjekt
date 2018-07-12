import numpy as np
import matplotlib.pylab as plt
from NAtoms import *
import time

def test():

    N = 100
    T = 5
    dt = 0.05
    nt = int(T/dt) + 1
    m = 1

    #r0 = np.random.randint(1, 5, size=3*N)
    r0 = np.random.uniform(0.0, 50.0, size = 3*N)
    #r0 = np.array([-2, 0, 0, 2, 0, 0, 0, -2, 0, 0, 2, 0])

    v0 = np.zeros(3*N)

    t_start = time.clock()

    r, v, t = simulate(r0, v0, m, T, dt, 'VelVerlet', wrt_file=True)

    t_stop = time.clock()

    print("Tot. time: %.3f s" % (t_stop - t_start))
    print("Time pr iter: %.3f s" % ((t_stop - t_start)/(nt-1)))

    KinEng, PotEng = Energy(r, v, t, m)

    plt.plot(t, (KinEng + PotEng))
    plt.title("Total energy as function of time")
    plt.xlabel(r'$t$')
    plt.ylabel(r'$E(t)$')
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    test()
