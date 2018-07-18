import numpy as np
import matplotlib.pylab as plt
from NAtoms import *
import time

def test():

    #N = 50
    L = 15
    T = 10
    dt = 0.01
    nt = int(T/dt) + 1
    m = 1

    #r0 = np.array([-2, 0, 0, 2, 0, 0, 0, -2, 0, 0, 2, 0])
    #r0 = np.random.uniform(0.0, L, size=(3*N))

    r0 = fcc(4, L)

    v0 = np.random.uniform(-2, 2, size=(len(r0)))
    #v0 = np.zeros(len(r0))
    #v0 = np.array([-1, 0, 0, 1, 0, 0])

    r, v, t = simulate(r0, v0, L, m, T, dt, 'VelVerlet', wrt_file = True)


    KinEng, PotEng = Energy(r, v, t, m)

    plt.plot(t, (KinEng + PotEng))
    plt.title("Total energy as function of time")
    plt.xlabel(r'$t$')
    plt.ylabel(r'$E(t)$')
    plt.grid(True)
    plt.show()

def test_vel():

    L = 10
    T = 10
    dt = 0.01
    nt = int(T/dt) + 1
    m = 1

    n = 20   # no. experiments
    v_ave = np.zeros(nt)

    for i in range(n):
        r0 = fcc(3, L)

        v0 = np.random.uniform(-2, 2, size=(len(r0)))

        r, v, t = simulate(r0, v0, L, m, T, dt, 'VelVerlet', wrt_file = False)

        v_corr = velocity_corr(v, t, plot=False)

        v_ave += v_corr

    v_ave /= n

    sigma = 3.405*10**(-8)
    tau = 2.1569*10**(-12)

    intgr = np.trapz(v_ave, t, dx = 0.1)/3
    intgr *= sigma**2/tau
    print("Diff: ", intgr)

    plt.plot(t, v_ave)
    plt.xlabel(r'$t$')
    plt.ylabel(r'$\langle v(0) \cdot v(t) \rangle$')
    plt.title("Velocity auto-correlation function")
    plt.grid(True)
    plt.show()

def test_rdf():

    #L = 15
    L = 8.525 # Corresponding to 500 atoms
    T = 10
    dt = 0.01
    nt = int(T/dt) + 1
    m = 1

    r0 = fcc(5, L)

    v0 = np.random.uniform(-2, 2, size=(len(r0)))

    r, v, t = simulate(r0, v0, L, m, T, dt, 'VelVerlet', wrt_file = False)

    rdf(r, L)


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
    #test_vel()
    test_rdf()
    #test()
