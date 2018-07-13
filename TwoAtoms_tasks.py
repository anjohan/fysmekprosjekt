import numpy as np
import matplotlib.pylab as plt
from TwoAtoms import *

def task_ab():

    plot_LJ()


def task_cd():

    m = 1

    r0 = np.array([0, 0, 0, 0, 0, 0], dtype=float)
    v0 = np.array([0, 0, 0, 0, 0, 0], dtype=float)

    for rx in [1.5, 0.95]:

        r0[3] = rx

        r, v, t = simulate(r0, v0, m, 10, 0.01, 'EulerCromer')

        d_abs = np.sqrt((np.sum((r[:, :3] - r[:, 3:])**2, axis=1)))

        plt.plot(t, d_abs)
        plt.title("Distance as function of time with initial distance r = %g" % (rx))
        plt.xlabel(r"$t$")
        plt.ylabel(r"$r$")
        plt.grid(True)
        plt.show()

        KinEng, PotEng = Energy(r, v, t, m, plot=False)

        plt.plot(t, (KinEng + PotEng))
        plt.title("Total energy as function of time for initial distance r = %g" % (rx))
        plt.xlabel(r'$t$')
        plt.ylabel(r'$E(t)$')
        plt.grid(True)
        plt.show()



def compare_methods():

    m = 1   # Mass set to 1 for simplicity

    r0 = np.array([0, 0, 0, 1.5, 0, 0], dtype=float)
    v0 = np.array([0, 0, 0, 0, 0, 0], dtype=float)

    methods = ['Euler', 'EulerCromer', 'VelVerlet']
    for method in methods:

        r, v, t = simulate(r0, v0, m, 10, 0.01, method)

        KinEng, PotEng = Energy(r, v, t, m, plot=False)

        plt.plot(t, (KinEng + PotEng), label=method)

    plt.xlabel(r'$t$')
    plt.ylabel(r'$E(t)$')
    plt.grid(True)
    plt.legend()
    plt.show()



def test():

    m = 1   # Mass set to 1 for simplicity

    #r0 = np.array([0, 0, 0, 2**(1/6), 1, 2**(1/6)])
    r0 = np.array([0, 0, 0, 1.5, 0, 0])
    v0 = np.array([0, 0, 0, 0, 0, 0])

    #d_abs = np.sqrt((np.sum((r[:, :3] - r[:, 3:])**2, axis=1)))

    #plt.plot(t, d_abs)
    #plt.show()

    r, v, t = simulate(r0, v0, m, 10, 0.01, 'VelVerlet', wrt_file=True)


if __name__ == "__main__":
    test()
    #task_ab()
    #task_cd()
    #compare_methods()
