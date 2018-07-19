import numpy as np
import matplotlib.pylab as plt
from TwoAtoms import *

"""Program that solves the tasks given in the Introduction and Two-Atom sections of the project.
The names of the functions below aim to solve the tasks belonging to the subsection
with the same/similar name."""

def potential_force():
    """Solves the 'Understanding the potential' and 'Forces and equations of motion' sections."""

    plot_LJ()


def motion_energy():
    """Solves the 'Motion' and first part of the 'Energy' sections."""

    r0 = np.array([0, 0, 0, 0, 0, 0], dtype=float)
    v0 = np.array([0, 0, 0, 0, 0, 0], dtype=float)

    for rx in [1.5, 0.95]:

        r0[3] = rx

        r, v, t = simulate(r0, v0, 5, 0.01, 'EulerCromer', wrt_file='2atom_'+str(rx)+'.xyz')

        d_abs = np.sqrt((np.sum((r[:, :3] - r[:, 3:])**2, axis=1)))

        plt.plot(t, d_abs)
        plt.title("Distance as function of time with initial distance r = %g" % (rx))
        plt.xlabel(r"$t$")
        plt.ylabel(r"$r$")
        plt.grid(True)
        plt.show()

        KinEng, PotEng = Energy(r, v, t, plot=True)

        plt.plot(t, (KinEng + PotEng))
        plt.title("Total energy as function of time for initial distance r = %g" % (rx))
        plt.xlabel(r'$t$')
        plt.ylabel(r'$E(t)$')
        plt.grid(True)
        plt.show()



def compare_methods():
    """Solves the last part of the 'Energy' section, which tasks to compare the integration methods"""

    r0 = np.array([0, 0, 0, 1.5, 0, 0], dtype=float)
    v0 = np.array([0, 0, 0, 0, 0, 0], dtype=float)

    methods = ['Euler', 'EulerCromer', 'VelVerlet']
    for method in methods:

        r, v, t = simulate(r0, v0, 5, 0.01, method)

        KinEng, PotEng = Energy(r, v, t, plot=False)

        plt.plot(t, (KinEng + PotEng), label=method)

    plt.title("Energy conservation comparison, dt = 0.01")
    plt.xlabel(r'$t$')
    plt.ylabel(r'$E(t)$')
    plt.grid(True)
    plt.legend()
    plt.show()

    #Values found through trial and error.
    for dt, method in zip([0.005, 0.075, 0.075], ['Euler', 'EulerCromer', 'VelVerlet']):

        r, v, t = simulate(r0, v0, 5, dt, method)
        KinEng, PotEng = Energy(r, v, t, plot=False)

        r_abs = np.sqrt((np.sum((r[:, :3] - r[:, 3:])**2, axis=1)))

        fig, (ax1, ax2) = plt.subplots(2, 1)
        ax1.plot(t, r_abs)
        ax1.set_title("Radial distance and energy with %s, dt=%g" %(method, dt))
        ax1.set_xlabel(r"$t$")
        ax1.set_ylabel(r"$r(t)$")
        ax1.grid(True)
        ax2.plot(t, (KinEng + PotEng))
        ax2.set_xlabel(r"$t$")
        ax2.set_ylabel(r"$E(t)$")
        ax2.grid(True)
        plt.show()


if __name__ == "__main__":
    potential_force()
    motion_energy()
    compare_methods()
