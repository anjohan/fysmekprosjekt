import numpy as np
import matplotlib.pylab as plt
from NAtoms import *
import time

"""Program that solves the tasks given in the 'Large Systems' section of the project.
The names of the functions below aim to solve the tasks belonging to the subsection
with the same/similar name."""

def plot_LJ_cutoff():
    """ Plots the potential and force described in the 'Implementation' subsection."""

    plot_LJ()

def verification():

    L = 1   #The value doesn't matter since we won't use boundary conditions
    T = 10
    dt = 0.01

    r0 = np.array([0, 0, 0, 1.5, 0, 0], dtype=float)
    v0 = np.zeros(len(r0), dtype=float)

    r, v, t = simulate(r0, v0, L, T, dt, 'VelVerlet', boundary = False, wrt_file = 'ver_2.xyz')

    for r_pert, word in zip([0, 0.1], ['without', 'with']):

        r0 = np.array([1, r_pert, 0, 0, 1, 0, -1, 0, 0, 0, -1, 0], dtype=float)
        v0 = np.zeros(len(r0), dtype=float)

        r, v, t = simulate(r0, v0, L, T, dt, 'VelVerlet', boundary = False, wrt_file = 'ver_'+word+'.xyz')

        KinEng, PotEng = Energy(r, v, t)

        plt.plot(t, (KinEng + PotEng))
        plt.title("Total energy as function of time %s perturbation" % (word))
        plt.xlabel(r'$t$')
        plt.ylabel(r'$E(t)$')
        plt.grid(True)
        plt.show()

def initialization():

    #N = 108 # 4, 32, 108, 256, 500, 864
    nBoxes = 3
    L = 20
    m = 1

    T = 1
    dt = 1

    r0 = fcc(nBoxes,L)
    v0 = np.zeros(len(r0))

    r, v, t = simulate(r0, v0, L, T, dt, 'VelVerlet', boundary=False, wrt_file = 'fcc.xyz')

def many_atoms():

    #N = 50
    L = 10
    T = 10
    dt = 0.01
    nt = int(T/dt) + 1

    r0 = fcc(3, L)

    v0 = np.zeros(len(r0))

    r, v, t = simulate(r0, v0, L, T, dt, 'VelVerlet', boundary=False, wrt_file= 'open.xyz')

    KinEng, PotEng = Energy(r, v, t)

    plt.plot(t, (KinEng + PotEng))
    plt.title("Total energy as function of time")
    plt.xlabel(r'$t$')
    plt.ylabel(r'$E(t)$')
    plt.grid(True)
    plt.show()



def boundary_conditions():

    L = 6.8
    T = 3
    dt = 0.01
    nt = int(T/dt) + 1

    r0 = fcc(4, L)

    v0 = np.random.uniform(-2, 2, size=(len(r0)))

    r, v, t = simulate(r0, v0, L, T, dt, 'VelVerlet', wrt_file= 'bc.xyz')
    return r,v,t


def test_vel():

    nBoxes = 5
    dens = 0.807
    N = 4*nBoxes**3

    L = (N/dens)**(1/3)
    T = 1.5
    dt = 0.005
    nt = int(T/dt) + 1

    n = 20   # no. experiments

    for i in range(n):
        r0 = fcc(nBoxes, L)

        #v0 = np.random.uniform(-1.5, 1.5, size=(len(r0)))
        v0= np.random.normal(0, 0.885, size=(len(r0)))

        r, v, t = simulate(r0, v0, L, T, dt, 'VelVerlet')

        v_corr = velocity_corr(v, t, plot=False)

        if i ==0:
            v_ave = v_corr
        else:
            v_ave += v_corr

    v_ave /= n
    t_corr = t[:len(v_ave)]

    sigma = 3.405*10**(-8)  # cm
    tau = 2.1569*10**(-12)  # s

    intgr = np.trapz(v_ave, t_corr, dx = 0.1)/3
    intgr *= sigma**2/tau
    print("Diffusion constant: ", intgr, "cm^2 sec^-1")

    plt.plot(t_corr*tau*10**(12), v_ave)
    plt.xlabel(r'$t [10^{12} s]$')
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

    r, v, t = simulate(r0, v0, L, T, dt, 'VelVerlet', wrt_file = False)

    rdf(r, L)


if __name__ == "__main__":
    #plot_LJ_cutoff()
    #verification()
    #initialization()
    #many_atoms()
    boundary_conditions()
    #test_vel()
    #test_rdf()
