import numpy as np
import matplotlib.pylab as plt
import sys

"""
This file contains all functions necessary to solve to the tasks for the 2-atom Lennard-Jones
problem. The tasks themselves will be solved in TwoAtoms_tasks.py, in which all the function
in this file will be imported.
"""

def LJ_pot(r, epsilon=1, sigma=1):
    # Function for computing the potential of the Lennard-Jones Potential

    return 4*epsilon*((sigma/r)**12 - (sigma/r)**6)

def LJ_force(r, epsilon=1, sigma=1):
    # Function for computing the forces of the Lennard-Jones Potential

    return 24*epsilon/r*(2*(sigma/r)**12 - (sigma/r)**6)


def plot_LJ():
    # Simple function for plotting the curves of the LJ potential and forces.
    # Simply for showing the shape of the curves, the r-array is a set of arbitrary r-values in an
    # interval near y = 0.

    r = np.linspace(0.2, 3, 1001)

    fig, (ax1, ax2) = plt.subplots(2, 1)
    ax1.set_title("Potential and Force curves of the LJ potential")
    ax1.plot(r, LJ_pot(r))
    ax1.set_ylim([-1, 3])
    ax1.set_xlabel(r"$r/\sigma$")
    ax1.set_ylabel(r"$V_{LJ}(r)$")
    ax1.grid(True)
    ax2.plot(r, LJ_force(r))
    ax2.set_ylim([-3, 3])
    ax2.set_xlabel(r"$r/\sigma$")
    ax2.set_ylabel(r"$F_{LJ}(r)$")
    ax2.grid(True)
    plt.show()



def intergrator(r, v, dt, m, intg):

    nt = v.shape[0]

    infile = open("in.xyz", 'w')

    # Time integration, Euler-Cromer temporarily
    for i in range(nt-1):

        d_vec = r[i, :3] - r[i, 3:]
        d = np.sqrt(np.sum((d_vec)**2))

        F = LJ_force(d)*d_vec/d

        # Forces between the two atoms are the same, but with opposite signs. Putting these in the same array.
        F = np.concatenate((F, -F))

        a = F/m

        if intg == 'Euler':
            v[i+1] = v[i] + a*dt
            r[i+1] = r[i] + v[i]*dt

        elif intg == 'EulerCromer':
            v[i+1] = v[i] + a*dt
            r[i+1] = r[i] + v[i+1]*dt

        elif intg == 'VelVerlet':
            r[i+1] = r[i] + v[i]*dt + 0.5*a*dt**2

            d_vec = r[i+1, :3] - r[i+1, 3:]
            d = np.sqrt(np.sum((d_vec)**2))

            F_next = LJ_force(d)*d_vec/d
            F_next = np.concatenate((F_next, -F_next))
            a_next = F_next/m

            v[i+1] = v[i] + 0.5*(a + a_next)*dt

        else:
            print('Error: The variable intg must either be "Euler", "EulerCromer" or "VelVerlet".\nExiting...')
            sys.exit(0)


        infile.write("2\n\n")
        for j in range(2):
            infile.write(str(r[i, 3*j:3*j +3])[1:-1] + "\n")

    infile.close()

    return r, v



def simulate(m, T, dt, intg):

    nt = int(T/dt) + 1   # Number of timesteps
    t = np.linspace(0, T, nt)

    r = np.zeros((nt, 6))
    r[0] = np.array([0, 0, 0, 2**(1/6), 1, 2**(1/6)])

    v = np.zeros((nt, 6))

    infile = open("in.xyz", 'w')

    r, v = intergrator(r, v, dt, m, intg)

    return r, v, t


def Energy(r, v, t, m, plot = True):

    nt = v.shape[0]
    KinEng = np.zeros(nt)
    PotEng = np.zeros(nt)

    for i in range(nt):
        for j in range(2):
            v_abs2 = np.sum(v[i, 3*j:3*j + 3]**2)
            KinEng[i] += 0.5*m*v_abs2

        d = np.sqrt(np.sum((r[i, :3] - r[i, 3:])**2))
        PotEng[i] = LJ_pot(d)

    if plot:
        fig, (ax1, ax2) = plt.subplots(2, 1)
        ax1.set_title("Kinetic and Potential energy of the 2-atom system")
        ax1.plot(t, KinEng)
        #ax1.set_ylim([-1, 3])
        ax1.set_xlabel(r"$t$")
        ax1.set_ylabel(r"$K(t)$")
        ax1.grid(True)
        ax2.plot(t, PotEng)
        #ax2.set_ylim([-3, 3])
        ax2.set_xlabel(r"$t$")
        ax2.set_ylabel(r"$U(t)$")
        ax2.grid(True)
        plt.show()

    return KinEng, PotEng



if __name__ == "__main__":
    #plot_LJ()
    #simulate(1, 10, 0.05)
