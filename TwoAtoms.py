import numpy as np
import matplotlib.pylab as plt
import sys

"""
This file contains all functions necessary to solve to the tasks for the 2-atom Lennard-Jones
problem. The tasks themselves will be solved in TwoAtoms_tasks.py, in which all the function
in this file will be imported.
"""

def LJ_pot(r, epsilon=1, sigma=1):
    """Function for computing the potential of the Lennard-Jones Potential.

    Args:
        r: float or array of floats with distances between atoms.
        epsilon: Characteristic energy, set to 1 as default.
        sigma: Characteristic length, set to 1 as default.

    Returns:
        The computed potential between the atoms given as input.
    """

    return 4*epsilon*((sigma/r)**12 - (sigma/r)**6)

def LJ_force(r, epsilon=1, sigma=1):
    """Function for computing the forces of the Lennard-Jones Potential.

    Args:
        r: float or array of floats with distances between atoms.
        epsilon: Characteristic energy, set to 1 as default.
        sigma: Characteristic length, set to 1 as default.

    Returns:
        The computed forces between the atoms given as input.
    """

    return 24*epsilon/r*(2*(sigma/r)**12 - (sigma/r)**6)


def plot_LJ():
    """ Simple function for plotting the curves of the LJ potential and forces.
    Simply for showing the shape of the curves, the r-array is a set of arbitrary r-values in an
    interval near y = 0.
    """

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



def simulate(r0, v0, m, T, dt, intg, wrt_file = False):

    """ A function for simulating the 2-atom model with an integration method of choice.

    Args:
        r0: Array with the initial condition for positions.
        v0: Array with the initial condition for velocities.
        m: mass of the atoms
        T: The end time where the integration will stop.
        dt: Size of the time step.
        intg: String telling what method to use; either 'Euler', 'EulerCromer' og 'VelVerlet'.
        wrt_file: Whether or not to write the position vector to an .xyz file. Default set to 'False'.

    returns:
        r: Array for posisions.
        v: Array for velocities.
        t: Array for the discrete time points.
    """

    nt = int(T/dt) + 1   # Number of timesteps
    t = np.linspace(0, T, nt)

    r = np.zeros((nt, 6))
    r[0] = np.array(r0)

    v = np.zeros((nt, 6))
    v[0] = np.array(v0)


    # Time integration
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

    if wrt_file:
        write_xyz(r)

    return r, v, t


def Energy(r, v, t, m, plot = True):
    """ Calculate both the kinetic and potential energies of the 2-atom model

    Args:
        r: Array for positions.
        v: Array for velocities.
        t: Array for the discrete time points.
        m: Mass of the atoms.
        plot: Boolean expression, whether to plot the results or not, 'True' as default.

    Returns:
        KinEng: Array for the kinetic energy at each timestep.
        PotEng: Array for the potential energy at each timestep.
    """

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

def write_xyz(r):
    """Write positions to .xyz file for visualization in external program (i.e. Ovito).

    Args:
        r: Array with positions for all atoms at all timesteps
    """

    nt = r.shape[0]
    N = int(r.shape[1]/3)

    infile = open("in.xyz", 'w')

    for i in range(nt):
        infile.write("%i\n\n" %(N))
        for j in range(N):
            infile.write(str(r[i, 3*j:3*j +3])[1:-1] + "\n")

    infile.close()
