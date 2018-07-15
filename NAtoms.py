import numpy as np
import matplotlib.pylab as plt
import sys

"""
This file contains all functions necessary to solve to the tasks for the 2-atom Lennard-Jones
problem. The tasks themselves will be solved in TwoAtoms_tasks.py, in which all the function
in this file will be imported.
"""

def LJ_pot(r, epsilon=1, sigma=1):
    """Function for computing the potential of the Lennard-Jones Potential. Cut-off at r = 3.0.

    Args:
        r: float or array of floats with distances between atoms.
        epsilon: Characteristic energy, set to 1 as default.
        sigma: Characteristic length, set to 1 as default.

    Returns:
        The computed potential between the atoms given as input.
    """
    if r < 3.0:
        U_corr = 4*epsilon*((sigma/3.0)**12 - (sigma/3.0)**6)
        return 4*epsilon*((sigma/r)**12 - (sigma/r)**6) - U_corr
    else:
        return 0

def LJ_force(r, epsilon=1, sigma=1):
    """Function for computing the forces of the Lennard-Jones Potential.

    Args:
        r: float or array of floats with distances between atoms.
        epsilon: Characteristic energy, set to 1 as default.
        sigma: Characteristic length, set to 1 as default.

    Returns:
        The computed forces between the atoms given as input.
    """
    if r < 3.0:
        return 24*epsilon/r*(2*(sigma/r)**12 - (sigma/r)**6)
    else:
        return 0


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

    """ A function for simulating the N-atom model with an integration method of choice. Stores all
    positions and velocities at every time steps in arrays. May run into memory issues.

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

    N = int(r0.shape[0]/3)

    r = np.zeros((nt, 3*N))
    r[0] = np.array(r0)

    v = np.zeros((nt, 3*N))
    v[0] = np.array(v0)

    # Time integration
    for i in range(nt-1):

        if intg == 'VelVerlet' and i != 0:
            a = a_next

        else:

            F = np.zeros(3*N)

            for j in range(N):
                for k in range(j+1, N):
                    d_vec = r[i, 3*j: 3*j+3] - r[i, 3*k: 3*k+3]
                    d = np.sqrt(np.sum((d_vec)**2))

                    F_jk = LJ_force(d)*d_vec/d
                    F[3*j: 3*j+3] += F_jk
                    F[3*k: 3*k+3] -= F_jk

            a = F/m

        if intg == 'Euler':
            v[i+1] = v[i] + a*dt
            r[i+1] = r[i] + v[i]*dt

        elif intg == 'EulerCromer':
            v[i+1] = v[i] + a*dt
            r[i+1] = r[i] + v[i+1]*dt

        elif intg == 'VelVerlet':
            r[i+1] = r[i] + v[i]*dt + 0.5*a*dt**2

            F_next = np.zeros(3*N)

            for j in range(N):
                for k in range(j+1, N):
                    d_vec = r[i+1, 3*j: 3*j+3] - r[i+1, 3*k: 3*k+3]
                    d = np.sqrt(np.sum((d_vec)**2))

                    F_jk = LJ_force(d)*d_vec/d
                    F_next[3*j: 3*j+3] += F_jk
                    F_next[3*k: 3*k+3] -= F_jk

            a_next = F_next/m

            v[i+1] = v[i] + 0.5*(a + a_next)*dt

        else:
            print('Error: The variable intg must either be "Euler", "EulerCromer" or "VelVerlet".\nExiting...')
            sys.exit(0)

    if wrt_file:
        write_xyz(r)

    return r, v, t

def simulate_sparse(r0, v0, m, T, dt, intg, filename):

    """ A function for simulating the N-atom model with an integration method of choice. Same as simulate(),
    but doesn't store positions and velocities in arrays, but writes them to file instead.

    Args:
        r0: Array with the initial condition for positions.
        v0: Array with the initial condition for velocities.
        m: mass of the atoms
        T: The end time where the integration will stop.
        dt: Size of the time step.
        intg: String telling what method to use; either 'Euler', 'EulerCromer' og 'VelVerlet'.
        filename: Name of the file the data will be written to. Needs an extension like ".xyz".

    returns:
        t: Array for the discrete time points.
    """

    nt = int(T/dt) + 1   # Number of timesteps
    t = np.linspace(0, T, nt)

    N = int(r0.shape[0]/3)

    r = r0
    v = v0

    infile = open(filename, 'w')
    infile.write("%i\n\n" %(N))

    for j in range(N):
        infile.write(str(r[3*j: 3*j +3])[1:-1] + " " + str(v[3*j : 3*j +3])[1:-1] + "\n")


    # Time integration
    for i in range(nt-1):

        infile.write("%i\n\n" %(N))

        if intg == 'VelVerlet' and i != 0:
            a = a_next

        else:

            F = np.zeros(3*N)

            for j in range(N):
                for k in range(j+1, N):
                    d_vec = r[3*j: 3*j+3] - r[3*k: 3*k+3]
                    d = np.sqrt(np.sum((d_vec)**2))

                    F_jk = LJ_force(d)*d_vec/d
                    F[3*j: 3*j+3] += F_jk
                    F[3*k: 3*k+3] -= F_jk

            a = F/m

        if intg == 'Euler':
            v_next = v + a*dt
            r_next = r + v*dt

        elif intg == 'EulerCromer':
            v_next = v + a*dt
            r_next = r + v_next*dt

        elif intg == 'VelVerlet':
            r_next = r + v*dt + 0.5*a*dt**2

            F_next = np.zeros(3*N)

            for j in range(N):
                for k in range(j+1, N):
                    d_vec = r_next[3*j: 3*j+3] - r_next[3*k: 3*k+3]
                    d = np.sqrt(np.sum((d_vec)**2))

                    F_jk = LJ_force(d)*d_vec/d
                    F_next[3*j: 3*j+3] += F_jk
                    F_next[3*k: 3*k+3] -= F_jk

            a_next = F_next/m

            v_next = v + 0.5*(a + a_next)*dt

        else:
            print('Error: The variable intg must either be "Euler", "EulerCromer" or "VelVerlet".\nExiting...')
            sys.exit(0)

        for j in range(N):
            infile.write(str(r_next[3*j : 3*j +3])[1:-1] + " " + str(v_next[3*j : 3*j +3])[1:-1]+ "\n")

        r = r_next
        v = v_next

    infile.close()

    return t


def Energy(r, v, t, m, plot = True):
    """ Calculate both the kinetic and potential energies of the N-atom model from the stored
    position and velocity arrays generated in simulate().

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
    N = int(v.shape[1]/3)
    KinEng = np.zeros(nt)
    PotEng = np.zeros(nt)

    for i in range(nt):
        for j in range(N):
            v_abs2 = np.sum(v[i, 3*j:3*j + 3]**2)
            KinEng[i] += 0.5*m*v_abs2

            for k in range(j+1, N):

                d = np.sqrt(np.sum((r[i, 3*j: 3*j +3] - r[i, 3*k: 3*k +3])**2))
                PotEng[i] += LJ_pot(d)

    if plot:
        fig, (ax1, ax2) = plt.subplots(2, 1)
        ax1.set_title("Kinetic and Potential energy of the " + str(N) + "-atom system")
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

def Energy_sparse(N, t, m, filename, plot = True):
    """ Calculate both the kinetic and potential energies of the N-atom model, reads positions
    and velocities from file generated by simulate_sparse().

    Args:
        N: Number of atoms
        t: Array for the discrete time points.
        m: Mass of the atoms.
        filename: Name of the file to read data from. Assumes that the file has the form:
            ----------------------------------
            <No. atoms>

            <r_1x r_1y r_1z v_1x v_1y v_1z>   # Positions and velocities of atom 1
            ...
            <r_Nx r_Ny r_Nz v_Nx v_Ny v_Nz>   # Positions and velocities of atom N
            ----------------------------------
            This snippet repeats for all time steps.
        plot: Boolean expression, whether to plot the results or not, 'True' as default.

    Returns:
        KinEng: Array for the kinetic energy at each timestep.
        PotEng: Array for the potential energy at each timestep.
    """
    from itertools import islice


    nt = t.shape[0]
    KinEng = np.zeros(nt)
    PotEng = np.zeros(nt)

    #infile = open(filename, 'r')
    with open(filename, 'r') as infile:
        for i in range(nt):
            r = np.empty(0)
            line_gen = islice(infile, 2, N+2)
            j = 0
            for line in line_gen:
                tmp = np.array(line.split(), dtype=np.float)
                r = np.concatenate((r, tmp[:3]))
                v = tmp[3:]

                KinEng[i] += 0.5*m*np.sum(v**2)

                for k in range(j):
                    d_vec = r[3*j : 3*j +3] - r[3*k : 3*k + 3]
                    d = np.sqrt(np.sum(d_vec**2))

                    PotEng[i] += LJ_pot(d)

                j += 1

    if plot:
        fig, (ax1, ax2) = plt.subplots(2, 1)
        ax1.set_title("Kinetic and Potential energy of the " + str(N) + "-atom system")
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
