import numpy as np
import matplotlib.pylab as plt


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

def simulate(T, dt):

    m = 1   # Mass set to 1 for simplicity
    nt = int(T/dt) + 1   # Number of timesteps

    r = np.zeros((nt, 6))
    r[0] = np.random.uniform(0, 5, size=6)

    v = np.zeros((nt, 6))
    v[0] = np.random.uniform(0, 1, size=6)

    infile = open("in.xyz", 'w')

    # Time integration, Euler-Cromer temporarily
    for i in range(nt-1):

        F = LJ_force(r[i, :3] - r[i, 3:])

        # Forces between the two atoms are the same, but with opposite signs
        F = np.concatenate((F, -F))

        a = F/m

        v[i+1] = v[i] + a*dt
        r[i+1] = r[i] + v[i+1]*dt
        infile.write("2\n\n")
        for j in range(2):
            infile.write(str(r[i, 3*j:3*j +3])[1:-1] + "\n")

        #infile.write("\n")

    infile.close()


if __name__ == "__main__":
    #plot_LJ()
    simulate(4, 0.1)
