import numpy as np
import matplotlib.pylab as plt
from TwoAtoms import *

def test():

    m = 1   # Mass set to 1 for simplicity

    r, v, t = simulate(m, 10, 0.01, 'VelVerlet')

    #d_abs = np.sqrt((np.sum((r[:, :3] - r[:, 3:])**2, axis=1)))

    #plt.plot(t, d_abs)
    #plt.show()

    KinEng, PotEng = Energy(r, v, t, m)

    plt.plot(t, (KinEng + PotEng))
    plt.show()

if __name__ == "__main__":
    test()
