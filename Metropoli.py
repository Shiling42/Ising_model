import random, math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pylab
from matplotlib.legend_handler import HandlerLine2D

global L
L=4
global N
N=L*L
S = [random.choice([1, -1]) for k in list(range(N))]#initialize spin configuration

 
def markovising(T,S):
    nsteps = 1000
    nbr = {i : ((i // L) * L + (i + 1) % L, (i + L) % N,
                (i // L) * L + (i - 1) % L, (i - L) % N) \
                                        for i in list(range(N))}
    beta = 1.0 / T
    for step in list(range(nsteps)):
        k = random.randint(0, N - 1)
        delta_E = 2.0 * S[k] * sum(S[nn] for nn in nbr[k])
        if random.uniform(0.0, 1.0) < math.exp(-beta * delta_E):
            S[k] *= -1
    M=np.abs(sum(S)/N)
    return M,S

#main
if __name__ == '__main__':
	list_T=np.linspace( 6,0.1,10)
	list_M_ave=[]
	repeat=300
	for i in list(list_T):
		M_av=0
		M_temp=[]
		for k in list(range(repeat)):
			S1 = [random.choice([1, -1]) for k in list(range(N))]
			M_temp,S=markovising(i,S)
			M_av+=M_temp
		M_av=M_av/(k+1)
		list_M_ave.append(M_av)
    plt.figure(1,figsize=(15,10), dpi=200)
    line1,=pylab.plot(list_T,  list_M_ave, 'k.',label="Metropolis algorithm")
    lin2,=pylab.plot(A,B, 'k-', lw=1,color=np.random.rand(3,1),label="Exact calculation")
    plt.title('Average Magnetazition per Spin($%i\\times%i$ lattice)' % (L, L))
    plt.xlabel('$T$', fontsize=20)
    plt.ylabel(r'$\langle|m|\rangle$', fontsize=20)
    #plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.legend(handler_map={line1: HandlerLine2D(numpoints=4)})
    #pylab.show()
    plt.savefig('enumaration.png')

