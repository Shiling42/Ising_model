import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pylab

#decide which spin to flip
def gray_flip(t, N):
    k = t[0]
    if k > N: return t, k
    t[k - 1] = t[k]
    t[k] = k + 1
    if k != 1: t[0] = 1
    return t, k#k is the index of spin, t is 


def ising(T,L):
    M_ave=[]
    E2_ave=[]
    Cv=[]
    E_ave=[]
    for t in list(T):
        N = L * L
        nbr = {i : ((i // L) * L + (i + 1) % L, (i + L) % N,
                    (i // L) * L + (i - 1) % L, (i - L) % N)
                                            for i in list(range(N))} 
        #find nearst neighbors
        beta=1/t
        S = [-1] * N#configuration of spins
        E = -2 * N#total energy
        En=[E]
        p=[np.exp(-beta*E)]
        M=[np.sum(S)/N]#average spin
        tau = list(range(1, N + 2))  
        for i in list(range(1, 2 ** N)):
            tau, k = gray_flip(tau, N)
            h = sum(S[n] for n in nbr[k - 1])#nearst neighbor of kth spin
            E += 2 * h * S[k - 1] #energy
            S[k - 1] *= -1#flip the spin
            p.append(np.exp(-beta*E))
            M.append(sum(S)/N)
            En.append(E)   
        Z=np.sum(p)
        M_ave.append(np.dot(np.abs(M),p)/Z)
        E_temp=np.dot(En,p)/Z/N
        E2_ave=np.dot(np.power(En,2),p)/Z/N**2
        Cv.append((E2_ave-E_temp**2)*beta**2)
        E_ave.append(E_temp)
        #print(M_ave,E_ave)
    return M_ave,E_ave,Cv


if __name__ == '__main__':
    list_T=list(np.linspace( 0.1,10,50))
    for m in list([2,3,4]):
        list_M_ave=[]
        list_E=[]
        list_Cv=[]
        #for i in list(list_T):
        list_M,list_E,list_Cv=ising(list_T,m)
     #   list_M_ave=np.gradient(list_M_ave,0.12)
        pylab.plot(list_T,  list_Cv, 'k-', lw=1,color=np.random.rand(3,1))
    pylab.xlabel('$T$', fontsize=20)
    pylab.ylabel('$|M|_{av}$', fontsize=20)
    #pylab.show()
    #pylab.savefig('enumaration.png')
    


