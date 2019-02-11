#Trabalho 2 de Modelagem de Reservatórios de Petróleo - 2018
#Autor: Leonardo Simões

import math
import matplotlib.pyplot as plt
import scipy.special as bessel

def plotaGrafico(X, P1, P2, P3, legendas, img):
    fig = plt.figure(1)
    plt.xlabel(legendas[0])
    plt.ylabel(legendas[1])
    plt.title(legendas[2])
    plt.plot(X, P1, label=legendas[3])
    plt.plot(X, P2, label=legendas[4])
    plt.plot(X, P3, label=legendas[5])
    plt.legend()
    plt.grid()
    plt.show()
    fig.savefig(img + '.png')

def dadosGrupo4():
    k = 30 * 10**-14
    fi = 0.15
    pin = 35 * 10**6 
    re = 4000
    Lz = 15
    B0 = 1.35
    u = 0.6 * 10**-3
    ct = 10 * 10**-10
    return (k, fi, pin, re, Lz, u, ct, B0)

def Ni(k, fi, u, ct):
    return k/(fi*u*ct)

def Vi():
    V = {}
    V[1]=0.08333
    V[2]=-32.08333
    V[3]=1279
    V[4]=-15623.666
    V[5]=84244.1666
    V[6]=-236957.5
    V[7]=375911.666
    V[8]=-340071.666
    V[9]=164062.5
    V[10]=-32812.5
    return V

def DeltaPs(s, rw, constantes, n, qw):
    (k, fi, pin, re, h, u, ct, B0) = constantes
    sn = math.sqrt(s/n)
    a = qw*u/(2*math.pi*k*h)
    b = s*sn*rw*(bessel.iv(1,sn*re)*bessel.kv(1,sn*rw)-bessel.iv(1,sn*rw)*bessel.kv(1,sn*re))
    c = bessel.iv(0,sn*rw)*bessel.kv(1,sn*re)+bessel.iv(1,sn*re)*bessel.kv(0,sn*rw)
    return a*c/b
   
def solucaoRegimeTransientePseudoPermanente(constantes, deltaT, Tmax, qsc1, qsc2, qsc3, rw1, rw2, rw3, u1, u2, u3):
    (k, fi, pin, re, h, u, ct, B) = constantes
    P1, P2, P3, P4, P5, P6, P7, P8, P9 = [], [], [], [], [], [], [], [], []
    T, LogT = [], []
    V = Vi()
    n = Ni(k, fi, u, ct)
    qw1, qw2, qw3 = qsc1 * B /86400, qsc2 * B /86400, qsc3 * B /86400  
    constantes1, constantes2, constantes3 = list(constantes), list(constantes), list(constantes)
    constantes1[5], constantes2[5], constantes3[5] = u1, u2, u3 
    t = 1000
    while t <= Tmax:
        T.append(t)
        LogT.append(math.log(t, 10))
        VDeltaP1 = VDeltaP2 = VDeltaP3 = VDeltaP4 = VDeltaP5 = VDeltaP6 = VDeltaP7 = VDeltaP8 = VDeltaP9 = 0
        for i in range(1,11):
            s = i*math.log(2)/t
            VDeltaP1 += V[i]*DeltaPs(s, rw1, constantes, n, qw1)
            VDeltaP2 += V[i]*DeltaPs(s, rw1, constantes, n, qw2)
            VDeltaP3 += V[i]*DeltaPs(s, rw1, constantes, n, qw3)
            VDeltaP4 += V[i]*DeltaPs(s, rw1, constantes, n, qw1)
            VDeltaP5 += V[i]*DeltaPs(s, rw2, constantes, n, qw1)
            VDeltaP6 += V[i]*DeltaPs(s, rw3, constantes, n, qw1)
            VDeltaP7 += V[i]*DeltaPs(s, rw1, constantes1, n, qw1)
            VDeltaP8 += V[i]*DeltaPs(s, rw1, constantes2, n, qw1)
            VDeltaP9 += V[i]*DeltaPs(s, rw1, constantes3, n, qw1)
        P1.append(pin - (math.log(2)/t * VDeltaP1))
        P2.append(pin - (math.log(2)/t * VDeltaP2))
        P3.append(pin - (math.log(2)/t * VDeltaP3))
        P4.append(pin - (math.log(2)/t * VDeltaP4))
        P5.append(pin - (math.log(2)/t * VDeltaP5))
        P6.append(pin - (math.log(2)/t * VDeltaP6))
        P7.append(pin - (math.log(2)/t * VDeltaP7))
        P8.append(pin - (math.log(2)/t * VDeltaP8))
        P9.append(pin - (math.log(2)/t * VDeltaP9))
        deltaT = deltaT*1.5
        t= t + deltaT
    legendas1 = ('log 10 (Tempo)', 'Pressão (MPa)', 'Pressão x Logaritmo do Tempo', 'qsc='+str(qsc1), 'qsc='+str(qsc2), 'qsc='+str(qsc3))
    plotaGrafico(LogT, P1, P2, P3, legendas1, 'Gráfico 1')
    legendas2 = ('log 10 (Tempo)', 'Pressão (MPa)', 'Pressão x Logaritmo do Tempo', 'rw='+str(rw1), 'rw='+str(rw2), 'rw='+str(rw3))
    plotaGrafico(LogT, P4, P5, P6, legendas2, 'Gráfico 2')
    legendas3 = ('log 10 (Tempo)', 'Pressão (MPa)', 'Pressão x Logaritmo do Tempo', 'u='+str(constantes1[5]), 'u='+str(constantes2[5]), 'u='+str(constantes3[5]))
    plotaGrafico(LogT, P7, P8, P9, legendas3, 'Gráfico 3')

#Função Principal - MAIN        
if __name__ == '__main__':
    constantes = dadosGrupo4()
    Tmax = 87550000
    deltaT = 86400
    qsc1, qsc2, qsc3 = 500, 600, 700
    rw1, rw2, rw3 = 0.20, 0.45, 0.6
    u1, u2, u3 = 0.4 * 10**-3, 0.5 * 10**-3, constantes[5]
    solucaoRegimeTransientePseudoPermanente(constantes, deltaT, Tmax, qsc1, qsc2, qsc3, rw1, rw2, rw3, u1, u2, u3)
