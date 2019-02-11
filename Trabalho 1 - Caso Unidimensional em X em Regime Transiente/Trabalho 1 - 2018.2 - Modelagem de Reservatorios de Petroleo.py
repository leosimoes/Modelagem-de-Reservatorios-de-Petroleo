#Trabalho 1 de Modelagem de Reservatórios de Petróleo - 2018
#Autor: Leonardo Simões

import math
import matplotlib.pyplot as plt

def plotaGrafico(X, P1, P2, P3, legendas):
    plt.xlabel(legendas[0])
    plt.ylabel(legendas[1])
    plt.title(legendas[2])
    plt.plot(X, P1, label=legendas[3])
    plt.plot(X, P2, label=legendas[4])
    plt.plot(X, P3, label=legendas[5])
    plt.legend()
    plt.grid()
    plt.show()

def dadosGrupo4():
    k = 30 * 10**-14
    fi = 0.15
    pin = 35 * 10**6 
    Ly = 80
    Lz = 15
    A = Ly * Lz
    B0 = 1.35
    u = 0.6 * 10**-3
    ct = 10 * 10**-10
    return (k, fi, pin, A, u, ct, B0)

def Ni(k, fi, u, ct):
    return k/(fi*u*ct)

def pressaoTransiente(x, t, constantes, qw):
    (k, fi, pin, A, u, ct, B) = constantes
    ni = Ni(k, fi, u, ct)
    return (pin - (qw*u)*(math.sqrt(4*ni*t/math.pi)*math.exp((-x*x)/(4*ni*t))-x*math.erfc(x/math.sqrt(4*ni*t)))/(k*A))/(10**6)

def solucaoRegimeTransiente(qsc1,qsc2,qsc3, t1, t2, t3, Xmax, deltaX, Tmax, deltaT, constantes):
    P1, P2, P3, P4, P5, P6, P7, P8, P9 = [], [], [], [], [], [], [], [], []
    qw1, qw2, qw3 = qsc1 * B /86400, qsc2 * B /86400, qsc3 * B /86400 
    X = [x for x in range(0,Xmax+deltaX,deltaX)]
    for x in X:
        P1.append(pressaoTransiente(x,t1,constantes,qw1))
        P2.append(pressaoTransiente(x,t2,constantes,qw1))
        P3.append(pressaoTransiente(x,t3,constantes,qw1))
        P4.append(pressaoTransiente(x,t1,constantes,qw1))
        P5.append(pressaoTransiente(x,t1,constantes,qw2))
        P6.append(pressaoTransiente(x,t1,constantes,qw3))
    T = [t for t in range(1,Tmax+deltaT,deltaT)]
    x = 50
    for t in T:
        P7.append(pressaoTransiente(x,t,constantes,qw1))
        P8.append(pressaoTransiente(x,t,constantes,qw2))
        P9.append(pressaoTransiente(x,t,constantes,qw3))
    legendas1 = ('Distância (m)', 'Pressão (MPa)', 'Pressão x Distância', 't='+str(t1), 't='+str(t2), 't='+str(t3))
    plotaGrafico(X, P1, P2, P3, legendas1)
    legendas2 = ('Distância (m)', 'Pressão (MPa)', 'Pressão x Distância', 'qsc='+str(qsc1), 'qsc='+str(qsc2), 'qsc='+str(qsc3))
    plotaGrafico(X, P4, P5, P6, legendas2)
    legendas3 = ('Tempo', 'Pressão (MPa)', 'Pressão x Tempo', 'qsc='+str(qsc1), 'qsc='+str(qsc2), 'qsc='+str(qsc3))
    plotaGrafico(T, P7, P8, P9, legendas3)

#Função Principal - MAIN        
if __name__ == '__main__':
    constantes = dadosGrupo4()
    B = constantes[6]
    qsc1, qsc2, qsc3 = 500, 600, 700
    t1, t2, t3 = 10000, 20000, 30000
    Xmax = 1500
    deltaX = 5
    Tmax = 870000
    deltaT = 86400
    solucaoRegimeTransiente(qsc1,qsc2,qsc3, t1, t2, t3, Xmax, deltaX, Tmax, deltaT, constantes)
