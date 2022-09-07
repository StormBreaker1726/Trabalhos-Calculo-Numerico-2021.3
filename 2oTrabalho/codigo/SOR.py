#Implementação seguindo o método SOR
import math as mt
from matplotlib import projections
import sympy as sp
from mpl_toolkits import mplot3d
import matplotlib as mp
from matplotlib import pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import Gnuplot

'''
n ==> número de pontos em cada direção
(a,b) ==> intervalo
tol ==> tolerância
k_Max ==> número máximo de iterações
fator_de_acleração ==> fator de aceleração a ser usado
'''

#seguindo o algoritmo esplicitado no roteiro da atividade

def SOR(n, a, b, tol, k_Max, fator_de_aceleracao):
    #alocação das matrizes u_n+1 e u_old e sua inicialização com zeros
    u = np.zeros((n, n))
    u_old = np.zeros((n,n)) 
    xy = np.linspace(a, b, n)
    
    [xx, yy] = np.meshgrid(xy, xy)

    h = (b - a)/(n - 1) #definição de h

    e = 1 #definição de e

    k = 0  #definição de k

    #Aplicação das condições de contorno
    for i in range(0, n):
        u[i][n-1] = 1
    
    while(e >= tol and k <= k_Max):
        k = k + 1 #atualização do k
        e = 0

        u_old = np.copy(u) #cópia da matriz

        for j in range(1, n-1):
            y = j*(h) 
            for i in range(1, n-1):
                x = i*(h) 
                S = 10*(x*x+y*y+5) #Cálculo do S
                u[i][j] = (1 - fator_de_aceleracao)*u_old[i][j] + (1/4)*fator_de_aceleracao*(u[i-1][j]+u_old[i+1][j]+u[i][j-1]+u_old[i][j+1]+mt.pow(h,2)*S)

                e = e+abs((u[i][j] - u_old[i][j]))
        e = e/((n-2)*(n-2))

    firg = plt.figure()
    ax = firg.gca(projection='3d')
    surf = ax.plot_surface(xx, yy, u)
    plt.title("Método SOR com fator de aceleração " +str(fator_de_aceleracao)+":")
    # plt.show()
    plt.savefig("Método SOR com fator de aceleração " +str(fator_de_aceleracao)+".png")
    return u, k, e #retorna a matriz da solução, o número de iterações e o erro
