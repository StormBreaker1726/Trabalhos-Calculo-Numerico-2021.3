import time
from matplotlib.pyplot import close
from numpy import meshgrid
from scipy.stats import multivariate_normal
from SOR import *
from Jacobi import *
from mpl_toolkits import mplot3d

N = 50

file = open("resultados_dos_testes.txt", 'w')

fator_de_aceleracao = np.array((1.00, 1.95, 1.99))

for i in range(len(fator_de_aceleracao)):
    file.write("\n\n\nResultados obtidos pelo método SOR com w = "+str(fator_de_aceleracao[i])+":\n")
    ini = time.time()
    u, k, e = SOR(N, 0, 1, mt.pow(10,-7), 10000, fator_de_aceleracao[i])
    fim = time.time()
    file.write("Tempo de cálculo = "+str(fim-ini)+" segundos\n")
    file.write("k = "+str(k))
    file.write("\nerro = "+str(e))
    


file.write("\n\n\nResultados obtidos pelo método de Jacobi:\n")
ini = time.time()
u, k, e = jacobi(N, 0, 1, mt.pow(10,-7), 10000)
fim = time.time()
file.write("Tempo de cálculo = "+str(fim-ini)+" segundos\n")
file.write("k = "+str(k))
file.write("\nerro = "+str(e))



file.close()