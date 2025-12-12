#Importando as bibliotecas
from numba import jit
import numpy as np
from constants import *
from functions import *

#Definindo o tamanho do grid
Nx = 257
Ny = 257
grid_num = Nx*Ny

#Definindo tempo de simulação.
tf = 100
t_interval = 10

#Cálculo de tau.
Re = 100
u_max = 1
delta_t = 1

ni = u_max*Ny/Re

tau = ni * a_s ** 2 + delta_t/2

omega = tau**(-1)

#Definindo velocidade e densidade inicial.
rho = np.ones(grid_num)
ux = np.zeros(grid_num)
uy = np.zeros(grid_num)

rho_plot = []
ux_plot = []
uy_plot = []

#Cálculo inicial de f de equilíbrio.
f_eq = np.zeros(Nx*Ny*Q)

for x in range(Nx):
    for y in range(Ny):
        for i in range(Q):
            g_id = grid_id(x, y, Nx)
            f_eq[pop_id(x, y, i, Nx, Q)] = equilibrium(rho[g_id], w_i[i], a_s, ux[g_id], 
                                                       c_ix[i], uy[g_id], c_iy[i])

for y in range(Ny):
    for x in range(Nx):
        print(x, y)
        if x==0 and y==0:
            print("Sudoeste")
            I_s = [0, 1, 2, 5]
            O_s = [0, 1, 2, 3, 4, 6, 7, 8]
        elif x==(Nx-1) and y==0:
            print("Sudeste")
            I_s = [0, 2, 3, 6]
            O_s = [0, 1, 2, 3, 4, 5, 7, 8]
        elif x==0 and y==(Ny-1):
            print("Noroeste")
            I_s = [0, 1, 4, 8]
            O_s = [0, 1, 2, 3, 4, 5, 6, 7]
        elif x==(Nx-1) and y==(Ny-1):
            print("Nordeste")
            I_s = [0, 3, 4, 7]
            O_s = [0, 1, 2, 3, 4, 5, 6, 8]
        elif y==0:
            print("Sul")
            I_s = [0, 1, 2, 3, 5, 6]
            O_s = [0, 1, 3, 4, 7, 8]
        elif y==(Ny-1):
            print("Norte")
            I_s = [0, 1, 3, 4, 7, 8]
            O_s = [0, 1, 2, 3, 5, 6]
        elif x==0:
            print("Oeste")
            I_s = [0, 1, 2, 4, 5, 8]
            O_s = [0, 2, 3, 4, 6, 7]
        elif x==(Nx-1):
            print("Leste")
            I_s = [0, 2, 3, 4, 6, 7]
            O_s = [0, 1, 2, 4, 5, 8]
        else:
            print("Centro")
            I_s = [0, 1, 2, 3, 4, 5, 6, 7, 8]