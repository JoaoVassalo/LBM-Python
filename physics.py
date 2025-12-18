from geometry import *
from stencil import *
import numpy as np

tf = 10000 #Tempo de simulação.
t_interval = 100 #Intervalo entre cada frame.

#Propriedades físicas.
Re = 1000 #Reynolds.
u_max = 0.0256 #Velocidade da tampa.
delta_t = 1 #Variação de tempo entre cada step.

#Cálculo de tau.
ni = u_max*Ny/Re
tau = ni * a_s ** 2 + delta_t/2
omega = tau**(-1)

