#Importando as bibliotecas
import numpy as np
from physics import *
from stencil import *
from geometry import *
from functions import *
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

rho, ux, uy, mxx, myy, mxy, rho_plot, ux_plot, uy_plot = init_domain()



for t in range(tf):
    

    #Colisão.
    for x in range(Nx):
        for y in range(Ny):
            g_id = grid_id(x, y)
            for i in range(Q):
                p_id = pop_id(g_id, i)
                col[p_id] = - (f_i[p_id] - f_eq[p_id])/tau
                f_col[p_id] = f_i[p_id] + col[p_id]

    #Propagação.
    for x in range(Nx):
        for y in range(Ny):
            for i in range(Q):

                x_to = (x + c_ix[i] + Nx) % Nx   #Contornos periódicos.
                y_to = (y + c_iy[i] + Ny) % Ny


                idx_to   = pop_id(grid_id(x_to, y_to), i)
                idx_from = pop_id(grid_id(x, y), i)

                f_i[idx_to] = f_col[idx_from]

    for x in range(Nx):
        for y in range(Ny):
            g_id = grid_id(x, y)

            if x==0 and y==0:   #Sudoeste
                I_s = [0, 3, 4, 7]
                O_s = [0, 1, 2, 5]
                rho[g_id] = rho_corner(f_i, I_s, O_s, g_id)
                ux[g_id] = 0
                uy[g_id] = 0
                mxy = 0

                for i in range(Q):
                    p_id = pop_id(g_id, i)
                    f_i[p_id] = w_i[i] * rho[g_id] * (1 + 
                                                      a_s**2 * (ux[g_id]*c_ix[i] + uy[g_id]*c_iy[i]) + 
                                                      a_s**4/2 * (ux[g_id]**2*(c_ix[i]**2 - 1/a_s**2) + uy[g_id]**2*(c_iy[i]**2 - 1/a_s**2)) + 
                                                      a_s**4 * (mxy*c_ix[i]*c_iy[i]) )

            elif x==(Nx-1) and y==0:    #Sudeste
                I_s = [0, 1, 4, 8]
                O_s = [0, 2, 3, 6]
                rho[g_id] = rho_corner(f_i, I_s, O_s, g_id)
                ux[g_id] = 0
                uy[g_id] = 0
                mxy = 0

                for i in range(Q):
                    p_id = pop_id(g_id, i)
                    f_i[p_id] = w_i[i] * rho[g_id] * (1 + 
                                                      a_s**2 * (ux[g_id]*c_ix[i] + uy[g_id]*c_iy[i]) + 
                                                      a_s**4/2 * (ux[g_id]**2*(c_ix[i]**2 - 1/a_s**2) + uy[g_id]**2*(c_iy[i]**2 - 1/a_s**2)) + 
                                                      a_s**4 * (mxy*c_ix[i]*c_iy[i]) )

            elif x==0 and y==(Ny-1):    #Noroeste
                I_s = [0, 2, 3, 6]
                O_s = [0, 1, 4, 8]
                ux[g_id] = u_max
                uy[g_id] = 0
                mxy_I = m_xy_I(f_i, I_s, g_id)
                mxy = m_xy_north(mxy_I, ux[g_id], I_s, O_s, omega)
                rho[g_id] = rho_north(mxy, f_i, ux[g_id], I_s, O_s, omega, g_id)

                for i in range(Q):
                    p_id = pop_id(g_id, i)
                    f_i[p_id] = w_i[i] * rho[g_id] * (1 + 
                                                      a_s**2 * (ux[g_id]*c_ix[i] + uy[g_id]*c_iy[i]) + 
                                                      a_s**4/2 * (ux[g_id]**2*(c_ix[i]**2 - 1/a_s**2) + uy[g_id]**2*(c_iy[i]**2 - 1/a_s**2)) + 
                                                      a_s**4 * (mxy*c_ix[i]*c_iy[i]) )

            elif x==(Nx-1) and y==(Ny-1):   #Nordeste
                I_s = [0, 1, 2, 5]
                O_s = [0, 3, 4, 7]
                ux[g_id] = u_max
                uy[g_id] = 0
                mxy_I = m_xy_I(f_i, I_s, g_id)
                mxy = m_xy_north(mxy_I, ux[g_id], I_s, O_s, omega)
                rho[g_id] = rho_north(mxy, f_i, ux[g_id], I_s, O_s, omega, g_id)

                for i in range(Q):
                    p_id = pop_id(g_id, i)
                    f_i[p_id] = w_i[i] * rho[g_id] * (1 + 
                                                      a_s**2 * (ux[g_id]*c_ix[i] + uy[g_id]*c_iy[i]) + 
                                                      a_s**4/2 * (ux[g_id]**2*(c_ix[i]**2 - 1/a_s**2) + uy[g_id]**2*(c_iy[i]**2 - 1/a_s**2)) + 
                                                      a_s**4 * (mxy*c_ix[i]*c_iy[i]) )


            elif y==0:  #Sul
                I_s = [0, 1, 3, 4, 7, 8]
                O_s = [0, 1, 2, 3, 5, 6]
                mxy_I = m_xy_I(f_i, I_s, g_id)
                ux[g_id] = 0
                uy[g_id] = 0
                mxy = m_xy_wall(mxy_I, I_s, O_s, omega)
                rho[g_id] = rho_wall(mxy, f_i, I_s, O_s, omega, g_id)

                for i in range(Q):
                    p_id = pop_id(g_id, i)
                    f_i[p_id] = w_i[i] * rho[g_id] * (1 + 
                                                      a_s**2 * (ux[g_id]*c_ix[i] + uy[g_id]*c_iy[i]) + 
                                                      a_s**4/2 * (ux[g_id]**2*(c_ix[i]**2 - 1/a_s**2) + uy[g_id]**2*(c_iy[i]**2 - 1/a_s**2)) + 
                                                      a_s**4 * (mxy*c_ix[i]*c_iy[i]) )

            elif y==(Ny-1): #Norte
                I_s = [0, 1, 2, 3, 5, 6]
                O_s = [0, 1, 3, 4, 7, 8]
                mxy_I = m_xy_I(f_i, I_s, g_id)
                ux[g_id] = u_max
                uy[g_id] = 0
                mxy = m_xy_north(mxy_I, ux[g_id], I_s, O_s, omega)
                rho[g_id] = rho_north(mxy, f_i, ux[g_id], I_s, O_s, omega, g_id)

                for i in range(Q):
                    p_id = pop_id(g_id, i)
                    f_i[p_id] = w_i[i] * rho[g_id] * (1 + 
                                                      a_s**2 * (ux[g_id]*c_ix[i] + uy[g_id]*c_iy[i]) + 
                                                      a_s**4/2 * (ux[g_id]**2*(c_ix[i]**2 - 1/a_s**2) + uy[g_id]**2*(c_iy[i]**2 - 1/a_s**2)) + 
                                                      a_s**4 * (mxy*c_ix[i]*c_iy[i]) )

            elif x==0:  #Oeste
                I_s = [0, 2, 3, 4, 6, 7]
                O_s = [0, 1, 2, 4, 5, 8]
                mxy_I = m_xy_I(f_i, I_s, g_id)
                ux[g_id] = 0
                uy[g_id] = 0
                mxy = m_xy_wall(mxy_I, I_s, O_s, omega)
                rho[g_id] = rho_wall(mxy, f_i, I_s, O_s, omega, g_id)

                for i in range(Q):
                    p_id = pop_id(g_id, i)
                    f_i[p_id] = w_i[i] * rho[g_id] * (1 + 
                                                      a_s**2 * (ux[g_id]*c_ix[i] + uy[g_id]*c_iy[i]) + 
                                                      a_s**4/2 * (ux[g_id]**2*(c_ix[i]**2 - 1/a_s**2) + uy[g_id]**2*(c_iy[i]**2 - 1/a_s**2)) + 
                                                      a_s**4 * (mxy*c_ix[i]*c_iy[i]) )

            elif x==(Nx-1): #Leste
                I_s = [0, 1, 2, 4, 5, 8]
                O_s = [0, 2, 3, 4, 6, 7]
                mxy_I = m_xy_I(f_i, I_s, g_id)
                ux[g_id] = 0
                uy[g_id] = 0
                mxy = m_xy_wall(mxy_I, I_s, O_s, omega)
                rho[g_id] = rho_wall(mxy, f_i, I_s, O_s, omega, g_id)

                for i in range(Q):
                    p_id = pop_id(g_id, i)
                    f_i[p_id] = w_i[i] * rho[g_id] * (1 + 
                                                      a_s**2 * (ux[g_id]*c_ix[i] + uy[g_id]*c_iy[i]) + 
                                                      a_s**4/2 * (ux[g_id]**2*(c_ix[i]**2 - 1/a_s**2) + uy[g_id]**2*(c_iy[i]**2 - 1/a_s**2)) + 
                                                      a_s**4 * (mxy*c_ix[i]*c_iy[i]) )

            else:   #Centro
                rho[g_id] = calc_rho(f_i, g_id)
                ux[g_id], uy[g_id] = calc_velocity(f_i, g_id, rho[g_id])



    for x in range(Nx):
        for y in range(Ny):
            for i in range(Q):
                g_id = grid_id(x, y)
                f_eq[pop_id(g_id, i)] = equilibrium(rho[g_id], w_i[i], a_s, ux[g_id], 
                                                        c_ix[i], uy[g_id], c_iy[i])


    if t % t_interval == 0:
        print("---------------------------------------", t, "---------------------------------------")
        rho_plot.append(rho.copy())
        ux_plot.append(ux.copy())
        uy_plot.append(uy.copy())

# Cria a figura e o eixo
fig, ax = plt.subplots()
img = ax.imshow(np.zeros((Ny, Nx)), origin='lower', cmap='viridis',
                vmin=np.min(rho_plot), vmax=np.max(rho_plot))
ax.set_title("Densidade ρ ao longo do tempo")
ax.set_xlabel("x")
ax.set_ylabel("y")

cb = plt.colorbar(img, ax=ax)
cb.set_label("ρ")

# Função que atualiza cada frame
def update(t):
    rho_t = np.array(rho_plot[t]).reshape((Ny, Nx))
    img.set_data(rho_t)
    ax.set_title(f"Densidade ρ — passo {t}")
    return [img]

# Cria a animação
ani = FuncAnimation(fig, update, frames=len(rho_plot), interval=100, blit=True)

# Caminho onde salvar o gif
output_path = r"C:\Users\08552591910\Pictures\Attempt_2.gif"

# Salva o GIF (usa PillowWriter)
writer = PillowWriter(fps=10)  # ajuste fps se quiser
ani.save(output_path, writer=writer, dpi=150)

print(f"GIF salvo em: {output_path}")
plt.close()

# Cria a figura e eixo
fig, ax = plt.subplots(figsize=(6,6))

# Inicializa com zeros
img = ax.imshow(np.zeros((Ny, Nx)), origin='lower', cmap='viridis',
                vmin=0, vmax=u_max)  # ajusta vmin/vmax se souber os limites
ax.set_title("Magnitude da velocidade |u|")
ax.set_xlabel("x")
ax.set_ylabel("y")
cb = plt.colorbar(img, ax=ax)
cb.set_label("|u|")

# Função de atualização
def update(t):
    U = np.array(ux_plot[t]).reshape((Ny, Nx))
    V = np.array(uy_plot[t]).reshape((Ny, Nx))
    magnitude = np.sqrt(U**2 + V**2)
    
    img.set_data(magnitude)
    ax.set_title(f"Magnitude da velocidade — passo {t}")
    return [img]

# Cria a animação
ani = FuncAnimation(fig, update, frames=len(rho_plot), interval=100, blit=True)

# Salva em GIF
output_path = r"C:\Users\08552591910\Pictures\Velocity_Field_Magnitude.gif"
writer = PillowWriter(fps=10)
ani.save(output_path, writer=writer, dpi=150)

print(f"GIF salvo em: {output_path}")
plt.close()