from stencil import *
from geometry import *
from physics import *

def init_domain():
    #Definindo velocidades, densidade e momentos iniciais.
    rho = np.ones(grid_num)
    ux = np.zeros(grid_num)
    uy = np.zeros(grid_num)
    mxx = np.zeros(grid_num)
    myy = np.zeros(grid_num)
    mxy = np.zeros(grid_num)

    rho_plot = []
    ux_plot = []
    uy_plot = []

    for x in range(Nx):
        for y in range(Ny):
            g_id = grid_id(x, y)
            for i in range(Q):
                mxx[g_id] += equilibrium(rho[g_id], ux[g_id], uy[g_id], i)*(c_ix[i]**2 - 1/(a_s**2))
                myy[g_id] += equilibrium(rho[g_id], ux[g_id], uy[g_id], i)*(c_iy[i]**2 - 1/(a_s**2))
                mxy[g_id] += equilibrium(rho[g_id], ux[g_id], uy[g_id], i)*(c_ix[i]*c_iy[i])
            mxx[g_id] *= 1/rho[g_id]
            myy[g_id] *= 1/rho[g_id]
            mxy[g_id] *= 1/rho[g_id]
            
    return rho, ux, uy, mxx, myy, mxy, rho_plot, ux_plot, uy_plot

def grid_id(x, y):
    return (x + Nx*y)

def pop_id(g_id, i):
    return (g_id*Q + i)

def equilibrium(rho, ux, uy, i):
    return rho*w_i[i]*( 1 + a_s**2*( ux*c_ix[i] + uy*c_iy[i] ) 
                    + a_s**4/2*( ( ux*c_ix[i] + uy*c_iy[i] )**2 
                                - ( ux**2 + uy**2 )/a_s**2 ) )

def m_xy_I(f_i, I_s, g_id):
    num = 0
    div = 0

    for i in I_s:
        p_id = pop_id(g_id, i)
        num += f_i[p_id] * c_ix[i]*c_iy[i]
        div += f_i[p_id]

    return (num/div)

def m_xy_north(m_xy_I, 
               I_s, 
               O_s):
    Is_up = 0
    Is_down = 0
    Os_down = 0
    Os_up = 0

    for i in I_s:
        Is_up += w_i[i]*c_ix[i]*c_iy[i]*(1 + 
                         a_s**2*u_max*c_ix[i] + 
                         a_s**4/2*u_max*(c_ix[i]**2 - 1/a_s**2))
        Is_down += w_i[i]*a_s**4*c_ix[i]**2*c_iy[i]**2

    for i in O_s:
        Os_down += w_i[i]*a_s**4*c_ix[i]*c_iy[i]
        Os_up += w_i[i]*(1 + 
                         a_s**2*u_max*c_ix[i] + 
                         a_s**4/2*u_max**2*(c_ix[i]**2 - 1/a_s**2))
        
    m_xy = (Is_up - m_xy_I*Os_up)/(m_xy_I*(1-omega)*Os_down - Is_down)
    return m_xy

def m_xy_wall(m_xy_I, 
              c_ix, 
              c_iy, 
              I_s, 
              O_s):
    Is_up = 0
    Is_down = 0
    Os_down = 0
    Os_up = 0
    for i in I_s:
        Is_up += w_i[i]*c_ix[i]*c_iy[i]
        Is_down += w_i[i]*a_s**4*c_ix[i]**2*c_iy[i]**2
    for i in O_s:
        Os_up += w_i[i]
        Os_down += w_i[i]*a_s**4*c_ix[i]*c_iy[i]
    m_xy = (Is_up - m_xy_I*Os_up)/(m_xy_I*(1-omega)*Os_down - Is_down)
    return m_xy
    
def rho_north(m_xy, 
              f_i, 
              ux, 
              I_s, 
              O_s,
              g_id):
    rho_I_rho = 0
    rho_I = 0

    for i in O_s:
        rho_I_rho += (w_i[i]*(1 + 
                             a_s**2*ux*c_ix[i] + 
                             a_s**4/2*ux**2*(c_ix[i]**2 - 1/a_s**2)) + 
                             w_i[i]*(1-omega)*a_s**4*m_xy*c_ix[i]*c_iy[i])
        
    for i in I_s:
        p_id = pop_id(g_id, i)
        rho_I += f_i[p_id]
    
    rho = rho_I/rho_I_rho
    return rho

def rho_wall(m_xy, 
             f_i, 
             I_s, 
             O_s, 
             g_id):
    rho_I_rho = 0
    rho_I = 0
    for i in O_s:
        rho_I_rho += (1-omega)*w_i[i]*a_s**4*m_xy*c_ix[i]*c_iy[i] + w_i[i]
    
    for i in I_s:
        p_id = pop_id(g_id, i)
        rho_I += f_i[p_id]
    
    rho = rho_I/rho_I_rho
    return rho

def rho_corner(f_i, I_s, O_s, g_id):
    rho_I = 0
    for i in I_s:
        p_id = pop_id(g_id, i)
        rho_I += f_i[p_id]
    
    sum_wi = 0
    for i in O_s:
        sum_wi += w_i[i]

    return rho_I/sum_wi

def calc_rho(f_i, g_id):
    rho = 0
    for i in range(Q):
        p_id = pop_id(g_id, i)
        rho += f_i[p_id]

    return rho

def calc_velocity(f_i, rho, g_id):
    u_x = 0
    u_y = 0
    for i in range(Q):
        p_id = pop_id(g_id, i)
        u_x += f_i[p_id]*c_ix[i]
        u_y += f_i[p_id]*c_iy[i]
    u_x *= 1/rho
    u_y *= 1/rho

    return u_x, u_y

def calc_momentum(f_i, rho, g_id):
    m_xx = 0
    m_yy = 0
    m_xy = 0
    for i in range(Q):
        p_id = pop_id(g_id, i)
        m_xx += f_i[p_id]*(c_ix**2 - 1/(a_s**2))
        m_yy += f_i[p_id]*(c_iy**2 - 1/(a_s**2))
        m_xy += f_i[p_id]*(c_ix*c_iy)
    m_xx *= 1/rho
    m_yy *= 1/rho
    m_xy *= 1/rho

    return m_xx, m_yy, m_xy
