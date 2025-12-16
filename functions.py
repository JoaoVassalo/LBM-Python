from constants import *

def grid_id(x, y, Nx):
    return (x + Nx*y)

def pop_id(g_id, i, Q):
    return (g_id*Q + i)

def equilibrium(rho, w_i, a_s, ux, c_ix, uy, c_iy):
    return rho*w_i*( 1 + a_s**2*( ux*c_ix + uy*c_iy ) 
                    + a_s**4/2*( ( ux*c_ix + uy*c_iy )**2 
                                - ( ux**2 + uy**2 )/a_s**2 ) )

def m_xy_I(f_i, I_s, c_ix, c_iy, Q, g_id):
    num = 0
    div = 0

    for i in I_s:
        p_id = pop_id(g_id, i, Q)
        num += f_i[p_id] * c_ix[i]*c_iy[i]
        div += f_i[p_id]

    return (num/div)

def m_xy_north(m_xy_I, 
               c_ix, 
               c_iy, 
               u_max, 
               I_s, 
               O_s, 
               a_s, 
               w_i, 
               omega):
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
              O_s, 
              a_s, 
              w_i, 
              omega):
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
              c_ix, 
              c_iy, 
              ux, 
              I_s, 
              O_s, 
              a_s, 
              w_i, 
              omega,
              Q,
              g_id):
    rho_I_rho = 0
    rho_I = 0

    for i in O_s:
        rho_I_rho += (w_i[i]*(1 + 
                             a_s**2*ux*c_ix[i] + 
                             a_s**4/2*ux**2*(c_ix[i]**2 - 1/a_s**2)) + 
                             w_i[i]*(1-omega)*a_s**4*m_xy*c_ix[i]*c_iy[i])
        
    for i in I_s:
        p_id = pop_id(g_id, i, Q)
        rho_I += f_i[p_id]
    
    rho = rho_I/rho_I_rho
    return rho

def rho_wall(m_xy, 
             f_i, 
             c_ix, 
             c_iy, 
             I_s, 
             O_s, 
             a_s, 
             w_i, 
             omega,
             Q,
             g_id):
    rho_I_rho = 0
    rho_I = 0
    for i in O_s:
        rho_I_rho += (1-omega)*w_i[i]*a_s**4*m_xy*c_ix[i]*c_iy[i] + w_i[i]
    
    for i in I_s:
        p_id = pop_id(g_id, i, Q)
        rho_I += f_i[p_id]
    
    rho = rho_I/rho_I_rho
    return rho

def rho_corner(f_i, I_s, O_s, Q, g_id):
    rho_I = 0
    for i in I_s:
        p_id = pop_id(g_id, i, Q)
        rho_I += f_i[p_id]
    
    sum_wi = 0
    for i in O_s:
        sum_wi += w_i[i]

    return rho_I/sum_wi

def calc_rho(f_i, Q, g_id):
    rho = 0
    for i in range(Q):
        p_id = pop_id(g_id, i, Q)
        rho += f_i[p_id]

    return rho

def calc_velocity(f_i, c_ix, c_iy, Q, g_id, rho):
    u_x = 0
    u_y = 0
    for i in range(Q):
        p_id = pop_id(g_id, i, Q)
        u_x += f_i[p_id]*c_ix[i]
        u_y += f_i[p_id]*c_iy[i]
    u_x *= 1/rho
    u_y *= 1/rho

    return u_x, u_y

def calc_momentum(f_i, c_ix, c_iy, Q, g_id, rho, a_s):
    m_xx = 0
    m_yy = 0
    m_xy = 0
    for i in range(Q):
        p_id = pop_id(g_id, i, Q)
        m_xx += f_i[p_id]*(c_ix**2 - 1/(a_s**2))
        m_yy += f_i[p_id]*(c_iy**2 - 1/(a_s**2))
        m_xy += f_i[p_id]*(c_ix*c_iy)
    m_xx *= 1/rho
    m_yy *= 1/rho
    m_xy *= 1/rho

    return m_xx, m_yy, m_xy
