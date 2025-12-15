def grid_id(x, y, Nx):
    return (x + Nx*y)

def pop_id(x, y, i, Nx, Q):
    return ((x+y*Nx)*Q + i)

def equilibrium(rho, w_i, a_s, ux, c_ix, uy, c_iy):
    return rho*w_i*( 1 + a_s**2*( ux*c_ix + uy*c_iy ) 
                    + a_s**4/2*( ( ux*c_ix + uy*c_iy )**2 
                                - ( ux**2 + uy**2 )/a_s**2 ) )

def m_xy_I(f_i, I_s, c_ix, c_iy, Q):
    num = 0
    div = 0
    for i in I_s:
        num += f_i * c_ix[i]*c_iy[i]
    for i in range(Q):
        div += f_i
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
                         a_s**2*u_max*c_ix + 
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
        Is_up += w_i[i]*c_ix*c_iy
        Is_down += w_i[i]*a_s**4*c_ix**2*c_iy**2
    for i in O_s:
        Os_up += w_i[i]
        Os_down += w_i[i]*a_s**4*c_ix*c_iy
    m_xy = (Is_up - m_xy_I*Os_up)/(m_xy_I*(1-omega)*Os_down - Is_down)
    return m_xy
    
def rho_north(m_xy, 
              f_i, 
              c_ix, 
              c_iy, 
              u_max, 
              I_s, 
              O_s, 
              a_s, 
              w_i, 
              omega):
    rho_I_rho = 0
    rho_I = 0
    for i in O_s:
        rho_I_rho += (w_i[i]*(1 + 
                             a_s**2*u_max*c_ix + 
                             a_s**4/2*u_max**2*(c_ix[i]**2 - 1/a_s**2)) + 
                             w_i[i]*(1-omega)*a_s**4*m_xy*c_ix[i]*c_iy[i])
    for i in I_s:
        rho_I += f_i[i]
    
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
             omega):
    rho_I_rho = 0
    rho_I = 0
    for i in O_s:
        rho_I_rho += (1-omega)*w_i*a_s**4*m_xy*c_ix*c_iy + w_i
    
    for i in I_s:
        rho_I += f_i[i]
    
    rho = rho_I/rho_I_rho
    return rho

def rho_corner(f_i, I_s):
    rho = 0
    for i in I_s:
        rho += f_i[i]
    return rho

