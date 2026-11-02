 
#%%
# import packages

import numpy as np
from numpy import newaxis as na


def calculate_dxdy_eq(dx_all,dy_all,k,a,c):

    # wake filter
    fil_wake = ((dy_all<=(0.5+k*dx_all))&(dy_all>=(-0.5-k*dx_all)))&(dx_all<np.inf)
    dx_all[~fil_wake] = np.inf
    dy_all[~fil_wake] = np.inf

    # calculate dx equivalent
    dx_eq = np.min(dx_all,axis=(-1))

    # distinguish between positive and negative dy
    fil_dy_pos = dy_all>=0
    fil_dy_neg = ~fil_dy_pos
    dy_all_pos = dy_all.copy()
    dy_all_neg = dy_all.copy()
    dy_all_pos[~fil_dy_pos] = np.inf
    dy_all_neg[~fil_dy_neg] = -np.inf

    # calculate dy equivalent for positive and negative dy
    delta_dx_eq = np.zeros_like(dx_all)
    delta_dx_eq[fil_wake] = dx_all[fil_wake]-np.tile(dx_eq[:,na],(1,dx_all.shape[1]))[fil_wake]
    dy_eq_pos = np.min(dy_all_pos+a*delta_dx_eq,axis=(-1))
    dy_eq_neg = -np.min(np.abs(dy_all_neg)+a*delta_dx_eq,axis=(-1))

    # calculate dy equivalent
    dy_eq = np.zeros_like(dx_eq)

    fil_dominant_side = np.zeros_like(dy_eq,dtype=bool)
    fil_wake_eq = np.any(fil_wake,axis=(1))

    fil_dominant_side[fil_wake_eq] = np.abs(dy_eq_pos[fil_wake_eq]-np.abs(dy_eq_neg[fil_wake_eq]))>c
    dy_eq[fil_dominant_side] = np.sign(-(dy_eq_pos[fil_dominant_side]-np.abs(dy_eq_neg[fil_dominant_side])))*np.minimum(dy_eq_pos[fil_dominant_side],np.abs(dy_eq_neg)[fil_dominant_side])

    return dx_eq,dy_eq




#%%
# weight functions for dx

def weight_dx_linear(dy,m,a):
    return np.maximum(np.ones_like(dy),1+m*np.abs(dy)-m*a)

def weight_dx_exp(dy,m,a):
    return np.exp(m*np.abs(dy)+np.exp(m*a)+1)

#%%
# weight functions for dy

def weight_dy_linear(dx,m):
    eps = 1e-5
    return np.maximum(np.ones_like(dx)*eps,1-m*dx)

def weight_dy_exp(dx,m):
    return np.exp(-m*dx)


#%%
# dx dy equivalent functions (calculation along the last dimension)

def calculate_dx_equivalent(dx_all,
                            dy_all,
                            weighting,
                            m_eq_x,
                            a_eq_x):

    if weighting == 'wake_filter':
        eps = 1e-8
        w = np.ones_like(dy_all)
        fil_wake = np.abs(dy_all)>(0.5+a_eq_x*dx_all)
        w[fil_wake] = 1/eps
    elif weighting == 'linear':
        w = weight_dx_linear(dy_all,m_eq_x,a_eq_x)
    elif weighting == 'exponential':
        w = weight_dx_exp(dy_all,m_eq_x,a_eq_x)
    else:
        raise TypeError('Weighting fucntion not existent.')  
    
    ind_dx_eq = np.argmin(dx_all*w,axis=-1)
    dx_eq = dx_all[np.arange(dx_all.shape[0]),ind_dx_eq]
    return dx_eq


def calculate_dy_equivalent(dx_all,
                            dy_all,
                            weighting,
                            m_eq_y):

    if weighting == 'constant':
        w = np.ones_like(dx_all)
    elif weighting == 'linear':
        w = weight_dy_linear(dx_all,m_eq_y)
    elif weighting == 'exponential':
        w = weight_dy_exp(dx_all,m_eq_y)
    else:
        raise TypeError('Weighting fucntion not existent.')  
    
    dy_eq = np.sum(dy_all*w,axis=(-1))/np.sum(w,axis=(-1))
    return dy_eq