
#%%
# import packages
import numpy as np
from numpy import newaxis as na


# %%
# old functions

def calculate_geometric_yaw(dx,dy,yaw_max,m_x,m_y,m_center_x,c_center_y,k):

    yaw_geom = np.sign(dy)*np.maximum((yaw_max-m_x*dx-m_y*np.abs(dy)),np.zeros_like(dx))

    fil_wake = (dy<=(0.5+k*dx))&(dy>=(-0.5-k*dx))
    yaw_geom[~fil_wake] = 0

    fil_ub = dy>(c_center_y+(m_center_x*dx)**2)
    fil_lb = dy<(-c_center_y-(m_center_x*dx)**2)
    yaw_geom[(~fil_ub)&(~fil_lb)] = 0

    return yaw_geom


#def calculate_geometric_helix(dx,dy,helix_amp_max,c_x,c_y,m_center_x,c_center_y,k):
def calculate_geometric_helix(dx,dy,helix_amp_max,k_x,x_0,m_center_x,c_center_y,k):

    #helix_amp_geom = helix_amp_max*np.ones_like(dx)
    #helix_amp_geom[dx>c_x] = 0
    #helix_amp_geom[np.abs(dy)>c_y] = 0

    helix_amp_geom = helix_amp_max*(1/(1+np.exp(k_x*(dx-x_0))))

    fil_wake = (dy<=(0.5+k*dx))&(dy>=(-0.5-k*dx))
    helix_amp_geom[~fil_wake] = 0

    fil_ub = dy>(c_center_y+(m_center_x*dx)**2)
    fil_lb = dy<(-c_center_y-(m_center_x*dx)**2)
    helix_amp_geom[(fil_ub)|(fil_lb)] = 0

    return helix_amp_geom


def calculate_geometric_switch(dx,dy,m_x,c_y):

    fil_ub = dy>(c_y+(m_x*dx)**2)
    fil_lb = dy<(-c_y-(m_x*dx)**2)

    fil_yaw = (fil_ub)|(fil_lb)
    fil_helix = (~fil_ub)&(~fil_lb)

    return fil_yaw,fil_helix



def calculate_dxdy(x,y,wd,diameter):
    
    # extend dimensions of x and y
    x_mat_1 = np.tile(np.reshape(x,(len(x),1)),(1,len(x)))
    x_mat_2 = np.tile(np.reshape(x,(1,len(x))),(len(x),1))
    y_mat_1 = np.tile(np.reshape(y,(len(y),1)),(1,len(y)))
    y_mat_2 = np.tile(np.reshape(y,(1,len(y))),(len(y),1))

    # calculate dx and dy
    d = np.sqrt((x_mat_1-x_mat_2)**2+(y_mat_1-y_mat_2)**2)
    theta = np.arctan2(y_mat_2-y_mat_1,x_mat_2-x_mat_1)
    gamma = wd*(np.pi/180)-(3/2)*np.pi+theta
    dx = d*np.cos(gamma)/diameter
    dy = d*np.sin(gamma)/diameter

    # assign inf in case of negative dx
    eps = 1e-3
    dx[dx<eps] = np.inf
    
    return dx,dy



def  ws_eff_filter(dx_all,dy_all,ws_eff_i,ws_rated):
    ws_eff_ii = np.tile(ws_eff_i[:,na],(1,len(ws_eff_i)))
    fil_ws_eff = ws_eff_ii>=ws_rated
    dx_all_fil = dx_all.copy()
    dy_all_fil = dy_all.copy()
    dx_all_fil[fil_ws_eff] = np.inf
    return dx_all_fil,dy_all_fil


#%%

coeff_dict_IEA22MW_EmpGauss = {
    'k_yaw': 0.01111,
    'yaw_max': 22.67,
    'yaw_m_x': 0.2651,
    'yaw_m_y': 15.00,
    'yaw_m_center_x': 1.000e-06,
    'yaw_c_center_y': 0.000,
    # 'k_helix': 0.02000,
    # 'helix_amp_max': 5.000,
    # 'helix_k_x': 1.000,
    # 'helix_x_0': 30.00,
    # 'helix_m_center_x': 0.01500,
    # 'helix_c_center_y': 0.3000,
    # 'sw_m_x': 0.01200,
    # 'sw_c_y': 0.02000,
    # 'k_farmEqFunc' : 0.05000,
    # 'a_farmEqFunc' : 0.000,
    # 'c_farmEqFunc' : 1.000,
}

class GeomWFFC():

    def __init__(self,coeff_dict=coeff_dict_IEA22MW_EmpGauss):

        self.coeff_dict = coeff_dict

        # yaw coefficients
        self.k_yaw = coeff_dict.get('k_yaw',coeff_dict_IEA22MW_EmpGauss['k_yaw'])
        self.yaw_max = coeff_dict.get('yaw_max',coeff_dict_IEA22MW_EmpGauss['yaw_max'])
        self.yaw_m_x = coeff_dict.get('yaw_m_x',coeff_dict_IEA22MW_EmpGauss['yaw_m_x'])
        self.yaw_m_y = coeff_dict.get('yaw_m_y',coeff_dict_IEA22MW_EmpGauss['yaw_m_y'])
        self.yaw_m_center_x = coeff_dict.get('yaw_m_center_x',coeff_dict_IEA22MW_EmpGauss['yaw_m_center_x'])
        self.yaw_c_center_y = coeff_dict.get('yaw_c_center_y',coeff_dict_IEA22MW_EmpGauss['yaw_c_center_y'])
        
        # # helix coefficients
        # self.k_helix = coeff_dict.get('k_helix',coeff_dict_IEA22MW_EmpGauss['k_helix'])
        # self.helix_amp_max = coeff_dict.get('helix_amp_max',coeff_dict_IEA22MW_EmpGauss['helix_amp_max'])
        # self.helix_k_x = coeff_dict.get('helix_k_x',coeff_dict_IEA22MW_EmpGauss['helix_k_x'])
        # self.helix_x_0 = coeff_dict.get('helix_x_0',coeff_dict_IEA22MW_EmpGauss['helix_x_0'])
        # self.helix_m_center_x = coeff_dict.get('helix_m_center_x',coeff_dict_IEA22MW_EmpGauss['helix_m_center_x'])
        # self.helix_c_center_y = coeff_dict.get('helix_c_center_y',coeff_dict_IEA22MW_EmpGauss['helix_c_center_y'])
        
        # # combined coefficients
        # self.sw_m_x = coeff_dict.get('sw_m_x',coeff_dict_IEA22MW_EmpGauss['sw_m_x'])
        # self.sw_c_y = coeff_dict.get('sw_c_y',coeff_dict_IEA22MW_EmpGauss['sw_c_y'])

        # farm equivalent function coefficients
        self.k_farmEqFunc = coeff_dict.get('k_farmEqFunc',coeff_dict_IEA22MW_EmpGauss['k_farmEqFunc'])
        self.a_farmEqFunc = coeff_dict.get('a_farmEqFunc',coeff_dict_IEA22MW_EmpGauss['a_farmEqFunc'])
        self.c_farmEqFunc = coeff_dict.get('c_farmEqFunc',coeff_dict_IEA22MW_EmpGauss['c_farmEqFunc'])


    def _calculate_geometric_yaw(self,dx,dy):
        yaw_geom = np.sign(dy)*np.maximum((self.yaw_max-self.yaw_m_x*dx-self.yaw_m_y*np.abs(dy)),np.zeros_like(dx))
        fil_wake = (dy<=(0.5+self.k_yaw*dx))&(dy>=(-0.5-self.k_yaw*dx))
        yaw_geom[~fil_wake] = 0
        fil_ub = dy>(self.yaw_c_center_y+(self.yaw_m_center_x*dx)**2)
        fil_lb = dy<(-self.yaw_c_center_y-(self.yaw_m_center_x*dx)**2)
        yaw_geom[(~fil_ub)&(~fil_lb)] = 0
        return yaw_geom


    # def _calculate_geometric_helix(self,dx,dy):
    #     helix_amp_geom = self.helix_amp_max*(1/(1+np.exp(self.helix_k_x*(dx-self.helix_x_0))))
    #     fil_wake = (dy<=(0.5+self.k_helix*dx))&(dy>=(-0.5-self.k_helix*dx))
    #     helix_amp_geom[~fil_wake] = 0
    #     fil_ub = dy>(self.helix_c_center_y+(self.helix_m_center_x*dx)**2)
    #     fil_lb = dy<(-self.helix_c_center_y-(self.helix_m_center_x*dx)**2)
    #     helix_amp_geom[(fil_ub)|(fil_lb)] = 0
    #     return helix_amp_geom


    # def _calculate_geometric_switch(self,dx,dy):
    #     fil_ub = dy>(self.sw_c_y+(self.sw_m_x*dx)**2)
    #     fil_lb = dy<(-self.sw_c_y-(self.sw_m_x*dx)**2)
    #     fil_yaw = (fil_ub)|(fil_lb)
    #     fil_helix = (~fil_ub)&(~fil_lb)
    #     return fil_yaw,fil_helix



    def _calculate_dxdy(self,x,y,wd,diameter):
        
        # extend dimensions of x and y
        x_mat_1 = np.tile(np.reshape(x,(len(x),1)),(1,len(x)))
        x_mat_2 = np.tile(np.reshape(x,(1,len(x))),(len(x),1))
        y_mat_1 = np.tile(np.reshape(y,(len(y),1)),(1,len(y)))
        y_mat_2 = np.tile(np.reshape(y,(1,len(y))),(len(y),1))

        # calculate dx and dy
        d = np.sqrt((x_mat_1-x_mat_2)**2+(y_mat_1-y_mat_2)**2)
        theta = np.arctan2(y_mat_2-y_mat_1,x_mat_2-x_mat_1)
        gamma = wd*(np.pi/180)-(3/2)*np.pi+theta
        dx = d*np.cos(gamma)/diameter
        dy = d*np.sin(gamma)/diameter

        # assign inf in case of negative dx
        eps = 1e-3
        dx[dx<eps] = np.inf
        
        return dx,dy



    def _ws_eff_filter(self,dx_ii,dy_ii,ws_eff_ik,ws_rated):
        ws_eff_iik = np.tile(ws_eff_ik[na,:,:],(ws_eff_ik.shape[0],1,1))
        fil_ws_eff = ws_eff_iik>=ws_rated
        dx_iik = np.tile(dx_ii[:,:,na],(1,1,ws_eff_ik.shape[1]))
        dy_iik = np.tile(dy_ii[:,:,na],(1,1,ws_eff_ik.shape[1]))
        dx_iik[fil_ws_eff] = np.inf
        return dx_iik,dy_iik
    

    def _calculate_dxdy_eq(self,
                           x,
                           y,
                           diameter,
                           ws_rated,
                           wd_array,
                           ws_array,
                           ws_eff_ilk):
        
        # initialize output
        dx_eq_ilk = np.zeros((len(x),len(wd_array),len(ws_array)))
        dy_eq_ilk = np.zeros((len(x),len(wd_array),len(ws_array)))
    
        # iterate for each wind direction
        for wd_ind in np.arange(len(wd_array)):

            # apply ws filter
            wd = wd_array[wd_ind]
            dx_ii,dy_ii = self._calculate_dxdy(x,y,wd,diameter)
            dx_iik,dy_iik = self._ws_eff_filter(dx_ii,dy_ii,ws_eff_ilk[:,wd_ind,:],ws_rated)

            # wake filter
            fil_wake = ((dy_iik<=(0.5+self.k_farmEqFunc*dx_iik))&(dy_iik>=(-0.5-self.k_farmEqFunc*dx_iik)))&(dx_iik<np.inf)
            dx_iik[~fil_wake] = np.inf
            dy_iik[~fil_wake] = np.inf

            # calculate dx equivalent
            dx_eq_ik = np.min(dx_iik,axis=(1))

            # distinguish between positive and negative dy
            fil_dy_pos = dy_iik>=0
            fil_dy_neg = ~fil_dy_pos
            dy_pos_iik = dy_iik.copy()
            dy_neg_iik = dy_iik.copy()
            dy_pos_iik[~fil_dy_pos] = np.inf
            dy_neg_iik[~fil_dy_neg] = -np.inf

            # calculate dy equivalent for positive and negative dy
            delta_dx_eq_iik = np.zeros_like(dx_iik)
            delta_dx_eq_iik[fil_wake] = dx_iik[fil_wake]-np.tile(dx_eq_ik[:,na,:],(1,dx_iik.shape[1],1))[fil_wake]
            dy_eq_pos_ik = np.min(dy_pos_iik+self.a_farmEqFunc*delta_dx_eq_iik,axis=(1))
            dy_eq_neg_ik = -np.min(np.abs(dy_neg_iik)+self.a_farmEqFunc*delta_dx_eq_iik,axis=(1))

            # calculate dy equivalent
            dy_eq_ik = np.zeros_like(dx_eq_ik)
            fil_dominant_side = np.zeros_like(dy_eq_ik,dtype=bool)
            fil_wake_eq = np.any(fil_wake,axis=(1))
            fil_dominant_side[fil_wake_eq] = np.abs(dy_eq_pos_ik[fil_wake_eq]-np.abs(dy_eq_neg_ik[fil_wake_eq]))>self.c_farmEqFunc
            dy_eq_ik[fil_dominant_side] = np.sign(-(dy_eq_pos_ik[fil_dominant_side]-np.abs(dy_eq_neg_ik[fil_dominant_side])))*np.minimum(dy_eq_pos_ik[fil_dominant_side],np.abs(dy_eq_neg_ik)[fil_dominant_side])

            # assign values
            dx_eq_ilk[:,wd_ind,:] = dx_eq_ik
            dy_eq_ilk[:,wd_ind,:] = dy_eq_ik

        return dx_eq_ilk,dy_eq_ilk



    def calculate_geomYaw(self,
                          x,
                          y,
                          diameter,
                          ws_rated,
                          wd_array,
                          ws_array,
                          ws_eff_ilk):
        
        dx_eq_ilk,dy_eq_ilk = self._calculate_dxdy_eq(x,y,diameter,ws_rated,wd_array,ws_array,ws_eff_ilk)
        yaw_geom_ilk = self._calculate_geometric_yaw(dx_eq_ilk,dy_eq_ilk)
        return yaw_geom_ilk


    # def calculate_geomHelix(self,
    #                         x,
    #                         y,
    #                         diameter,
    #                         ws_rated,
    #                         wd_array,
    #                         ws_array,
    #                         ws_eff_ilk):
        
    #     dx_eq_ilk,dy_eq_ilk = self._calculate_dxdy_eq(x,y,diameter,ws_rated,wd_array,ws_array,ws_eff_ilk)
    #     helix_amp_geom_ilk = self._calculate_geometric_helix(dx_eq_ilk,dy_eq_ilk)
    #     return helix_amp_geom_ilk


    # def calculate_geomCombined(self,
    #                            x,
    #                            y,
    #                            diameter,
    #                            ws_rated,
    #                            wd_array,
    #                            ws_array,
    #                            ws_eff_ilk):
        
    #     dx_eq_ilk,dy_eq_ilk = self._calculate_dxdy_eq(x,y,diameter,ws_rated,wd_array,ws_array,ws_eff_ilk)
    #     yaw_geom_ilk = self._calculate_geometric_yaw(dx_eq_ilk,dy_eq_ilk)
    #     helix_amp_geom_ilk = self._calculate_geometric_helix(dx_eq_ilk,dy_eq_ilk)
    #     fil_yaw_geom,fil_helix_geom = self._calculate_geometric_switch(dx_eq_ilk,dy_eq_ilk)
    #     yaw_geom_ilk[~fil_yaw_geom] = 0 
    #     helix_amp_geom_ilk[~fil_helix_geom] = 0 
    #     return yaw_geom_ilk,helix_amp_geom_ilk



# %%
