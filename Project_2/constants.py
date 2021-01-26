import numpy as np

class constants:
    def __init__(self,T,rho):
        # T and rho is the only inputs
        self.T_sol = T
        self.rho_sol = rho

        # pre defined constants
        self.m_u = 1.6605e-27                           # kg
        self.T_9 = self.T_sol*1e-9                      # 10^9 K
        self.N_A = 6.0221e23                            # avogadros constant

        # mass fraction of different elements
        self.X = 0.7
        self.Y_3He = 1e-10
        self.Y = 0.29
        self.Z_7Li = 1e-7
        self.Z_7Be = 1e-7
        self.Z_14N = 1e-11

        # electron density
        #self.n_etot = (self.rho_sol/self.m_u) * (self.X+(2/3)*self.Y_3He
        #                +self.Y/2+(3/7)*self.Z_7Li+(4/7)*self.Z_7Be+self.Z_14N/2)
        self.n_etot = self.rho_sol*(1+self.X)/(2*self.m_u)

        # lambda values needed to find the reaction rates
        self.T_9str1 = self.T_9/(1+4.95e-2*self.T_9)
        self.T_9str2 = self.T_9/(1+0.759*self.T_9)
        self.l_pp = (4.01e-15*self.T_9**(-2/3)*np.exp(-3.380*self.T_9**(-1/3))*(1+0.123*self.T_9**(1/3)+1.09*self.T_9**(2/3)+0.938*self.T_9))/self.N_A/1e6
        self.l_33 = (6.04e10*self.T_9**(-2/3)*np.exp(-12.276*self.T_9**(-1/3))*(1+0.034*self.T_9**(1/3)-0.522*self.T_9**(2/3)-0.124*self.T_9+0.353*self.T_9**(4/3)+0.213*self.T_9**(5/3)))/self.N_A/1e6
        self.l_34 = (5.61e6*self.T_9str1**(5/6)*self.T_9**(-3/2)*np.exp(-12.826*self.T_9str1**(-1/3)))/self.N_A/1e6
        self.l_e7 = (1.34e-10*self.T_9**(-1/2)*(1-0.537*self.T_9**(1/3)+3.86*self.T_9**(2/3)+0.0027*self.T_9**(-1)*np.exp(2.515e-3*self.T_9**(-1))))/self.N_A/1e6
        # electron capture limit
        if self.T_sol<1e6 and self.l_e7>1.57e-7/(self.n_etot*self.N_A): self.l_e7=1.57e-7/(self.n_etot*self.N_A*1e6)
        self.l_17mrk = (1.096e9*self.T_9**(-2/3)*np.exp(-8.472*self.T_9**(-1/3))-4.83e8*self.T_9str2**(5/6)*self.T_9**(-3/2)*np.exp(-8.472*self.T_9str2**(-1/3))+1.06e10*self.T_9**(-3/2)*np.exp(-30.442*self.T_9**(-1)))/self.N_A/1e6
        self.l_17 = (3.11e5*self.T_9**(-2/3)*np.exp(-10.262*self.T_9**(-1/3))+2.53e3*self.T_9**(-3/2)*np.exp(-7.306*self.T_9**(-1)))/self.N_A/1e6
        self.l_p14 = (4.9e7*self.T_9**(-2/3)*np.exp(-15.228*self.T_9**(-1/3)-0.092*self.T_9**2)*(1+0.027*self.T_9**(1/3)-0.778*self.T_9**(2/3)-0.149*self.T_9+0.261*self.T_9**(4/3)+0.127*self.T_9**(5/3))+2.37e3*self.T_9**(-3/2)*np.exp(-3.011*self.T_9**(-1))+2.19e4*np.exp(-12.53*self.T_9**(-1)))/self.N_A/1e6

        # realesed energy to the thermal bath
        self.Q_pp = 1.177*1.6022e-13
        self.Q_pd = 5.494*1.6022e-13
        self.Q_33 = 12.86*1.6022e-13
        self.Q_34 = 1.586*1.6022e-13
        self.Q_e7 = 0.049*1.6022e-13
        self.Q_17mrk = 17.346*1.6022e-13
        self.Q_17 = 0.137*1.6022e-13
        self.Q_p14 = 7.297*1.6022e-13
        self.Q_8 = 8.367*1.6022e-13
        self.Q_8mrk = 2.995*1.6022e-13
        self.Q_p12 = 1.944*1.6022e-13
        self.Q_13 = 1.513*1.6022e-13
        self.Q_p13 = 7.551*1.6022e-13
        self.Q_p14 = 7.297*1.6022e-13
        self.Q_15 = 1.757*1.6022e-13
        self.Q_p15 = 4.966*1.6022e-13

        # realesed neutrino energy in reactions
        self.Q_nu_pp = 0.265*1.6022e-13
        self.Q_nu_e7 = 0.815*1.6022e-13
        self.Q_nu_8 = 6.711*1.6022e-13
        self.Q_nu_13 = 0.707*1.6022e-13
        self.Q_nu_15 = 0.997*1.6022e-13
