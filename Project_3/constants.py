

class constants():
    def __init__(self):
        self.G = 6.6742e-11                             # N m^2 kg^-2
        self.M_sun = 1.989e30                           # kg
        self.R_sun = 6.9634e8                             # m
        self.T_sun = 5778                               # K (surface temp)
        self.rho_sun = 2e-4                             # kg m^-3 (surface density)
        self.P_sun = 1.8e8                              # Pa (surface pressure)
        self.gamma = 1 + 2/3                            # cp/cv (assuming dof=3)
        self.m_u = 1.6605e-27                           # kg
        self.k_b = 1.3806e-23                           # m^2 kg s^-2 K^-1

        self.g = self.G * self.M_sun / self.R_sun**2    # m s^-2
