# visulaliser
import FVis3 as FVis
import matplotlib.pyplot as plt
import constants
import numpy as np
from scipy import signal


class convection2D:

    def __init__(self):

        """
        define variables
        """

        self.X = 12e6               # m
        self.Y = 4e6                # m
        self.nx = 300
        self.ny = 100
        self.dx = self.X/self.nx
        self.dy = self.Y/self.ny
        self.dt = 1e-16

        self.tot_t = 0              # total time of convection
        tot_time_steps = 2          # next and previous timestep

        self.T = np.zeros((self.nx, self.ny))
        self.P = np.zeros((self.nx, self.ny))
        self.rho = np.zeros((tot_time_steps, self.nx, self.ny))
        self.e = np.zeros((tot_time_steps, self.nx, self.ny))

        self.u = np.zeros((tot_time_steps, self.nx, self.ny))
        self.w = np.zeros((tot_time_steps, self.nx, self.ny))

        self.C = constants.constants()
        self.mu = 0.61
        self.nab = 2/5 + 0.00001             # double logarithmic gradient



    def initialise(self, perturbation=False):

        """
        initialise temperature, pressure, density and internal energy
        """
        C = self.C

        self.P[:,0] = C.P_sun
        self.T[:,0] = C.T_sun
        self.e[0,:,0] = self.P[:,0] / (C.gamma-1)
        self.rho[0,:,0] = self.P[:,0]*self.mu*C.m_u / (C.k_b*self.T[:,0])



        # initialize in the y-direction
        for i in range(self.ny-1):
            self.P[:,i+1] = self.P[:,i] - self.rho[0,:,i]*(-C.g)*self.dy
            H_P = self.P[:,i] / (self.rho[0,:,i]*(-C.g))
            self.T[:,i+1] = self.T[:,i] - self.nab*self.T[:,i]*self.dy / H_P
            self.e[0,:,i+1] = self.P[:,i+1] / (C.gamma-1)
            self.rho[0,:,i+1] = self.P[:,i+1]*self.mu*C.m_u / (C.k_b*self.T[:,i+1])

        # temperature Gaussian perturbation
        if perturbation==True:
            mat = np.zeros((self.nx, self.ny))
            zx = signal.gaussian(self.nx, std=10)
            zy = signal.gaussian(self.ny, std=10)
            #zx = signal.gaussian(self.nx, std=self.nx/25)
            #zy = signal.gaussian(self.ny, std=self.ny/25)


            for i in range(self.nx):
                for j in range(self.ny):
                    mat[i,j] = zx[i]*zy[j]

            self.T = self.T + 0.6*mat*self.T
            self.rho[0,:,:] = self.P*self.mu*C.m_u / (C.k_b*self.T)



        self.T = np.flip(self.T, axis=1)
        self.P = np.flip(self.P, axis=1)
        self.e[0,:,:] = np.flip(self.e[0,:,:], axis=1)
        self.rho[0,:,:] = np.flip(self.rho[0,:,:], axis=1)

        x = int(self.nx/2)
        print("P[%i,0]:     %.2e" % (x, self.P[x,0]))
        print("P[%i,%i]:    %.2e" % (x, self.ny-1, self.P[x,self.ny-1]))
        print("T[%i,0]:     %.2e" % (x, self.T[x,0]))
        print("T[%i,%i]:    %.2e" % (x, self.ny-1, self.T[x,self.ny-1]))
        print("e[%i,0]:     %.2e" % (x, self.e[0,x,0]))
        print("e[%i,%i]:    %.2e" % (x, self.ny-1, self.e[0,x,self.ny-1]))
        print("rho[%i,0]:   %.2f" % (x, self.rho[0,x,0]))
        print("rho[%i,%i]:  %.2f" % (x, self.ny-1, self.rho[0,x,self.ny-1]))
        print("")



    def timestep(self):

        """
        calculate timestep
        """

        p = 0.1            # some small number

        rel_u = np.zeros((self.nx, self.ny))
        rel_w = np.zeros((self.nx, self.ny))
        rel_rhou = np.zeros((self.nx, self.ny))
        rel_rhow = np.zeros((self.nx, self.ny))
        rel_rho = np.zeros((self.nx, self.ny))
        rel_e = np.zeros((self.nx, self.ny))

        rel_rho[:,1:-1] = np.abs(self.drho_dt[:,1:-1] / self.rho[1,:,1:-1])
        rel_e[:,1:-1] = np.abs(self.de_dt[:,1:-1] / self.e[1,:,1:-1])
        rel_x = np.abs(self.u[1,:,1:-1] / self.dx)
        rel_y = np.abs(self.w[1,:,1:-1] / self.dy)

        u0, u1 = np.where(np.abs(np.log10(np.abs(self.u[1,:,1:-1]))+2) < 3)
        w0, w1 = np.where(np.abs(np.log10(np.abs(self.w[1,:,1:-1]))+2) < 3)
        rel_rhou[u0,u1+1] = np.abs(self.drhou_dt[u0,u1+1]/(self.rho[1,u0,u1+1]*self.u[1,u0,u1+1]))
        rel_rhow[w0,w1+1] = np.abs(self.drhow_dt[w0,w1+1]/(self.rho[1,w0,w1+1]*self.w[1,w0,w1+1]))

        rel_u = np.nan_to_num(rel_u, nan=np.inf)
        rel_w = np.nan_to_num(rel_w, nan=np.inf)

        delta = np.max([np.max(rel_rho), np.max(rel_u), np.max(rel_w), np.max(rel_e), np.max(rel_x), np.max(rel_y)])

        self.dt = p/delta
        dt2 = self.dt
        dt_lim_high = 1e-1
        dt_lim_low = 1e-3
        if delta==np.inf:
            self.dt = dt_lim_high
        if self.dt > dt_lim_high:
            self.dt = dt_lim_high
        if self.dt < dt_lim_low:
            self.dt = dt_lim_low
        #print("Delta: %.2e    dt: %.2e    calc dt: %.2e    total time: %.2f" % (delta, self.dt, dt2, self.tot_t))



    def boundary_conditions(self):

        """
        boundary conditions for energy, density and velocity
        """

        C = self.C

        # Vertical boundary: vertical velocity
        self.w[0,:,0] = 0
        self.w[0,:,-1] = 0

        # Vertical boundary: horizontal velocity
        self.u[0,:,0] = (4*self.u[0,:,1] - self.u[0,:,2]) / 3
        self.u[0,:,-1] = (4*self.u[0,:,-2] - self.u[0,:,-3]) / 3


        # Vertical boundary: density and energy
        self.e[0,:,0] = self.P[:,0] / (C.gamma-1)
        self.e[0,:,-1] = self.P[:,-1] / (C.gamma-1)
        self.rho[0,:,0] = self.e[0,:,0]*(C.gamma-1)*self.mu*C.m_u / (C.k_b*self.T[:,0])
        self.rho[0,:,-1] = self.e[0,:,-1]*(C.gamma-1)*self.mu*C.m_u / (C.k_b*self.T[:,-1])



    def central_x(self,func):

        """
        central difference scheme in x-direction
        """

        der_func = np.zeros((self.nx, self.ny))

        forward_func = np.roll(func, -1, axis=0)
        backward_func = np.roll(func, 1, axis=0)

        der_func = (forward_func - backward_func) / (2*self.dx)
        return der_func



    def central_y(self,func):

        """
        central difference scheme in y-direction
        """

        der_func = np.zeros((self.nx, self.ny))

        forward_func = np.roll(func, -1, axis=1)
        backward_func = np.roll(func, 1, axis=1)

        der_func = (forward_func - backward_func) / (2*self.dy)
        return der_func



    def upwind_x(self,func,u):

        """
        upwind difference scheme in x-direction
        """

        der_func = np.zeros((self.nx, self.ny))

        forward_func = np.roll(func, -1, axis=0)
        backward_func = np.roll(func, 1, axis=0)

        pos = np.where(u >= 0)
        neg = np.where(u < 0)

        der_func[pos[0],pos[1]] = (func[pos[0],pos[1]] - backward_func[pos[0],pos[1]]) / self.dx
        der_func[neg[0],neg[1]] = (forward_func[neg[0],neg[1]] - func[neg[0],neg[1]]) / self.dx
        return der_func



    def upwind_y(self,func,u):

        """
        upwind difference scheme in y-direction
        """

        der_func = np.zeros((self.nx, self.ny))

        forward_func = np.roll(func, -1, axis=1)
        backward_func = np.roll(func, 1, axis=1)

        pos = np.where(u >= 0)
        neg = np.where(u < 0)

        der_func[pos[0],pos[1]] = (func[pos[0],pos[1]] - backward_func[pos[0],pos[1]]) / self.dy
        der_func[neg[0],neg[1]] = (forward_func[neg[0],neg[1]] - func[neg[0],neg[1]]) / self.dy
        return der_func



    def hydro_solver(self):

        """
        hydrodynamic equations solver
        """

        C = self.C


        # solving continuity equation
        du_dx_c = self.central_x(self.u[0,:,:])
        dw_dy_c = self.central_y(self.w[0,:,:])
        drho_dx = self.upwind_x(self.rho[0,:,:], self.u[0,:,:])
        drho_dy = self.upwind_y(self.rho[0,:,:], self.w[0,:,:])
        self.drho_dt = - self.rho[0,:,:]*(du_dx_c + dw_dy_c) \
                       - self.u[0,:,:]*drho_dx - self.w[0,:,:]*drho_dy
        self.drho_dt[:,0] = 0
        self.drho_dt[:,-1] = 0
        self.rho[1,:,:] = self.rho[0,:,:] + self.drho_dt*self.dt


        # solving the horizontal momentum equation
        du_dx_hm = self.upwind_x(self.u[0,:,:], self.u[0,:,:])
        dw_dy_hm = self.upwind_y(self.w[0,:,:], self.u[0,:,:])
        drhou_dx = self.upwind_x(self.rho[0,:,:]*self.u[0,:,:], self.u[0,:,:])
        drhou_dy = self.upwind_y(self.rho[0,:,:]*self.u[0,:,:], self.w[0,:,:])
        dP_dx = self.central_x(self.P)
        rho_u = self.rho[0,:,:]*self.u[0,:,:]
        self.drhou_dt = - rho_u*(du_dx_hm+dw_dy_hm) - self.u[0,:,:]*drhou_dx \
                        - self.w[0,:,:]*drhou_dy - dP_dx
        self.drhou_dt[:,0] = 0
        self.drhou_dt[:,-1] = 0
        self.u[1,:,:] = (rho_u + self.drhou_dt*self.dt) / self.rho[1,:,:]


        # solving the vertical momentum equation
        du_dx_vm = self.upwind_x(self.u[0,:,:], self.w[0,:,:])
        dw_dy_vm = self.upwind_y(self.w[0,:,:], self.w[0,:,:])
        drhow_dx = self.upwind_x(self.rho[0,:,:]*self.w[0,:,:], self.u[0,:,:])
        drhow_dy = self.upwind_y(self.rho[0,:,:]*self.w[0,:,:], self.w[0,:,:])
        dP_dy = self.central_y(self.P)
        rho_g = self.rho[0,:,:]*(-C.g)
        rho_w = self.rho[0,:,:]*self.w[0,:,:]
        self.drhow_dt = - rho_w*(du_dx_vm+dw_dy_vm) - self.u[0,:,:]*drhow_dx \
                        - self.w[0,:,:]*drhow_dy - dP_dy + rho_g
        self.drhow_dt[:,0] = 0
        self.drhow_dt[:,-1] = 0
        self.w[1,:,1:-1] = (rho_w[:,1:-1] + self.drhow_dt[:,1:-1]*self.dt) / (self.rho[1,:,1:-1])

        # solving the energy equation
        du_dx_e = self.central_x(self.u[0,:,:])
        dw_dy_e = self.central_y(self.w[0,:,:])
        de_dx = self.upwind_x(self.e[0,:,:], self.u[0,:,:])
        de_dy = self.upwind_y(self.e[0,:,:], self.w[0,:,:])
        self.de_dt = (-self.e[0,:,:]-self.P)*(du_dx_e+dw_dy_e) - self.u[0,:,:]*de_dx \
                      - self.w[0,:,:]*de_dy
        self.de_dt[:,0] = 0
        self.de_dt[:,-1] = 0
        self.e[1,:,1:-1] = self.e[0,:,1:-1] + self.de_dt[:,1:-1]*self.dt

        # updating pressure and temperature
        self.T[:,1:-1] = self.e[1,:,1:-1]*(C.gamma-1)*self.mu*C.m_u / (self.rho[1,:,1:-1]*C.k_b)
        self.P[:,1:-1] = (C.gamma-1)*self.e[1,:,1:-1]


        self.timestep()
        self.tot_t += self.dt

        self.rho[0,:,:] = self.rho[1,:,:]
        self.e[0,:,:] = self.e[1,:,:]
        self.u[0,:,:] = self.u[1,:,:]
        self.w[0,:,:] = self.w[1,:,:]
        self.boundary_conditions()

        return self.dt



if __name__ == '__main__':
    print("Perturbation on/off?")
    pert = input()

    vis = FVis.FluidVisualiser()
    solver = convection2D()
    if pert == "on":
        solver.initialise(perturbation=True)
    elif pert == "off":
        solver.initialise()
    else:
        msg = "Did not initialize solver! Need to insert perturbation to on or off!"
        assert False, msg

    print("How long will the simulation last in real time seconds?")
    sim_time = input()

    try:
        sim_time = float(sim_time)
    except:
        msg = "Input must be a real number!"
        assert False, msg

    print("What parameter do you want to plot? Type 'help' for listing of parameters to plot.")
    parameter = input()

    if parameter == "help":
        print("'rho'  -   Mass density")
        print("'drho' -   Mass density contrast")
        print("'u'    -   Horizontal velocity")
        print("'w'    -   Vertical velocity")
        print("'e'    -   Internal energy density")
        print("'de'   -   Internal energy density contrast")
        print("'es'   -   Specific internal energy")
        print("'P'    -   Pressure")
        print("'dP'   -   Pressure contrast")
        print("'T'    -   Temperature")
        print("'dT'   -   Temperature contrast")
        print("'v'    -   Speed")
        print("'ru'   -   Horizontal momentum density")
        print("'rw'   -   Vertical momentum density")
        print("'rv'   -   Momentum density")
        print("'eu'   -   Horizontal energy flux")
        print("'ew'   -   Vertical energy flux")
        print("'ev'   -   Energy flux")
        parameter = input()


    solver.boundary_conditions()
    vis.save_data(sim_time, solver.hydro_solver, rho=solver.rho[0,:,:], e=solver.e[0,:,:], \
                        u=solver.u[0,:,:], w=solver.w[0,:,:], T=solver.T, P=solver.P, sim_fps=5)

    vis.animate_2D(parameter, matrixLike=False, extent=[0,12,0,4], units={"Lx": "Mm", "Lz": "Mm"})
