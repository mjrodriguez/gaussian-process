import numpy as np
import linedg as linedg

class time_integration:
    def __init__(self, params, linedg):
        self.__linedg = linedg
        self.__CourantNumber = params.CourantNumber()
        self.__integrator = params.TimeIntegration()

    
    def Evolve(self, dt, t, uold):
        if (self.__integrator == "rk4"):
            unew = self.RK4(dt,t,uold)
        return unew
    def RK4(self,dt,t,uold):
        k1 = dt*self.__linedg.AssembleElement(uold)
        k2 = dt*self.__linedg.AssembleElement(uold + 0.5*k1)
        k3 = dt*self.__linedg.AssembleElement(uold + 0.5*k2)
        k4 = dt*self.__linedg.AssembleElement(uold + k3)

        unew = uold + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0
        return unew
 
    def ComputeDt(self,u):
        umax = np.amax(u)
        return self.__CourantNumber*self.__linedg.m.dx()/umax/(self.__linedg.order()**2)
