import numpy as np
import Interpolation as Interpolation

class Limit:
    def __init__(self, params, linedg):
        self.__p = params
        self.__lim = params.LimitSolution()
        self.__nnodes = params.nnodes()
        self.__nquads = params.nquads()
        self.__nels   = params.nels()
        self.__dx     = linedg.m.dx()
        self.ip = Interpolation.Interpolation(self.__nnodes, self.__nquads)

    def minmod(self,ubar):
        s = np.sum(np.sign(ubar))/ubar.size
        if (np.abs(s) == 1):
            m = s*np.amin(np.abs(ubar))
        else:
            m = 0
        return m

    def FindElementsToLimit(self):
        # Compute Cell averages
        self.ubar = np.zeros([self.__nels+2,1])
        for i in range(self.__nels+2):
            self.ubar[i] = np.sum( np.matmul( self.ip.W(),np.matmul(self.ip.B(),self.u[i,:]) ) )
        
        # compute interface fluxes
        ids = np.zeros([self.__nels,1])
        for i in range(1,self.__nels+1):
            vl = self.ubar[i] - self.minmod( np.array([ self.ubar[i]-self.u[i-1,0], self.ubar[i] - self.ubar[i-1], self.ubar[i+1] - self.ubar[i]]) )
            vr = self.ubar[i] + self.minmod( np.array([ self.ubar[i]-self.u[i-1,-1], self.ubar[i] - self.ubar[i-1], self.ubar[i+1] - self.ubar[i] ]) )

            
            if ( (np.abs(vl-self.u[i,0]) < 1e-8) and (np.abs(vr-self.u[i,-1]) < 1e-8) ):
                j = i-1
                ids[j] = 1
        return ids

    def LimitSolution(self,u):
        self.u_og = u
        self.u = np.zeros([self.__nels+2, self.__nnodes])
        self.u[1:self.__nels+1,:] = u

        if (self.__p.BoundaryConditions() == "periodic"):
            self.u[0,:] = u[-1,:]
            self.u[-1,:] = u[0,:]
        else:
            self.u[0,:] = np.ones([self.__nnodes])*self.__p.LeftBC()
            self.u[-1,:] = np.ones([self.__nnodes])*self.__p.RightBC()
        
        ids = self.FindElementsToLimit()
        ids[0] = 1
        self.__el_lim_ids = np.nonzero(ids > 0)
        if (self.__el_lim_ids[0].size != 0):
            self.PI1()
        # print(self.__el_lim_ids[0])
    
    def PI1(self):
        self.dubar = np.zeros([self.__nels+2,1])
        for i in range(self.__nels+2):
            self.dubar[i] = np.sum( np.matmul( self.ip.W(), np.matmul(self.ip.D(), self.u[i,:])))/self.__dx
        
