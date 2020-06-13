import numpy as np
import Interpolation as interpolation

class linedg:
    def __init__(self, mesh, P, u):
        self.P = P
        self.__dim   = 1
        self.__neqs  = P.neqs()
        self.__order = P.order()
        self.__nnodes = self.__order+1
        self.__nels  = P.nels()
        self.__nquads = P.nquads()
        self.ip = interpolation.Interpolation(self.__nnodes, self.__nquads)
        self.m  = mesh
        self.mass = np.matmul(np.transpose(self.ip.B()), np.matmul(self.ip.W(), self.ip.B()) )
        self.invMass = np.linalg.inv(self.mass)
        self.u = u

        self.__q = np.zeros([self.__nels,self.__nnodes])


    # def InitialCondition(self):
    #     u = np.ones([self.__nels, self.__nnodes])
    #     return u
    def q(self):
        return self.__q

    def Flux(self,u):
        # Compute flux
        F = 0.5*np.multiply(u,u)
        return F*self.m.J()*self.m.invJ()

    def SetBC(self):
        if (self.P.BoundaryConditions() == "periodic"):
            self.__leftBC = self.u[self.__nels-1, self.__nnodes-1]
            self.__rightBC = self.u[0,0]
        elif (self.P.BoundaryConditions() == "nonperiodic"):
            self.__leftBC = P.LeftBC()
            self.__rightBC = P.RightBC()

    def AssembleElement(self):
        self.SetBC()
        q = np.zeros([self.__nnodes,1])
        fstar = np.zeros([2,1])
        for iel in range(self.__nels):
            # u at quadrature points
            U = np.matmul( self.ip.B(), self.u[iel,:])
            # Compute reference flux
            F = self.Flux(U)
            #Compute volume term
            q = -np.matmul( np.transpose(self.ip.D()), np.matmul(self.ip.W(), F) )
            if (iel == 0):
                fstar[0] = self.RiemannSolver(self.__leftBC,self.u[iel,0])
                fstar[-1] = self.RiemannSolver(self.u[iel,-1], self.u[iel+1,0])
            elif (iel == self.__nels-1):
                fstar[0] = self.RiemannSolver(self.u[iel-1,-1], self.u[iel,0])
                fstar[-1] =self.RiemannSolver(self.u[iel,-1], self.__rightBC)
            else:
                fstar[0] = self.RiemannSolver(self.u[iel-1,-1], self.u[iel,0])
                fstar[-1] = self.RiemannSolver(self.u[iel,-1], self.u[iel+1,0])
            
            q[0] -= fstar[0]
            q[-1] += fstar[1]
            self.__q[iel,:] = np.matmul(self.invMass, q)
    
    ######################################################
    #
    # Riemann Solver Functions
    #
    #######################################################

    def RiemannSolver(self, uleft, uright):
        rs = self.P.RiemannSolver()
        if (rs == "upwind"):
            c = np.max(self.u)
            fhat = self.Upwind(c,uleft,uright)
        return fhat

    def Upwind(self, c, uleft, uright):
        if (c >= 0 ):
            ustar = uleft
        elif (c < 0):
            ustar = uright
        fstar = self.Flux(ustar)
        return fstar
        


        

    

    



    