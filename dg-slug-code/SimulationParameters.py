import numpy as np

class Parameters:
    def __init__(self, order=2, quads=3*2+1, numOfEls=30, domain=np.array([0,1]), courantNumber=0.8, maxTime=0.25, timeIntegration="rk4" ,boundaryCondition="periodic", riemannSolver="godunov", limitBool="False"):
        self.__order = order
        self.__p     = order+1
        self.__quads = 3*self.__order+1
        self.__nK    = numOfEls
        self.__domain = domain
        self.__courantNumber     = courantNumber
        self.__maxTime = maxTime
        self.__timeIntegration = timeIntegration
        self.__boundaryCondition = boundaryCondition
        self.__limit  = limitBool
        self.__riemannSolver = riemannSolver
        

    def Order(self):
        return self.__order
    def Nnodes(self):
        return self.__p
    def Nquads(self):
        return self.__quads
    def Nels(self):
        return self.__nK
    def Domain(self):
        return self.__domain
    def CourantNumber(self):
        return self.__courantNumber
    def MaxTime(self):
        return self.__maxTime
    def BoundaryConditions(self):
        return self.__boundaryCondition
    def LimitSolution(self):
        return self.__limit
    def RiemannSolver(self):
        return self.__riemannSolver