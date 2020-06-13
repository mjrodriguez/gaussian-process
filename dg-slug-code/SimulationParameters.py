import numpy as np

class Parameters:
    def __init__(self, neqs=1, order=2, quads=3*2+1, numOfEls=30, domain=np.array([0,1]), courantNumber=0.8, maxTime=0.25, timeIntegration="rk4" ,boundaryCondition="periodic", riemannSolver="godunov", limitBool="False"):
        self.__order = order
        self.__neqs  = neqs
        self.__p     = order+1
        self.__quads = quads
        self.__nK    = numOfEls
        self.__domain = domain
        self.__courantNumber     = courantNumber
        self.__maxTime = maxTime
        self.__timeIntegration = timeIntegration
        self.__boundaryCondition = boundaryCondition
        self.__limit  = limitBool
        self.__riemannSolver = riemannSolver
        if (boundaryCondition == "nonperiodic"):
            self.__leftBC = 0
            self.__rightBC = 1
        
    def neqs(self):
        return self.__neqs
    def order(self):
        return self.__order
    def nnodes(self):
        return self.__p
    def nquads(self):
        return self.__quads
    def nels(self):
        return self.__nK
    def domain(self):
        return self.__domain
    def LeftBC(self):
        return self.__leftBC
    def RightBC(self):
        return self.__rightBC
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