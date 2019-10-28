import numpy as np
import Interpolation as ip


class Mesh:
    def __init__(self,domain,numOfEls,numOfNodes,numOfQuads):
        self.__xmin = domain[0]
        self.__xmax = domain[1]
        self.__K    = numOfEls
        self.__nN   = numOfNodes
        self.__nQ   = numOfQuads
        self.__ip   = ip.Interpolation(self.__nN, self.__nQ)
        self.__hx   = (self.__xmax-self.__xmin)/self.__K
        
    def Xmin(self):
        return self.__xmin
    def Xmax(self):
        return self.__xmax

    #number of elements    
    def K(self):
        return self.__K
    def Hx(self):
        return self.__hx
    def X(self, nk, inode):
        return self.__hx*nk + self.__ip.Nodes()[inode]*self.__hx
        
        

