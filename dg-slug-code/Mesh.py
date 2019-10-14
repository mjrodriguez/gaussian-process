import numpy as np
import Interpolation as ip


class Mesh:
    def __init__(self,domain,numOfEls,numOfNodes,numOfQuads):
        self.__xmin = domain[0]
        self.__xmax = domain[1]
        self.__K    = numOfEls
        self.__nN     = numOfNodes
        self.__nQ     = numOfQuads

    
    def __ConstructGlobalCoordinates(self):
        self.__hx = (self.__xmax-self.__xmin)/self.__K
        

