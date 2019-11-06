import numpy as np

class Kernel:
    def K(self,x,y,hyper):
        top = 2*(hyper[0]+hyper[1]*x*y)
        bottom = np.sqrt(  (1+2*(hyper[0]+hyper[1]*x*x))*(1+2*(hyper[0]+hyper[1]*y*y))  )
        k = 2/np.pi*np.arcsin( top/bottom )

        return k
