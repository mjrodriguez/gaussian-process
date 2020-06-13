import numpy as np
import matplotlib.pyplot as plt
import SimulationParameters as sp
import Mesh as m
import linedg
from scipy import sparse
from mpl_toolkits.mplot3d import Axes3D
import sys
np.set_printoptions(linewidth=250)

if __name__ == "__main__":
    neqs = 1
    order = 2
    nquads = 2*order
    nels = 10
    domain = np.array([-1,1])
    maxtime = 0.25
    cfl = 0.8

    params = sp.Parameters(neqs, order, nquads, nels, domain, cfl, maxtime, "rk4", "periodic", "upwind", False)
    mesh = m.Mesh(params.domain(), params.nels(), params.nnodes(), params.nquads())
    u = 0.5 + np.sin(np.pi*mesh.X())
    print(mesh.X())
    print(u)
    ldg = linedg.linedg(mesh, params, u)
    ldg.AssembleElement()
    # print("q = ", ldg.q())

    


    

