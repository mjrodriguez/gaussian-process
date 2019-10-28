import numpy as np
import Interpolation as ip
import Mesh as mesh
import SimulationParameters as sp

if __name__ == "__main__":

    params = sp.Parameters()
    interp = ip.Interpolation(params.Nnodes(), params.Nquads())
    m      = mesh.Mesh(params.Domain(), params.Nels(), params.Nnodes(), params.Nquads())

    