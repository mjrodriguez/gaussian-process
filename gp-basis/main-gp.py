import numpy as np
import matplotlib.pyplot as plt
import gl

# Define the kernel function
def kernel(a, b, param):
    sqdist = np.sum(a**2,1).reshape(-1,1) + np.sum(b**2,1) - 2*np.dot(a, b.T)
    return np.exp(-.5 * (1/param) * sqdist)

if __name__ == "__main__":

    # Creating Data Set
    order = 4
    n = order+1
    xnodes = gl.lglnodes(order)
    xquad  = gl.lglnodes(2*order)
    B,D = gl.lagint(xnodes[0],xquad[0])

    Xtrain = xnodes[0].reshape(-1,1)
    # Xtest = np.linspace(0, 1, order+1).reshape(-1,1)
    Xtest = Xtrain
    
    print("---- xtrain ----")
    print(Xtrain)
    print("---- xtest ----")
    print(Xtest)
    
    # Lagrange polynomial
    ytrain = np.sum(B,axis=0)
    print("---- ytrain ----")
    print(ytrain)

    # Gaussian Process
    param = 0.1
    K_ss = kernel(Xtest, Xtest, param) 

    # Get cholesky decomposition (square root) of the
    # covariance matrix
    # Apply the kernel function to our training points
    K = kernel(Xtrain, Xtrain, param)
    L = np.linalg.cholesky(K + 0.00005*np.eye(len(Xtrain)))

    # Compute the mean at our test points.
    K_s = kernel(Xtrain, Xtest, param)
    Lk = np.linalg.solve(L, K_s)
    mu = np.dot(Lk.T, np.linalg.solve(L, ytrain)).reshape((n,))

    # Compute the standard deviation so we can plot it
    s2 = np.diag(K_ss) - np.sum(Lk**2, axis=0)
    stdv = np.sqrt(s2)
    # Draw samples from the posterior at our test points.
    L = np.linalg.cholesky(K_ss + 1e-6*np.eye(n) - np.dot(Lk.T, Lk))
    f_post = mu.reshape(-1,1) + np.dot(L, np.random.normal(size=(n,3)))

    plt.plot(Xtrain, ytrain, 'bs', ms=8)
    plt.plot(Xtest, f_post)
    plt.gca().fill_between(Xtest.flat, mu-2*stdv, mu+2*stdv, color="#dddddd")
    plt.plot(Xtest, mu, 'r--', lw=2)
    plt.axis([0,1, 0, 3])
    plt.title('Three samples from the GP posterior')
    plt.show()



    
