import numpy as np
import numpy.linalg as la


def LinearCG(A, b, x0, eps=1e-5):
    xk = x0
    rk = b - np.dot(A, xk)
    pk = rk
    rk_norm = la.norm(rk)
    
    num_iter = 0
    while rk_norm > eps:
        apk = np.dot(A, pk)
        rkrk = np.dot(rk, rk)
        
        alpha = rkrk / np.dot(pk, apk)
        xk += alpha * pk
        rk -= alpha * apk
        beta = np.dot(rk, rk) / rkrk
        pk = rk + beta * pk
        
        num_iter += 1
        rk_norm = la.norm(rk)
        print('Iteration: {} \t x = {} \t residual = {:.4f}'.
              format(num_iter, xk, rk_norm))
    
    print('\nSolution: \t x = {}'.format(xk))
