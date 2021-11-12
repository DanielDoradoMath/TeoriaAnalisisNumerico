import numpy as np

# N_max = 70
# b = np.array([1,1,1,1])

# r = b - 

# i = 0

# while i < N_max and :
#     i += 1

    


def A_prod (x,y, A):
    prod = 0
    for i in range(len(x)):
        for j in range(len(y)):
            prod += y[j] * A[j][i] * x[i]
    return prod


def LinearCG(A, b, x0, eps=1e-5):
    xk = x0
    rk = np.dot(A, xk) - b
    pk = -rk
    rk_norm = np.linalg.norm(rk)
    
    num_iter = 0
    while rk_norm > eps:
        apk = np.dot(A, pk)
        rkrk = np.dot(rk, rk)
        
        alpha = rkrk / np.dot(pk, apk)
        xk = xk + alpha * pk
        rk = rk + alpha * apk
        beta = np.dot(rk, rk) / rkrk
        pk = -rk + beta * pk
        
        num_iter += 1
        rk_norm = np.linalg.norm(rk)
        print('Iteration: {} \t x = {} \t residual = {:.4f}'.
              format(num_iter, xk, rk_norm))
    
    print('\nSolution: \t x = {}'.format(xk))
