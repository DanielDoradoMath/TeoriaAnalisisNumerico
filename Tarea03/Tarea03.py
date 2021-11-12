import numpy as np
import numpy.linalg as la


def solve_system (xmin, xmax, ymin, ymax, f, g, c, n, m, eps, N_max, omega, imet):
    xk = np.random.random(n*m)
    b = constantVector(xmin, xmax, ymin, ymax, f, g, n, m)
    diag = getDiagonal(xmin, xmax, ymin, ymax, c, n, m)
    k = (ymax-ymin)/(m+1)
    h = (xmax-xmin)/(n+1)
    p = (k/h)**2
    A_dot = lambda vec: A_left(vec, diag, p, n)
    rk = b - A_dot(xk)
    rk_norm = la.norm(rk)


    if imet == 1:
        pass
    elif imet == 2:
        pass
    else:
        num_it = 0
        pk = rk
        while (num_it < N_max) and (rk_norm > eps):
            apk = A_dot(pk)
            rkrk = np.dot(rk,rk)

            alpha = rkrk / np.dot(pk, apk)
            xk += alpha * pk
            rk -= alpha * apk
            beta = np.dot(rk, rk) / rkrk
            pk = rk + beta * pk
            
            num_it += 1
            rk_norm = la.norm(rk)


def getDiagonal(xmin, xmax, ymin, ymax, c, n, m):
    """
    Obtiene los términos de la diagonal de la matriz A.
    """
    x = np.linspace(xmin, xmax, n+2)
    y = np.linspace(ymin, ymax, m+2)
    k = (ymax-ymin)/(m+1)
    h = (xmax-xmin)/(n+1)
    p = (k/h)**2
    alpha = np.zeros(n*m)
    for i in range(1,n+1):
        for j in range(1,m+1):
            alpha[(j-1)*n+i-1] = 2*p + 2 + (k**2)*c(x[i],y[j])
    return alpha


def constantVector(xmin, xmax, ymin, ymax, f, g, n, m):
    """
    Obtiene el vector constante b para el sistema Ax = b.
    """
    x = np.linspace(xmin, xmax, n+2)
    y = np.linspace(ymin, ymax, m+2)
    k = (ymax-ymin)/(m+1)
    h = (xmax-xmin)/(n+1)
    p = (k/h)**2
    b = np.zeros(m*n)
    for i in range(1,n+1):
        for j in range(1,m+1):
            b[(j-1)*n+i-1] = k**2*f(x[i],y[j])
            if i == 1:
                b[(j-1)*n+i-1] += p*g(x[i-1], y[j])
            if j == 1:
                b[(j-1)*n+i-1] += g(x[i], y[j-1])
            if i == n:
                b[(j-1)*n+i-1] += p*g(x[i+1], y[j])
            if j == m:
                b[(j-1)*n+i-1] += g(x[i], y[j+1])
    return b

def A_left(x, diag, p, n):
    s = x.size
    prod = np.zeros(s)
    for i, comp in enumerate(x):
        r = diag[i]*comp
        if i>1:
            r += -p*x[i-1]
        if i<s-1:
            r += -p*x[i+1]
        if i+n<s:
            r -= x[i+n]
        if i-n>0:
            r -= x[i-x]
        prod[i] = r
    return prod
