import numpy as np
import numpy.linalg as la


def solvePDE (xmin, xmax, ymin, ymax, f, g, c, n, m, eps, N_max, omega, imet):
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
        num_it = 0
        while (num_it < N_max) and (rk_norm > eps):
            rk_norm = 0
            xk = -A_left(xk, np.zeros(xk.size), p, n) + b
            for i in range(xk.size):
                xk[i] *= 1/(diag[i])
            
            rk = b - A_dot(xk)
            rk_norm = la.norm(rk)
            num_it += 1
        printResult(num_it, N_max, xk, rk_norm, n, m)
    
    elif imet == 2:
        s = xk.size
        num_it=0
        while (num_it < N_max) and (rk_norm > eps):
            for i in range(s):
                r = 0
                if i>=1:
                    r -= p*xk[i-1]
                if i<s-1:
                    r -= p*xk[i+1]
                if i+n<s:
                    r -= xk[i+n]
                if i-n>=0:
                    r -= xk[i-n]
                q = b[i] - r
                q *= omega/diag[i]
                xk[i] = (1-omega)*xk[i] + q
            
            rk = b - A_dot(xk)
            rk_norm = la.norm(rk)
            num_it+=1
        printResult(num_it, N_max, xk, rk_norm, n, m)

    elif imet == 3:
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
        printResult(num_it, N_max, xk, rk_norm, n, m)
    else:
        print('Método no implementado.')


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
        if i>0:
            r -= p*x[i-1]
        if i<s-1:
            r -= p*x[i+1]
        if i+n<s:
            r -= x[i+n]
        if i-n>=0:
            r -= x[i-n]
        prod[i] = r
    return prod


def printResult(num_it, N_max, vec, res_norm, n, m):
    if num_it > N_max:
        print('No hubo convergencia en el número de iteraciones máximo.')
        print('Norma del residuo final: {:.10f}.'.format(res_norm))
        print('Tabla de resultados:')
        print('|   i   |   j   |  k(i,j)  |  mu(i,j)  |')
        print('|-------|-------|----------|-----------|')
        for j in range(1, m+1):
            for i in range(1, n+1):
                k = (j-1)*n+i
                print('|{:^7d}|{:^7d}|{:^10d}|{:^11.3f}|'.format(i, j, k, vec[k-1]))
    else:
        print('Número de iteraciones realizadas: {}'.format(num_it))
        print('Norma del residuo final: {:.10f}.'.format(res_norm))
        print('Tabla de resultados:')
        print('|   i   |   j   |  k(i,j)  |  mu(i,j)  |')
        print('|-------|-------|----------|-----------|')
        for j in range(1, m+1):
            for i in range(1, n+1):
                k = (j-1)*n+i
                print('|{:^7d}|{:^7d}|{:^10d}|{:^11.3f}|'.format(i, j, k, vec[k-1]))
