import numpy as np


def solve_system (xmin, xmax, ymin, ymax, f, g, c, n, m, eps, N_max, omega, imet):
    xk = np.random.random(n*m)
    b = constantVector(xmin, xmax, ymin, ymax, f, g, n, m)
    alpha = getConstants(xmin, xmax, ymin, ymax, c, n, m)
    k = (ymax-ymin)/(m+1)
    h = (xmax-xmin)/(n+1)
    p = (k/h)**2
    rk = 2
    rk_norm = np.linalg.norm()


def getConstants(xmin, xmax, ymin, ymax, c, n, m):
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


# def getInitValueMatrix(xmin, xmax, ymin, ymax, g, n, m):
#     """
#     Dado un rectángulo y una función definida sobre su frontera, esta función
#     retorna una matriz de Numpy de tamaño n x m con los valores de g en los
#     bordes y valores aleatorios al interior.
#     Esta matriz da los valores de x0 para los métodos iterativos.
#     """
#     x = np.linspace(xmin, xmax, n+2)
#     y = np.linspace(ymin, ymax, m+2)
#     u = np.ones((m+2, n+2))
#     for i in range(n+2):
#         for j in range(m+2):
#             if i*j == 0 or i==n+1 or j==m+1:
#                 u[j,i] = g(x[i], y[j])
#             else:
#                 u[j,i] = np.random.random()
#     return u

# def getVector(u):
#     """
#     Dada una matriz de valores, obtiene el vector x que aproxima la solución 
#     al sistema Ax=b.
#     """
#     n = u[0].size - 2
#     m = u[:,0].size - 2
#     x = np.zeros(n*m)
#     for j in range(1,m+1):
#         for i in range(1,n+1):
#             x[(j-1)*n+i-1] = u[j,i]
#     return x


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
