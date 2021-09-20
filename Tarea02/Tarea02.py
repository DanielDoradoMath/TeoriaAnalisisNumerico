"""
Tarea 2 - Teoría del Análisis Numérico

Realizado por: José Darío Flórez & Daniel Dorado Toro

"""
import numpy as np
import sympy as sy
# import sys


def jacobian(m1, m2, m3):
    """
    Calcula el jacobiano de la transformación que envía la región triangular
    cuyos vertices son (0,0), (1,0) y (0,1) a la región triangular con
    vértices en m1, m2 y m3.
    m1, m2, m3 = 2-tuplas con los puntos del triángulo
    """
    x1 = m1[0]
    x2 = m2[0]
    x3 = m3[0]

    y1 = m1[1]
    y2 = m2[1]
    y3 = m3[1]

    J = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)
    J = np.abs(J)
    return J


def integrateOverT(func, method='num'):
    """
    Integra la función func sobre la región triangular con vértices en (0,0),
    (1,0) y (0,1).
    """
    if method == 'num':
        V = 1/2
        t = 1/3
        A = (9*V)/40
        r = (6-np.sqrt(15))/21
        s = (9+2*np.sqrt(15))/21
        u = (6+np.sqrt(15))/21
        v = (9-2*np.sqrt(15))/21
        B = V*((155-np.sqrt(15))/1200)
        C = V*((155+np.sqrt(15))/1200)

        r = A*func(t, t) + B*func(r, r) + B*func(s, r) + B*func(r, s)
        r += C*func(u, u)+C*func(u, v) + C*func(v, u)

        return r
    elif method == 'symb':
        V = 1/2
        t = 1/3
        A = (9*V)/40
        r = (6-sy.sqrt(15))/21
        s = (9+2*sy.sqrt(15))/21
        u = (6+sy.sqrt(15))/21
        v = (9-2*sy.sqrt(15))/21
        B = V*((155-sy.sqrt(15))/1200)
        C = V*((155+sy.sqrt(15))/1200)

        r = A*func(t, t) + B*func(r, r) + B*func(s, r) + B*func(r, s)
        r += C*func(u, u)+C*func(u, v) + C*func(v, u)

        return r
    else:
        print('Método no implementado')


def redef(func, m1, m2, m3):
    """
    Redefine la función en términos de u y v.
    """
    x1 = m1[0]
    x2 = m2[0]
    x3 = m3[0]

    y1 = m1[1]
    y2 = m2[1]
    y3 = m3[1]

    a = x2 - x1
    b = x3 - x1

    c = y2-y1
    d = y3-y1

    return lambda u, v: func(a*u + b*v + x1, c*u + d*v + y1)


def integrate(func, m1, m2, m3, method='num'):
    """
    Integra la función func sobre la región triangular con vértices en m1, m2
    y m3
    """
    fT = redef(func, m1, m2, m3)
    return integrateOverT(fT, method)


def f(x, i, j):
    return (x ** i * (1-x)**(j+1))/(j+1)


def twoVarPoly(x, y, i, j):
    return (x**i) * (y**j)


x = sy.Symbol("x")

for i in range(6):
    j = 0
    while j <= 5 and i+j < 6:
        iint = sy.integrate(f(x, i, j), (x, 0, 1))
        poly = lambda x, y: twoVarPoly(x, y, i, j)
        iint_aprox = integrateOverT(poly)
        print('Función a integrar: (x^{})*(y^{})'.format(str(i), str(j)))
        print('Integral exacta: ' + str(iint))
        print('Integral por cuadratura: ' + str(iint_aprox))
        print('---')
        j += 1

iint = sy.integrate(f(x, 6, 0), (x, 0, 1))
poly = lambda x, y: twoVarPoly(x, y, 6, 0)
iint_approx = integrateOverT(poly)
print('Función a integrar: x^6')
print('Integral exacta: ' + str(iint))
print('Integral por cuadratura: ' + str(iint_aprox))


def integrand(x, y):
    return (np.sin(3*x+y**2))**2


iint = integrate(integrand, (-1, 2), (2, -1), (3, 3))
print('Función a integrar: sin^2(3*x + y^2)')
print('Dominio de integración: triángulo con vértices (-1,2), (2,-1) y (3,3)')
print('Integral por cuadratura: ' + str(iint))
