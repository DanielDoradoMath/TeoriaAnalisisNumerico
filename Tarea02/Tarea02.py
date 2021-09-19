"""
Tarea 1 - Teoría del Análisis Numérico

Realizado por: José Darío Flórez & Daniel Dorado Toro

Solución al problema 2.

"""
import numpy as np
import sympy as sp
import sys


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
        pass
    elif method == 'symb':
        pass
    else:
        print('Método no implementado')


def redef(func, m1, m2, m3):
    """
    Redefine la función en términos de u y v.
    """
    pass


def integrate(func, m1, m2, m3, method='num'):
    """
    Integra la función func sobre la región triangular con vértices en m1, m2
    y m3
    """
    pass
