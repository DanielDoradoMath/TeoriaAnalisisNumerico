"""
Tarea 1 - Teoría del Análisis Numérico

Realizado por: José Darío Flórez & Daniel Dorado Toro

Solución al problema 2.

"""

import numpy as np
import sys


# Parámetros a ingresar por consola

f = eval(sys.argv[1])   # Parámetro 1: función a interpolar
a = float(sys.argv[2])   # Parámetro 2: inicio del intervalo
b = float(sys.argv[3])   # Parámetro 3: final del intervalo
n = int(sys.argv[4])   # Parámetro 4: número de puntos de interpolación
x = float(sys.argv[5])   # Parámetro 5: punto a evaluar

"""
Para ingresar la función se debe usar la sintaxis lambda de Python.
Ejemplo:
Para interpolar coseno: python Tarea01.py "lambda x: np.cos(x)" <a> <b> <n> <x>

Las comillas son necesarias para que la función sea tomada como un solo
parámetro y no varios.
"""


pts_chev = []    # Define una lista para los puntos de Chebyshev
pts_equi = []    # Define una lista para los puntos equidistantes


# Cálculo de los puntos equidistantes

punto = a
for i in range(n):
    pts_equi.append(punto)
    punto += ((b - a) / n)

# Cálculo de los puntos de Chebyshev

punto = 0
for i in range(n):
    punto = np.cos(((2*i+1)*np.pi)/(2*n))
    punto = ((b-a)*punto+a+b)/2
    pts_chev.append(punto)

# Interpolación
res_e = 0    # Resultado de interpolación con puntos equidistantes
res_c = 0    # Resultado de interpolación con puntos de Chebyshev

for i in range(len(pts_equi)):
    prod_e = 1    # Producto para los puntos equidistantes
    prod_c = 1    # Producto para los puntos de Chebyshev
    for j in range(len(pts_equi)):
        if j != i:
            prod_e *= (x-pts_equi[j])/(pts_equi[i]-pts_equi[j])
            prod_c *= (x-pts_chev[j])/(pts_chev[i]-pts_chev[j])
    res_e += f(pts_equi[i])*prod_e
    res_c += f(pts_chev[i])*prod_c

# Cálculo del error de interpolación
err_e = np.abs(res_e-f(x))
err_c = np.abs(res_c-f(x))

# Respuesta en consola:
print('Punto evaluado (x): {}'.format(x))
print('Valor de la función: {:.6f}'.format(f(x)))
print('Resultado de interpolación (equidistantes): {}'.format(res_e))
print('Error de interpolación (equidistantes): {}'.format(err_e))
print('Resultado de interpolación (Chebyshev): {}'.format(res_c))
print('Error de interpolación (Chebyshev): {}'.format(err_c))
