import numpy as np
import math
import matplotlib.pyplot as plt
import sympy as sp
from scipy. integrate import odeint

def Jac(f, x, dx=1e-6):
    n = len(x)
    func = f(x)
    jac = np.zeros((n, n))
    for j in range(n):
        Dxj = (abs(x[j]) * dx if x[j] != 0 else dx)
        x_plus = [(xi if k != j else xi + Dxj) for k, xi in enumerate(x)]
        jac[:, j] = (f(x_plus) - func)/Dxj
    return jac, func

def Newton(f, x, eps = 1.0e-2):
    max_iter = 1000
    a = np.array([[1, 2, 8], [3, 5, 0], [3, 2, 1]])
    b = np.array([1, 2, 6])
    j = np.linalg.solve(a, b)
    
    for i in np.arange( 0, max_iter, 1):
        J, fn = Jac( f, x )
        if np.sqrt( np.dot( fn, fn ) / x.size ) < eps:
            return x
        dx = np.linalg.solve( J, np.array(fn.transpose()) )
        x = x - dx
    
def Iteration( f, y0, tSTART, tEND, tau, alpha ):
    def F(y_next):
        return y_next - tau * alpha * f(t[i], y_next) - y[i] - tau * ( 1. - alpha) * f( t[i], y[i] )
    t = np.arange( tSTART, tEND, tau )
    y = np.ones((t.size, 3))
    y[0] = y0
    for i in np.arange(0, t.size-1, 1):
        y_next = np.array([0,0,0])
        y_next = y[ i ] + tau * f( t[ i ], y[ i ] )
        y[ i + 1 ] = Newton( F, y_next )
    return t, y

def f(t, u):
    x = u[ 0 ]
    y = u[ 1 ]
    a = u[ 2 ]
    DerX = x * (1 - 0.5 * x - 2. * y / (7. * a * a)) 
    DerY = x * (2. * a - 3.5 * a * a * x - 0.5 * y)
    DerA = (2 - 7. * a * x) / 100
    return np.array([DerX, DerY, DerA])

    
def main():
    plt.figure(figsize = (15, 8))
    plt.title("Решение системы неявным методом Адамса 2-го порядка:")
    tSTART= 0.
    tEND = 0.1
    tau = 0.001
    y0 = np.array([ 1.5, 10., 0.1 ])
    alpha = 0.5
    t, y = Iteration( f, y0, tSTART, tEND, tau, alpha )

    for n in np.arange(0,3,1):
      u = []
      for i in range(int(y.size / 3)):
        u.append(y[i][n])
      plt.plot( t, u, label = str(n))

    plt.legend()
    plt.grid()
    plt.show()

if __name__ == '__main__':
    main()