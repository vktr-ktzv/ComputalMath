import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
import math


def Matrix(N, a): # A
    MatrixArray = np.zeros(shape = (N, N))
    for i in range(N):
        MatrixArray[i][i] = 2
        if (i != (N - 1)):
            MatrixArray[i][i + 1] = -1 - a
        if (i != 0):
            MatrixArray[i][i - 1] = -1 + a
    return MatrixArray

def RightVector(N, a): # f
    f = np.zeros(shape = (N, 1))
    f[0][0] = 1 - a
    f[N - 1][0] = 1 + a
    return f

def Decompos(N, A): # A = L + D + U
    L = np.zeros(shape = (N, N))
    D = np.zeros(shape = (N, N))
    U = np.zeros(shape = (N, N))
    for i in range(N):
        for j in range(i):
            L[i][j] = A[i][j]
    for i in range(N):
        for j in range(i):
            U[j][i] = A[j][i]
    for i in range(N):
        D[i][i] = A[i][i]
    return (L, D, U)

def NewX(LDinv, f, U, Xk): # новый итерационный элемент 
    Xk1 = np.matrix(LDinv) * np.matrix(f - np.dot(U, Xk))
    return Xk1

def solveZeid(N, Mtx, f, Eps): # решение 
    L, D, U = Decompos(N, Mtx)
    Cnt = 0
    LDinv = np.linalg.inv(L + D)
    Xk = np.zeros(shape = (N, 1))
    Xk1 = NewX(LDinv, f, U, Xk)
    while (norm((Xk1 - Xk), ord=2) > Eps):
        Cnt += 1
        Xk = Xk1
        Xk1 = NewX(LDinv, f, U, Xk)
    return (np.reshape(Xk1, (1, N)), Cnt)


N = 15
a = 0.1
Eps = 1e-6
A = Matrix(N, a)
f = RightVector(N, a)
X, Cnt = solveZeid(N, A, f, Eps)


print("Решение при n =", N, ", a = ", a, ":\n X = \n", X, "\n")

plt.figure(100)
plt.title("зависимость итераций от a для n = 15")

plt.xlabel('a')
plt.ylabel('N iter')

Alphas = np.arange(-1.4, 1.5, 0.1)
Cnt = len(Alphas)
Y = np.zeros(Cnt)
for i in range(Cnt):
    
    A = Matrix(N, Alphas[i])
    f = RightVector(N, Alphas[i])
    X, Cnt = solveZeid(N, A, f, Eps)
    Y[i] = Cnt

plt.scatter(Alphas, Y)
plt.plot(Alphas, Y, 'r')

plt.grid()
plt.show()