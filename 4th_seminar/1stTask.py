import numpy as np
import math
import matplotlib.pyplot as plt

def MatrixA(N): #A
    MatrixArray = np.zeros(shape = (N, N))
    for i in range(N - 1):
        MatrixArray[i][i] = 1
        for j in range(i + 1, N):
            MatrixArray[i][j] = -1
    for j in range(N):
        MatrixArray[N - 1][j] = 1
    return MatrixArray

def RightVector(N): #f
    return np.ones(shape = (N, 1))

def LU(N, A): #разложение LU
    L = np.zeros(shape = (N, N))
    U = A
    for i in range(N):
        for j in range(i, N):
            L[j][i] = U[j][i] * 1.0 / U[i][i]

    for k in range(1, N):
        for i in range(k-1, N):
            for j in range(i, N):
                L[j][i] = U[j][i] * 1.0 / U[i][i]
        for i in range(k, N):
            for j in range(k-1, N):
                U[i][j] = U[i][j] - L[i][k-1] * U[k-1][j]
    return (L, U)

def Norm3(A, Atr): #Third norm
    Mult = np.dot(A, Atr)
    MaxL = 0
    eigvals = np.linalg.eigvals(Mult)
    for l in eigvals:
        if (l > MaxL):
            MaxL = l
    return math.sqrt(MaxL)

def Mu(N): #число обусловленности 
    A = MatrixA(N)
    Atr= np.linalg.inv(A)
    return Norm3(A, np.transpose(A)) * Norm3(Atr, np.transpose(Atr))


def solveU(N, U, f): # решение UX = f
    Sol= np.zeros(N)
    Sol[N - 1] = f[N - 1, 0] * 1.0 / U[N - 1][N - 1]
    for n in range(N - 2, -1, -1):
        Sum = 0
        for i in range(n + 1, N):
            Sum += U[n][i] * Sol[i]
        Sol[n] = (f[n, 0] - Sum) / U[n][n]
    return Sol

def solve(N): # решение матрицы N на N
    f = RightVector(N)
    A = MatrixA(N)
    L, U = LU(N, A)
    f = np.dot(np.linalg.inv(L), f)
    return solveU(N, U, f)

plt.figure(100)
plt.title("Число обусловленности")

plt.xlabel('n')
plt.ylabel('Mu')

N = 10
X = np.arange(N)
Y = [Mu(i) for i in X]

A, B = np.polyfit(X, Y, 1)
YL = [A * x + B for x in X]

for i in range(1, N + 1):
    print("n = ", i, " A : \n", MatrixA(i), "\n Решение X =", solve(i), "\n")

plt.scatter(X, Y)
plt.plot(X, YL, 'y')

plt.grid()
plt.show()