import numpy as np
import math
import matplotlib.pyplot as plt
from numpy.linalg import norm

def Matrix(): #A
    MatrixArray = np.zeros(shape = (2, 2))
    MatrixArray[0][0] = 0.78
    MatrixArray[0][1] = 0.563
    MatrixArray[1][0] = 0.457
    MatrixArray[1][1] = 0.33
    return MatrixArray

def Norm3(A, Atr):
    Mult = np.dot(Atr, A)
    MaxL = 0
    lambdas = np.linalg.eigvals(Mult)
    for l in lambdas:
        if (l > MaxL):
            MaxL = l
    return math.sqrt(MaxL)

def Mu(): #число обусловленности 
    Mtx = Matrix()
    RMtx = np.linalg.inv(Mtx)
    return Norm3(Mtx, np.transpose(Mtx)) * Norm3(RMtx, np.transpose(RMtx))

def NewX(A, Xk, f): #новый элемент для следующей итерации в мпи
    Xk1 = (np.eye(2) - np.matrix(A)) * Xk + f
    return Xk1

def Solution(f, Epsilon): # Ax = f = f0 + dfN = RightVector(N)
    A = Matrix()
    Xk = np.matrix('1; 1')
    Xk1 = NewX(A, Xk, f)
    while (norm((Xk1 - Xk), ord=2) > Epsilon):
        Xk = Xk1
        Xk1 = NewX(A, Xk, f)
    return Xk1

Eps = 1e-6
X0 = np.matrix('1; -1')
NormX0 = norm(X0, ord=2)
Mu_ = Mu()

f0 = np.matrix('0.217; 0.127')
NormF0 = norm(f0, ord=2)
df1 = np.matrix('0; 0.0005')
df2 = np.matrix('0.0001; 0')
df3 = np.matrix('0.001; 0.0006')
F = [f0, f0 + df1, f0 + df2, f0 + df3]

for i in range(1, 4):
    X = Solution(F[i], Eps)
    dX = norm(X - X0, ord=2)
    print("f",i, ":\nОтвет X = (", format(X[0, 0], '.3f'), ",", format(X[1, 0], '.3f'), ")", "\n||X - X0|| = ", format(dX,  '.3f'))
    print("Проверка неравенства:", format(dX / NormX0, '.2f'), " <= ", format(Mu_ * norm(F[i] - f0, ord=2) / NormF0, '.2f'), "\n")
