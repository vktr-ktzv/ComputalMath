import numpy as np
import math

from numpy.linalg import norm

def GetMatrix(N): #A
    MatrixArray = np.zeros(shape = (N, N))
    
    for i in range(N):
        for j in range(N):
            MatrixArray[i, j] = 10 / (i + j + 1)
            
    return MatrixArray
            
def GetRightPart(N): #f
    Matrix = GetMatrix(N)
    f = np.zeros(shape = (N, 1))
    
    for i in range(N):
        for j in range(N):
            f[j, 0] += Matrix[i, j]
            
    return f

def LLDecomp(A): # A = L*L^T  
    N = len(A)
    L = np.zeros(shape = (N, N))
    
    for i in range(N):
        for k in range(i+1):
            s = sum(L[i][j] * L[k][j] for j in range(k))
            if (i != k):
                 L[i][k] = (1.0 / L[k][k] * (A[i][k] - s))
            else:
                L[i][k] = math.sqrt(A[i][i] - s)
                
    print("\nL:\n", L)
    return L

def SolveUX(U, f):# solve UX = f
    N = len(U)
    StdSolution = np.zeros(N)
    StdSolution[N - 1] = f[N - 1, 0] * 1.0 / U[N - 1][N - 1]
    for n in range(N - 2, -1, -1):
        s = 0
        for i in range(n + 1, N):
            s += U[n][i] * StdSolution[i]
        StdSolution[n] = (f[n, 0] - s) / U[n][n]
    return StdSolution

def SolveAX(A, f): # AX = f, метод Холецкого, Функция из задания
    L = LLDecomp(A)
    v = np.dot(np.linalg.inv(L), f)
    return SolveUX(np.transpose(L), v)

#########################################################

N = 6
A = GetMatrix(N)
f = GetRightPart(N)
X = SolveAX(A, f)

print("\n\nРешение Методом Холецкого\n", np.reshape(X, (len(A), 1)))

StdSolution = np.linalg.solve(A, f)
print("\n\nРешение с помощью np.linalg.solve\n", StdSolution)

print("\n Норма разницы решений", norm((np.reshape(X, (len(A), 1)) - StdSolution), 2))

