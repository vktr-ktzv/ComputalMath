import scipy
from scipy import integrate
import numpy as np
import math
import matplotlib.pyplot as plt

def MakeStepPoints(Start : float, Stop : float, N : int) -> float:
    # Возвращает массив из N точек c одинаковым интервалом
    Points = []
    Step = (Stop - Start) / (N - 1)
    for PointNum in range(N):
        Points.append(Start + PointNum * Step)
    return Points

def Getlk(t : float, k : int, ArgValues : list) -> float: #возвращает lk(t)
    # ArgValues - массив t
    n = len(ArgValues)
    lk = 1
    for j in range(n):
        if k != j:
            lk *= ((t - ArgValues[j]) / (ArgValues[k] - ArgValues[j]))  
        else: 1    
    return lk

def SimpsonIntegral(Functions : list, a : float, b : float) -> float:
    # Считает интеграл от поточечно заданной функции со значениями Functions[i] на отрезке [a, b] методом Симпсона
    Int = 0
    N = len(Functions)
    k = int(N / 2)
    h = (b - a) / N
    for i in range(1, k):
        F1 = Functions[2*i]
        F2 = Functions[2*i - 1]
        F3 = Functions[2*i - 2]
        Int += h / 3.0 * (F1 + 4 * F2 + F3)
    return Int

def GetDividedDifference(k : int, n : int, ArgValues : list, FunctionValues : list) -> float:
    # Возвращает разделенную разность 
    # ArgValues - массив значений t
    # FunctionValues - массив значений функции 
    if (k == n):
        return FunctionValues[0]
    F2 = GetDividedDifference(k + 1, n, ArgValues, FunctionValues[1:])
    F1 = GetDividedDifference(k, n - 1, ArgValues, FunctionValues[:-1])
    t2 = ArgValues[n]
    t1 = ArgValues[k]
    DivDiff = (F2 - F1) / (t2 - t1)
    return DivDiff


def GetNewtonPoly(ArgValues : list, FunctionValues : list) -> list:
    # Возвращает коэффициенты полинома Ньютона, являющиеся разделенными разностями
    n = len(ArgValues)
    DivDiffs = []
    for k in range(n):
        F = GetDividedDifference(0, k, ArgValues, FunctionValues)
        DivDiffs.append(F)
    return DivDiffs

def GetNewtonPolyValue(t : float, ArgValues : list, NewtonPol : list) -> float:
    # Возвращает значение полинома Ньютона в точке t
    Result = 0
    n = len(ArgValues)
    for k in range(n):
        Mult = 1
        for i in range(k):
            Mult *= (t - ArgValues[i])
        Result += Mult * NewtonPol[k]
    return Result

def GetNewtonValues(Args : list, ArgValues : list, FunctionValues : list) -> list:
    # Возвращает масcив значений полинома Ньютона 
    NewtonValues = []
    NewtonPol = GetNewtonPoly(ArgValues, FunctionValues)
    for Arg in Args:
        NewtonValues.append(GetNewtonPolyValue(Arg, ArgValues, NewtonPol))
    return NewtonValues


def g(X : float) -> float:
    return 1

def f(X : float) -> float:
    return math.cos(X * math.pi)

def K(X : float, S : float) -> float:
    return 0.2 / (0.04 + (X - S)**2)

def calculateTrapezoidByPoints(Functions : list, a : float, b : float) -> float:
    # Считает интеграл от поточечно заданной функции со значениями Functions[i] на отрезке [a, b] методом трапеций
    Int = 0
    N = len(Functions)
    h = (b - a) / N
    for i in range(1, N):
        F1 = Functions[i]
        F2 = Functions[i - 1]
        Int += h / 2.0 * (F1 + F2)
    return Int

def LUdecomposition(A : list) -> (list, list):
    # Разложение матрицы в произведение  A = LU
    # U - верхнетреугольная, L - нижнетреугольная
    N = len(A)
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

def solveU(U : list, F : list) -> list:
    # Решает СЛАУ UX = F
    # U - верхнетреугольная матрица
    N = len(U)
    Solution = np.zeros(N)
    Solution[N - 1] = F[N - 1, 0] * 1.0 / U[N - 1][N - 1]
    for n in range(N - 2, -1, -1):
        Sum = 0
        for i in range(n + 1, N):
            Sum += U[n][i] * Solution[i]
        Solution[n] = (F[n, 0] - Sum) / U[n][n]
    return Solution

def getSolution(Matrix : list, F : list) -> list:
    # Решает СЛАУ
    N = len(Matrix)
    L, U = LUdecomposition(Matrix)
    F = np.dot(np.linalg.inv(L), F)
    return solveU(U, F)

def getFunctionValues(Kfunc : "function", g : "function", f : "function", Lambda : float, a : float, b : float, N : int, ArrayX0 : list) -> (list, list):
    # Возвращает 
    # Nodes    - массив из N + k точек xk, полученных равным делением отрезка [a, b] на N - 1 частей и дополненный k точками
    # FunctionValues - массив значений функции u(x) в этих точках
    
    Nodes = MakeStepPoints(a, b, N)
    
    # Дополняем интервал так, чтобы масимальное из X0 попадало внутрь интервала интерполяции
    h = (b - a) / (N - 1)
    K = int((max(ArrayX0) - b) / h) + 1
    AddNodes = [b + i * h for i in range(1, K + 1)]
    Nodes += AddNodes
    
    # Ищем веса квадратурной формулы
    Weights = []
    for k in range(1, N + 1):
        Args = np.arange(a, b, 0.001)
        BaseLagranValues = [Getlk(i, k - 1, Nodes) for i in Args]
        Weights.append(SimpsonIntegral(BaseLagranValues, a, b))
        
    # Получаем матрицу СЛАУ Matrix и столбец правых частей F
    Matrix = np.zeros(shape = (N + K, N + K))
    for i in range(N + K):
        for j in range(N + K):
            if (j < N):
                Matrix[i, j] = -Lambda * Weights[j] * Kfunc(Nodes[i], Nodes[j])
    for i in range(N + K):
        Matrix[i][i] += g(Nodes[i])
        
    F = np.zeros(shape = (N + K, 1))
    for i in range(N + K):
        F[i, 0] = f(Nodes[i])
        
    # Решаем СЛАУ
    FunctionValues = getSolution(Matrix, F)
    
    return Nodes, FunctionValues

def getFunctionValue(X0 : float, ArgValues : list, FunctionValues : list) -> float:
    NumArg = 0
    while (ArgValues[NumArg] < X0):
        NumArg += 1
    return FunctionValues[NumArg]
    
        

a = -1
b = 1
Lambda = -1
# Число узлов
ArrayN = [3, 4, 5, 6]
# xi в которых хотим получить значение u(x)
ArrayX0 = [1.1, 1.25, 1.5]

plt.figure(figsize = (10, 10))
plt.title("Интерполяционный многочлен Ньютона для функции u(x)")

for N in ArrayN:
    ArgValues, FunctionValues = getFunctionValues(K, g, f, Lambda, a, b, N, ArrayX0)
    plt.scatter(ArgValues, FunctionValues, marker = "^")
    Args = np.arange(-1, max(ArgValues) + 0.01, 0.01)
    NewtonVals = GetNewtonValues(Args, ArgValues, FunctionValues)
    print("\n\nЗначения в точках с помощью построения полинома по n =", N, "узлам:\n")
    for X0 in ArrayX0:
        print("X0 =", X0, ", u(X0) =", format(getFunctionValue(X0, Args, NewtonVals), '.4f'), "\n")
    
    NameGraph = "Ньютон n =" + str(N)
    plt.plot(Args, NewtonVals, label = NameGraph)

plt.grid()
plt.legend()
plt.show()


