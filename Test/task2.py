import numpy as np
import math
import matplotlib.pyplot as plt

def getl(t : float, k : int, ArgValues : list) -> float:
    # Находит базисную функцию Лагранжа lk(t)
    # ArgValues - массив из n значений аргумента [t0, ... tn]
    n = len(ArgValues)
    lk = 1
    for j in range(n):
        Denom = ArgValues[k] - ArgValues[j]
        lk *= ((t - ArgValues[j]) / Denom) if k != j else 1
    return lk

def getLagrangePolinom(t : float, ArgValues : list, FunctionValues : list) -> float:
    # Возвращает значение полинома Лагранжа в точке t
    n = len(ArgValues)
    Value = 0
    for k in range(n):
        Value += (getl(t, k, ArgValues) * FunctionValues[k])
    return Value

def getLagrangeValues(Args : list, ArgValues : list, FunctionValues : list) -> list:
    # Возвращает массив значений полинома Лагранжа для каждого значения аргумента из Args
    # ArgValues      - массив из n значений аргумента t -- [t0, ... tn]
    # FunctionValues - массив из n значений функции в точках [t0, ... tn]
    Values = []
    for Arg in Args:
        Values.append(getLagrangePolinom(Arg, ArgValues, FunctionValues))
    return Values

LagrangeArgValues      = [6, 4.8, 2.85, 2.3, 1.92, -2.5]
LagrangeFunctionValues = [-2, -1.58, -0.75, -0.33, 0.08, 0.5]
X0 = 0

plt.figure(figsize = (10, 10))
plt.title("Интерполяционный многочлен Лагранжа")

Args = np.arange(-2.5, 6, 0.01)
LagrangeVals = getLagrangeValues(Args, LagrangeArgValues, LagrangeFunctionValues)

for k in range(6):
    print("k =", k, ", lk =", getl(X0, k, LagrangeArgValues))
plt.plot(Args, LagrangeVals, 'y')
plt.scatter(LagrangeArgValues, LagrangeFunctionValues, marker = "^")
    
plt.grid()
plt.show()

