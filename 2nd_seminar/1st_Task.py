import numpy as np
import math 
import matplotlib.pyplot as plt

def derivative (f, x, n, d):
    sum = 0

    for k in range (n + 1):
        sum += math.comb(n, k) * f(x + k * d) * ((-1)**(k+1))

    sum = sum / (d**n)
    return sum

def acc (f, x, deltax, deltay):

    n = 0
    u = f(0)
    while (abs(u - f(x)) >=  deltay):
        n += 1
        u += derivative(f, 0, n, deltax) * (x ** n) / (math.factorial(n))
    return n

def mistake_graph (f, x, deltax):
    n = 0
    u = f(0)

    mist = []
    iterations = 7
    arg_list = np.arange (1, iterations + 1, 1)

    for i in range (iterations):
        n += 1
        u += derivative (f, 0, n, deltax) * (x ** n) / (math.factorial(n))
        print ("N =", n, "u =", u, "|u - u*| =", abs(u - f(x))) 
        mist.append(abs(f(x) - u))

    plt.plot(arg_list, mist)
    plt.xlabel("N")
    plt.ylabel("Mistake")

    plt.show()

#####################

mistake_graph(np.sin, 0.5, 0.001)
print(acc(np.sin, 0.5, 0.001, 0.001))



