import numpy as np 

def Trapezoid(f, a: float, b: float, n: int) -> float:
    assert a < b
    
    h = (b - a) / float(n)
    ans = 0.0
    
    for i in range (0, n):
        x_i = a + i * h
        f_i = f(x_i)
        f_i1 = f(x_i + h)
        
        ans += (f_i + f_i1)
        
    
    return (ans / 2.0) * h 

a = 0.0
b = 10.0
n = 100000

# разложим интеграл по частям 
# Int_0_10[sin(x) / sqrt(x)] = [2 * sqrt(x) * sin(x)]_0_10 - # Int_0_10[2 * sqrt(x) * cos(x)]
def Integrate_By_Parts(x: float) -> float:
    return 2.0 * np.sqrt(x) * np.sin(x)


def f(x: float) -> float:
    return 2.0 * np.sqrt(x) * np.cos(x)

Ans = (Integrate_By_Parts(10.0) - Integrate_By_Parts(0.0)) - Trapezoid(f, a, b, n)


print("Ans: ", Ans)
print("сходится с правильынм ответом")
