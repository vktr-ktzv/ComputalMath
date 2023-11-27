import numpy as np 

Eps_tail = 1e-4 # Tail cotoff 
Eps = 1e-4 # Integral calculation precision

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

# |cos(x)| / (2 + x^2) <= 1 / (2 + x^2)
# Int[1 / (2 + x^2)] = arctan(x / sqrt(2)) / sqrt(2)
# Tail = arctan(x / sqrt(2)) / sqrt(2) | Cutoff -> +inf
# Eps_tail = [np.pi / 2 - arctan(Cutoff / sqrt(2))] / sqrt(2)
# np.pi / 2 - sqrt(2) * Eps_tail = arctan(Cutoff / sqrt(2))
Cutoff = np.sqrt(2.0) * np.tan(np.pi / 2.0 - np.sqrt(2.0) * Eps_tail)

a = 0.0
b = Cutoff

def f(x: float) -> float: 
    return np.cos(x) / (2.0 + x ** 2.0)

# получено из вольфрама 
f_second_der_max_abs = 1.0

# Eps >= (b - a) / 8.0 * f_second_der_max_abs * ((b - a) / n) ** 2.0 
# 8.0 * Eps / ((b - a) ** 3.0 * f_second_der_max_abs) >= n ** -2.0
n2 = int((8.0 * Eps / ((b - a) ** 3.0 * f_second_der_max_abs)) ** (-0.5))
print(n2, " - вот так много итераций получается")
h2 = (b - a) / float(n2)

# точное значение из вольфрама
TrueAns= np.exp(-np.sqrt(2.0)) * np.pi / (2.0 * np.sqrt(2.0))



TrapAns = Trapezoid(f, 0, Cutoff, int(n2))

print("h :", h2)
print("Trapezoid: ", TrapAns)
print ("TrueAns:", TrueAns)
print("|TrueAns - Trapezoid|:", abs(TrapAns - TrueAns))