import math
import numpy as np 
import matplotlib.pyplot as plt

def function(t, u):
    return 1 / (1 - u * math.cos(t))

def CustomMethod(f, u0, t0, t1, h = 1e-3):
    t_list = [t0 + i * h for i in range(int((t1-t0) / h))]
    u_list = [u0, u0 + f(t0, u0) * h]

    for i in range(len(t_list) - 2):
        i_temp = i + 2
        u_list.append(u_list[i_temp - 2] + h * (f(t_list[i_temp - 1], u_list[i_temp - 1]) + f(t_list[i_temp - 2], u_list[i_temp - 2])))

    return t_list, u_list

def count_R(u_list : list):
    R_num = len(u_list) - 1
    R_list = []

    for i in range(R_num):
        i_temp = i + 1
        R_list.append(u_list[i_temp] / u_list[i_temp - 1])

    return R_list

def sustainabilityCheck(R_list : list):
    for R in R_list:
        if abs(R) > 1:
            print("Не устойчивый")
            return
    
    print("Устойчивый")

def create_plot(x, y, title, x_label = 'x', y_label = 'y'):
    plt.figure(figsize = [12, 5])
    plt.plot(x, y, "r")
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.grid()
    plt.show()

def main():
    u0 = -1.5
    t0 = -1.5
    t1 = 0.55

    t_list, u_list = CustomMethod(function, u0, t0, t1, 1e-6)
    create_plot(t_list, u_list, 'u(t)', 't', 'u')

    

    R_list = count_R(u_list)
    sustainabilityCheck(R_list)
    
    n_list = [n for n in range(len(R_list))]
    create_plot(n_list, R_list, 'R(n)', 'n', 'R')
    
    #print(R_list)
    

    return 0

if (__name__ == '__main__'):
    main()