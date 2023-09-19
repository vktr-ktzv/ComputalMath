import numpy as np
import matplotlib.pyplot as plt  
### insert N
BiggestN = 30


###############3

def inverse_matrix(matrix_origin, N):
    
    # Склеиваем 2 матрицы: слева - первоначальная, справа - единичная
    
    m = np.hstack((matrix_origin, np.eye(N)))
    
    for nrow, row in enumerate(m):
        # nrow равен номеру строки
        # row содержит саму строку матрицы
        divider = row[nrow] # диагональный элемент
        # делим на диагональный элемент:
        row /= divider
        # теперь вычитаем приведённую строку из всех нижележащих строк:
        for lower_row in m[nrow+1:]:
            factor = lower_row[nrow] # элемент строки в колонке nrow
            lower_row -= factor*row # вычитаем, чтобы получить ноль в колонке nrow
    # обратный ход:
    for k in range(N - 1, 0, -1):
        for row_ in range(k - 1, -1, -1):
            if m[row_, k]:
                # 1) Все элементы выше главной диагонали делаем равными нулю
                m[row_, :] -= m[k, :] * m[row_, k]
    return m[:,N:].copy()

def Norm3(matrix):
    aig1, eigv1 =  np.linalg.eig(matrix*matrix.transpose())
    return np.sqrt(max(np.abs(aig1)))
 
def Norm1(matrix):
    sumLines = matrix.sum(axis=1, dtype=float)
    return float(max(abs(sumLines)))

def Norm2(matrix):
    sumColns = matrix.sum(axis=1, dtype=float)
    return float(max(abs(sumColns)))

def Mu (Norm1, Norm2):
    return Norm1*Norm2

########################
NormArr1 = np.zeros(BiggestN)
NormArr2 = np.zeros(BiggestN)
NormArr3 = np.zeros(BiggestN)
########################
i = 0
N = 2
while N <= BiggestN:
    A = []
    for i in range(N):
        A.append([0]*N)
        
    #fill matrix
    for i in range (N-1):
        A[i][i] = -2
        A[i][i+1] = 1

    A[N-1][N-1] = -2
    ###
    B = inverse_matrix(A, N)

    A1 = np.matrix(A)
    B1 = np.matrix(B)
    

    Nr1 = float(Norm1(A1))
    Nr2 = float(Norm1(B1))
    NormArr1[N-1] = Mu(Nr1, Nr2) 
    print(Nr1, Nr2)
    
    Nr1 = float(Norm2(A1))
    Nr2 = float(Norm2(B1))
    NormArr2[N-1] = Mu(Nr1, Nr2)  

    Nr1 = float(Norm3(A1))
    Nr2 = float(Norm3(B1))
    NormArr3[N-1] = Mu(Nr1, Nr2)   

    i += 1
    N +=1
####################

print(NormArr1)
print(NormArr2)
print(NormArr3)


plt.plot(NormArr1)
plt.plot(NormArr3)
plt.show()








