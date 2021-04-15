
# # Чижов Пётр Сергеевич; МКР

from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve, inv
import numpy as np
from numpy.linalg import norm
from math import pi, sin, cos, sqrt
import matplotlib.pyplot as plt

def f(x, y):
    return cos(pi * x) * cos(pi * y)

def laplas_f(x, y, dx, dy):
    return (dx + dy) * pi * pi * cos(pi * x) * cos(pi * y)

def interpolation(mass, i, j, x, y, h):
    temp1 = mass[(i + 1) * N + j] * (x - i*h) / h + mass[i * N + j] * ((i + 1) * h - x) / h 
    temp2 = mass[(i + 1) * N + j + 1] * (x - i*h) / h + mass[i * N + j + 1] * ((i + 1) * h - x) / h
    return temp1 * ((j + 1) * h - y) / h + temp2 * (y - j * h) / h

def cubature(mass, i, j, h):
    baricentric = np.array([[0.124949503233232, 0.437525248383384, 0.437525248383384],
                            [0.437525248383384, 0.124949503233232, 0.437525248383384],
                            [0.437525248383384, 0.437525248383384, 0.124949503233232],
                            [0.797112651860071, 0.165409927389841, 0.037477420750088],
                            [0.797112651860071, 0.037477420750088, 0.165409927389841],
                            [0.165409927389841, 0.797112651860071, 0.037477420750088],
                            [0.165409927389841, 0.037477420750088, 0.797112651860071],
                            [0.037477420750088, 0.165409927389841, 0.797112651860071],
                            [0.037477420750088, 0.797112651860071, 0.165409927389841]])
    #print (baricentric.shape)
    baricentric = baricentric.T
    #print(baricentric.shape)
    weight1 = 0.205950504760887
    weight2 = 0.063691414286223
    
    
    #coordinate transition matrices for first triangle
    conv_matr1 = np.array([[i * h, i * h, (i + 1) * h],
                          [j * h, (j + 1) * h, j * h],
                          [1, 1, 1]])
    
    #coordinate transition matrices for second triangle
    conv_matr2 = np.array([[i * h, (i + 1) * h, (i + 1) * h],
                          [(j + 1) * h, j * h, (j + 1) * h],
                          [1, 1, 1]])
    
    
    temp = np.matmul(conv_matr1, baricentric)
    
    temp = temp.T
    #if (i == 50 and j == 51):
    #    print(temp)
    
    summ1 = 0.0
    summ2 = 0.0
    iteratr = 0
    
    for vec in temp:
        iteratr += 1
        if (iteratr < 4):
            summ1 += weight1 * (interpolation(mass, i, j, vec[0], vec[1], h) - f(vec[0], vec[1]))**2
        else:
            summ1 += weight2 * (interpolation(mass, i, j, vec[0], vec[1], h) - f(vec[0], vec[1]))**2

            
    temp = np.matmul(conv_matr2, baricentric)
    iteratr = 0
    
    temp = temp.T
    
    for vec in temp:
        iteratr += 1
        if (iteratr < 4):
            summ2 += weight1 * (interpolation(mass, i, j, vec[0], vec[1], h) - f(vec[0], vec[1]))**2
        else:
            summ2 += weight2 * (interpolation(mass, i, j, vec[0], vec[1], h) - f(vec[0], vec[1]))**2 
            
    return (summ1 + summ2) * h**2       
            
N = 129            
#N = 65
#N = 33
dx = 1.0
#dy = 1.0
#dy = 10.0   #eps
dy = 100.0
h = 1.0 / (N - 1)
flag = 2 # approximation order of the Neumann boundary condition (1 or 2)
row = []
col = []
data = []

b = np.zeros(N * N)
analitic_sol = np.zeros(N * N)

def add_elem(r, c, d):
    row.append(r)
    col.append(c)
    data.append(d)
    

for i in range(0, N):
    for j in range(0, N):
        analitic_sol[i * N + j] = f(i * h, j * h)
        
        # init Neumann border crossing
        if (i == 0 and j == 0):
            if flag == 1:
                # order 1 of the Neumann boundary condition
                add_elem(i * N + j, i * N + j, (dx + dy) / h**2) # -2
                add_elem(i * N + j, i * N + j + 1, -dy / h**2)   #  1
                add_elem(i * N + j, (i + 1) * N + j, -dx / h**2) #  1
                #b[i * N + j] = 0.0  #analitic solution for g_N(0, 0)
                b[i * N + j] = -2 * pi * cos(pi * i * h) * sin(pi * j * h) / h + 0.5 * laplas_f(i * h, j * h, dx, dy)
            else:
                # order 2 of the Neumann boundary condition
                add_elem(i * N + j, i * N + j, (dx + dy) / h**2)
                add_elem(i * N + j, i * N + j + 1, -dy / h**2)
                add_elem(i * N + j, (i + 1) * N + j, -dx / h**2)
                b[i * N + j] = 0.5 * laplas_f(i * h, j * h, dx, dy)
            continue
            
        # right Dirichlet border    
        if i == N - 1: 
            add_elem(i * N + j, i * N + j, 1.0)
            b[i * N + j] = f(i * h, j * h)
            continue
            
        # up Dirichlet border    
        if j == N - 1:  
            add_elem(i * N + j, i * N + j, 1.0)
            b[i * N + j] = f(i * h, j * h)
            continue
        
        #left Neumann border
        if i == 0: 
            if flag == 1:
                add_elem(i * N + j, i * N + j, -dx / h)
                add_elem(i * N + j, (i + 1) * N + j, dx / h)
                #b[i * N + j] = 0.0  #analitic solution for Neumann funct
                b[i * N + j] = -pi * sin(pi * i * h) * cos(pi * j * h)
            else:
                #order 2 of the Neumann boundary condition
                add_elem(i * N + j, i * N + j - 1, -dy / (2 * h**2))
                add_elem(i * N + j, i * N + j, (dx + dy) / h**2)
                add_elem(i * N + j, i * N + j + 1, -dy / (2 * h**2))
                add_elem(i * N + j, (i + 1) * N + j, -dx / h**2)
                b[i * N + j] = 0.5 * laplas_f(i * h, j * h, dx, dy)
            continue
        #####
        #down Neumann border
        if j == 0: 
            if flag == 1:
                add_elem(i * N + j, i * N + j, -dy / h)
                add_elem(i * N + j, i * N + j + 1, dy / h)
                b[i * N + j] = -pi * cos(pi * i * h) * sin(pi * j * h) #equals null 
            else:
                #order 2 of the Neumann boundary condition
                add_elem(i * N + j, (i - 1) * N + j, -dx / (2 * h**2))
                add_elem(i * N + j, i * N + j, (dx + dy) / h**2)
                add_elem(i * N + j, i * N + j + 1, -dy / h**2)
                add_elem(i * N + j, (i + 1) * N + j, -dx / (2 * h**2))
                b[i * N + j] = 0.5 * laplas_f(i * h, j * h, dx, dy)
            
            continue
        
        add_elem(i * N + j, (i - 1) * N + j, -dx / h**2)
        add_elem(i * N + j, i * N + j - 1, -dy / h**2)
        
        add_elem(i * N + j, i * N + j, 2 *(dx + dy) / h**2)
        add_elem(i * N + j, i * N + j + 1, -dy / h**2)
        
        add_elem(i * N + j, (i + 1) * N + j, -dx / h**2)
        
        b[i * N + j] = laplas_f(j * h, i * h, dx, dy)
        


#get the matrix in csr format (input: a[row_ind[k], col_ind[k]] = data[k])
matr = csr_matrix((data, (row, col)), shape=(N * N, N * N)).astype('float')

x = spsolve(matr, b, permc_spec='NATURAL')


#problem s flag = 1
#print('L2norm: ', norm(x - analitic_sol, ord = 2))
print('Ch_norm: ', norm(x - analitic_sol, ord = np.inf))

L2_norm = 0.0
for i in range(N - 1):
    for j in range(N - 1):
        L2_norm += cubature(x, i, j, h) 

print('L2_norm: ', sqrt(L2_norm))        
##########


mtx = []
for i in range(N):
    temp = []
    for j in range(N):
        temp.append(np.fabs(x[i * N + j] - analitic_sol[i * N + j]))
    mtx.append(temp)

    
string = 'N_' + str(N - 1) + 'Eps_' + str(dy) + '.png'
fig, ax = plt.subplots()
    
#ax.pcolormesh(mtx)    
#ax.imshow(mt)

plt.imshow(mtx, origin="lower", interpolation='nearest')
plt.colorbar()
fig.set_size_inches(5, 5)
fig.savefig(string, dpi=150)

plt.show()



mtx = []
for i in range(N):
    temp = []
    for j in range(N):
        temp.append(analitic_sol[i * N + j])
    mtx.append(temp)

    
string = 'analitic_sol' + str(N - 1) + '.png'
fig, ax = plt.subplots()
    
#ax.pcolormesh(mtx)    
#ax.imshow(mt)

plt.imshow(mtx, origin="lower", interpolation='nearest')
plt.colorbar()
fig.set_size_inches(5, 5)
fig.savefig(string, dpi=150)

plt.show()       


mtx = []
for i in range(N):
    temp = []
    for j in range(N):
        temp.append(x[i * N + j])
    mtx.append(temp)

    
string = 'approximate_sol' + str(N - 1) + '.png'
fig, ax = plt.subplots()
    
#ax.pcolormesh(mtx)    
#ax.imshow(mt)

plt.imshow(mtx, origin="lower", interpolation='nearest')
plt.colorbar()
fig.set_size_inches(5, 5)
fig.savefig(string, dpi=150)

plt.show()
