# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#write functions for the shifts instead because ya girl's lazy
#A is numpy matrix array object

import numpy as np


# B = np.array([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12], [13, 14, 15, 16]])
# A = np.array([[1, 2], [3, 4], [5, 6], [7, 8]]).T

# C = np.array([[1, 2, 3],[4, 5, 6],[7, 8, 9]])
# D = np.array([[1, 2], [3,4], [5, 6]])

# A = np.array([[1, 2, 3],[4, 5, 6],[7, 8, 9]])
# B = np.array([[1], [2], [3]])

A_test = np.array([[0.62476263, 0.79810269, 0.98722467, 0.63232729],
       [0.7962554 , 0.33747336, 0.11919426, 0.11562367],
       [0.41263048, 0.94949498, 0.78006648, 0.66001042],
       [0.66239528, 0.33640294, 0.90699305, 0.51702071]])

B_test = np.array([[0.6522363 , 0.41598821],
       [0.58782856, 0.71551526],
       [0.64035368, 0.03207762],
       [0.08844245, 0.75489649]])

def sigma(A):
    n= len(A)
    m = len(A[0])
    M = np.zeros((n, m))
    for i in range(n):
        for j in range(m):
            M[i, j] = A[i, (i+j)%m]
    return M

def tau(A):
    n= len(A)
    m = len(A[0])
    M = np.zeros((n, m))
    for i in range(n):
        for j in range(m):
            M[i, j] = A[(i+j)%n, j]
    return M

def phi(A):
    n= len(A)
    m = len(A[0])
    M = np.zeros((n, m))
    for j in range(m):
        M[:,j] = A[:,(j+1)%(m)]
    return M

def psi(A):
    n= len(A)
    m = len(A[0])
    M = np.zeros((n, m))
    for i in range(n):
        M[i,:] = A[(i+1)%n, :]
    return M

# def hprod(A, B): #hermitian product (python has this with just *, maybe pointless)
#     if len(A) == len(B):
#         if len(A[0]) == len(B[0]):
#             n = len(A)
#             m = len(A[0])
#             M = np.zeros((n, m))
#             for i in range(n):
#                 for j in range(m):
#                     M[i, j] = A[i, j] * B[i, j]
#     return M


def matrixmultsquare(A, B):
    n = len(A)
    temp1 = sigma(A)
    temp2 = tau(B)
    sumM = temp1*temp2
    for k in range(n-1):
        temp1 = phi(temp1)
        temp2 = psi(temp2)
        sumM = sumM + temp1*temp2
    return sumM

#RESUME HERE, MAKE 2 SEPARATE CASES!!!!!
def appendzeros1(A, l, d): #this is for lxd times dxd case
    k = l
    if l<d:
        if len(A) == len(A[0]):
            return A
        while d%k != 0:
            k = k+1
        A = np.vstack((A, np.zeros((abs(k-l), len(A[0])))))
    elif l>d:
        if len(A) == len(A[0]): #takes care of B
            A = np.hstack((A, np.zeros((abs(l-d), len(A))).T))
            A = np.vstack((A, np.zeros((len(A[0])-len(A), len(A[0])))))
        else:
            A = np.hstack((A, np.zeros((abs(l-d), len(A))).T))
    return A

def appendzeros2(A, l, d): #this is for dxd times lxd case
    k = l
    if l<d:
        if len(A) == len(A[0]):
            return A
        while d%k != 0:
            k = k+1
        A = np.hstack((A, np.zeros((abs(k-l), len(A))).T))
    elif l>d:
        if len(A) == len(A[0]):
            A = np.hstack((A, np.zeros((abs(l-d), len(A))).T))
            A = np.vstack((A, np.zeros((len(A[0])-len(A), len(A[0])))))
        else:
            A = np.vstack((A, np.zeros((abs(l-d), len(A[0])))))
    return A


def matrixmult(A, B):
    #first determine if it is lxd times dxd or dxd times dxl
    if (len(B) == len(B[0])): #this is A = lxd times B = dxd case
        #need to determine if l<d (easier case)
        l = len(A)
        k = l
        d = len(B)
        if (len(A) == len(A[0])):
            return matrixmultsquare(A, B)
        if l<d: #so l<d
            if d%l != 0: #so l does not divide d, we need to append 0s
                A = appendzeros1(A, l, d)
                l = len(A)
                if l == d:
                    return matrixmultsquare(A, B)[range(k)]
        elif l>d: #so l>d, we need to append 0s to make d larger
            A = appendzeros1(A, l, d)
            B = appendzeros1(B, l, d)
            #print(A)
            #print(B)
            return matrixmultsquare(A, B)[:, range(d)]
        temp1 = sigma(A)
        temp2 = tau(B)
        sumM = 0
        for j in range(d//l):
            for i in range(l):
                sumM = sumM + (temp1)*(temp2[range(0, l)])
                temp1 = phi(temp1)
                temp2 = psi(temp2)
        return sumM
    elif (len(A) == len(A[0])): #this is A = dxd times B = dxl case
        l = len(B[0])
        k = l
        d = len(A)
        if l<d:
            if d%l != 0:
                B = appendzeros2(B, l, d)
                l = len(B[0])
                if l == d:
                    return matrixmultsquare(A, B)[:, range(k)]
        elif l>d:
            A = appendzeros2(A, l, d)
            B = appendzeros2(B, l, d)
            return matrixmultsquare(A, B)[range(d)]
        temp1 = sigma(A)
        temp2 = tau(B)
        sumM = 0
        for i in range(d//l):
            for j in range(l):
                sumM = sumM + (temp1[:, range(0, l)])*(temp2)
                temp1 = phi(temp1)
                temp2 = psi(temp2)
        return sumM
    else:
        print("Bad matrix!")

# print(A_test.dot(B_test))
# print(matrixmult(A_test, B_test))

# print(A.dot(B))
# print(matrixmult(A, B))
#print(C.dot(D))
#print(matrixmult(C, D))
#print(matrixmult(A,A))

#verify for n = 1 to 15, sample 1000 from each

def verification1():
    for n in range(1, 16):
        A = np.random.rand(n, n)
        for m in range(1, 16):
            B = np.random.rand(n, m)
            T1 = A.dot(B)
            T2 = matrixmult(A, B)
            for i in range(n):
                for j in range(m):
                    if (round(T1[i][j] - T2[i][j], 13)) != 0.0:
                        return 0
    return 1

def verification2():
    for n in range(1, 16):
        B = np.random.rand(n, n)
        for m in range(1, 16):
            A = np.random.rand(m, n)
            T1 = A.dot(B)
            T2 = matrixmult(A, B)
            for i in range(m):
                for j in range(n):
                    if (round(T1[i][j] - T2[i][j], 13)) != 0.0:
                        return 0
    return 1

if verification1() == 1:
    print("Works!")
else:
    print("Sucks")

if verification2() == 1:
    print("Works also!")
else:
    print("Sucks too")

# verification2()

# A = verification()[0]
# B = verification()[1]
# val = verification()[2]

