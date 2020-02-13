# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#write functions for the shifts instead because ya girl's lazy
#A is numpy matrix array object

import numpy as np

A = np.array([[1, 2, 3]])
B = np.array([[1, 2, 3],[4, 5, 6],[7, 8, 9]])

C = np.array([[1, 2, 3],[4, 5, 6],[7, 8, 9]])
D = np.array([[1, 2], [3,4], [5, 6]])

# A = np.array([[1, 2, 3],[4, 5, 6],[7, 8, 9]])
# B = np.array([[1], [2], [3]])

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

def appendzeros(A, l, d):
    k = l
    if len(A) == len(A[0]):
        A = np.hstack((A, np.zeros((l-d, len(A))).T))
        A = np.vstack((A, np.zeros((len(A[0])-len(A), len(A[0])))))
    elif l<d:
        while d%k != 0:
            k = k+1
        A = np.vstack((A, np.zeros((abs(k-l), len(A[0])))))
    else:
        A = np.hstack((A, np.zeros((abs(l-d), len(A))).T))
    return A

def matrixmult(A, B):
    #first determine if it is lxd times dxd or dxd times dxl
    if (len(B) == len(B[0])): #this is A = lxd times B = dxd case
        #need to determine if l<d (easier case)
        l = len(A)
        k = l
        d = len(B)
        if l<d: #so l<d
            if d%l != 0: #so l does not divide d, we need to append 0s
                A = appendzeros(A, l, d)
                l = len(A)
                if l == d:
                    return matrixmultsquare(A, B)[range(k)]
        elif l>d: #so l>d, we need to append 0s to make d larger
            A = appendzeros(A, l, d)
            B = appendzeros(B, l, d)
            # print(A)
            # print(B)
            return matrixmultsquare(A, B)[:, range(d)]
        temp1 = sigma(A)
        temp2 = tau(B)
        sumM = (temp1)*(temp2[range(0, l)])
        for j in range(d//l-1):
            for i in range(l):
                temp1 = phi(temp1)
                temp2 = psi(temp2)
                sumM = sumM + (temp1)*(temp2[range(0, l)])
        return sumM
    elif (len(A) == len(A[0])): #this is A = dxd times B = dxl case
        l = len(B[0])
        k = l
        d = len(A)
        if l<d:
            if d%l != 0:
                B = appendzeros(B, d, l)
                l = len(B[0])
                if l == d:
                    return matrixmultsquare(A, B)[:, range(k)]
        elif l>d:
            A = appendzeros(A, d, l)
            B = appendzeros(B, d, l)
            return matrixmultsquare(A, B)[range(d)]
        temp1 = sigma(A)
        temp2 = tau(B)
        sumM = (temp1[:, range(0, l)])*(temp2)
        for i in range(d//l-1):
            for j in range(l):
                temp1 = phi(temp1)
                temp2 = psi(temp2)
                sumM = sumM + (temp1[:, range(0, l)])*(temp2)
        return sumM
    else:
        print("Bad matrix!")


# print(A.dot(B))
# print(matrixmult(A, B))
print(C.dot(D))
print(matrixmult(C, D))
#print(matrixmult(A,A))

#verify for n = 1 to 15, sample 1000 from each
    