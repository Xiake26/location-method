# encoding=utf-8
# cython: language_level=3
from libc.math cimport sqrt,isnan,abs
def f1_cython(double[:, :] u, double h, double[:, :] S, double dt, int node1, int node2, int K):
    cdef int i, j, kk, order
    cdef double uymin, uxmin, ubar, diff
    cdef double[:, :] u2 = u
    cdef int[4][3] Order_j = [[0, node1, 1], [node1-1, -1, -1], [0, node1, 1], [node1-1, -1, -1]]
    cdef int[4][3] Order_i = [[node2-1, -1, -1], [0, node2, 1], [0, node2, 1], [node2-1, -1, -1]]
    for kk in range(K):
        for order in range(4):
            for i in range(Order_i[order][0], Order_i[order][1], Order_i[order][2]):
                for j in range(Order_j[order][0], Order_j[order][1], Order_j[order][2]):
                    if not isnan(u[i, j]):
                        if i > 0 and i < node2-1:
                            if u[i-1, j] < u[i+1, j]:
                                uymin = u[i-1, j]
                            else:
                                uymin = u[i+1, j]
                        elif i == 0:
                            uymin = u[1, j]
                        else:
                            uymin = u[node2-2, j]
                        if j > 0 and j < node1-1:
                            if u[i, j-1] < u[i, j+1]:
                                uxmin = u[i, j-1]
                            else:
                                uxmin = u[i, j+1]
                        elif j == 0:
                            uxmin = u[i, 1]
                        else:
                            uxmin = u[i, node1-2]
                        if abs(uxmin - uymin) >= h * S[i, j]:
                            ubar = min(uxmin, uymin) + h * S[i, j]
                        else:
                            ubar = (uxmin + uymin + sqrt(2 * S[i, j] ** 2 * h ** 2 - (uxmin - uymin) ** 2)) / 2
                        u[i, j] = min(u[i, j], ubar)

        if kk > 0:
            diff = 0
            for i from 0 <= i < node2 by 1:
                for j from 0 <= j < node1 by 1:
                    diff += (u[i, j] - u2[i, j]) ** 2
            if sqrt(diff) < dt:
                break
            else:
                u2[:, :] = u

    return u