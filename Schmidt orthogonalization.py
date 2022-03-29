import numpy as np


def orthogonal(n):  # Schmidt orthogonalization
    mat = np.zeros((n, n))
    matric = np.random.rand(n ** 2).reshape((n, n))
    mat[0, :] = matric[0, :]
    for i in range(1, n):
        l1 = matric[i, :]
        ln = l1
        for j in range(i):  # 获得正交基
            v = mat[j, :]
            ln -= np.dot(np.dot(l1, v) / np.dot(v, v), v)
        mat[i, :] = ln
    for _ in range(n):  # 获得单位基
        mat[_, :] = mat[_, :] / np.linalg.norm(mat[_, :])
    if (np.dot(mat.T, mat).astype(int) - np.eye(n) == 0).all():  # 检查是否为正交矩阵
        return mat
    else:
        return orthogonal(n)


def projection(vec, matric, n):
    try:
        x = np.zeros(n)
        for i in range(n):
            x[i] = np.dot(matric[i, :], vec)
        return x


    except Exception as e:
        return e


n = 5
W = orthogonal(n)
vec = np.random.rand(n)
proj = projection(vec, W, n)
print("生成的一组单位正交基W为：", W, "\n生成的随机向量的投影向量为：", proj)
