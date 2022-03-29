import numpy as np
import random

def rand_vector(n=3):  # 生成两个模长相同（norm=1）的复向量xi，yita
    a1, a2, b1, b2 = np.random.random(n), np.random.random(n), np.random.random(n), np.random.random(n)
    xi = a1 + 1j * a2
    yita = b1 + 1j * b2
    xi = xi.reshape((n, 1)) / np.linalg.norm(xi)
    yita = yita.reshape((n, 1)) / np.linalg.norm(yita)
    return xi, yita


def givens(n: int, xi, yita):  # 生成相应的givens矩阵
    U = np.identity(n, dtype="complex_")
    xi0, yita0 = xi, yita  # 储存初始xi和yita的值，用于验算
    G1 = []
    G2 = []
    for i in range(1, n):
        g1 = np.identity(n, dtype="complex_")
        g2 = np.identity(n, dtype="complex_")
        ta1 = np.abs(-xi[i] / xi[0])
        ta2 = np.abs(-yita[i] / yita[0])
        ph1 = -xi[i] / xi[0] / ta1
        ph2 = -yita[i] / yita[0] / ta2
        c1, c2 = 1 / (1 + ta1 ** 2) ** (1 / 2), 1 / (1 + ta2 ** 2) ** (1 / 2)
        s1, s2 = c1 * ta1 * ph1, c2 * ta2 * ph2
        g1[0, 0], g1[i, 0], g1[0, i], g1[i, i] = c1, s1, -s1.conjugate(), c1  # 获取相应的旋转元，生成每一步旋转的矩阵
        g2[0, 0], g2[i, 0], g2[0, i], g2[i, i] = c2, s2, -s2.conjugate(), c2
        xi = np.dot(g1, xi)  # 同时对xi,yita进行旋转
        yita = np.dot(g2, yita)
        G1.append(g1)  # 保存每一步旋转xi的矩阵和对yita旋转的矩阵复共轭
        G2.append(g2.T.conjugate())
    for j in range(len(G1)):  # 取单位阵U，将保存的g1作用于U，得到xi的变换矩阵
        U = np.dot(G1[j], U)
    phase = np.identity(n, dtype="complex_")  # 两向量旋转后在复空间仍有自由度，求复空间的相位因子进行旋转，使得转后的xi=yita
    phase[0, 0] = yita[0, 0] / xi[0, 0]
    U = np.dot(phase, U)  # 相位因子作用于矩阵U

    for k in range(len(G2) - 1, -1, -1):  # 将yita旋转矩阵的复共轭倒序作用于U，得到最终的U
        U = np.dot(G2[k], U)
    if np.linalg.norm(np.dot(U.T.conjugate(), U) - np.identity(n)) <= 10 ** (-10) and np.linalg.norm(
            np.dot(U, xi0) - yita0) <= 10 ** (-10):  # 检验该矩阵U为幺正矩阵，并且是作用于xi得到yita的变换矩阵
        print("该矩阵是幺正矩阵,并且作用于xi得到yita")
        return U
    else:
        print("该矩阵不为幺正矩阵，或不为xi->yita的变换矩阵")


def householder(n: int, xi, yita):  # householder矩阵
    phase = np.dot(xi.T.conjugate(), yita) / np.linalg.norm(np.dot(xi.T.conjugate(), yita))  # 得到相位因子e^(i\theta)
    w = (phase * xi - yita) / np.linalg.norm(phase * xi - yita)  # 得到W
    U = (np.identity(n) - 2 * w.T.conjugate() * w) * phase  # 得到变换矩阵U
    if np.linalg.norm(np.dot(U.T.conjugate(), U) - np.identity(n)) <= 10 ** (-10) and np.linalg.norm(
            np.dot(U, xi) - yita) <= 10 ** (-10):  # 检验该矩阵U为幺正矩阵，并且是作用于xi得到yita的变换矩阵
        print("该矩阵是幺正矩阵,并且作用于xi得到yita")
        return U
    else:
        print("该矩阵不为幺正矩阵，或不为xi->yita的变换矩阵")


if __name__ == '__main__':
    for i in range(3):
        n=random.randint(2,100)
        print("随机选择的n为：",n)
        n=2
        xi, yita = rand_vector(n)
        print("xi的值为:",xi,"\n yita的值为:",yita)
        U_g = givens(n, xi, yita)
        print("得到的Givens矩阵为：",U_g)
        U_h = householder(n, xi, yita)
        print("得到的householder矩阵为：",U_h)
