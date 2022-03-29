### Givens&Householder

**代码详见**[复数域Givens&Householder](https://github.com/witchbunting/linear-algebra-2-codes/blob/main/%E5%A4%8D%E6%95%B0%E5%9F%9FGivens%26Householder.py)

###### Givens矩阵

* 题目要求构造矩阵幺正矩阵$U$满足$U\xi=\eta$
* 对模长相同的$\xi$和$η$同时进行Givens旋转,使得除向量第一项其他均为零，即：
    $$\xi'=G_{n-1}G_{n-2}⋯G_{2}G_{1}\xi$$$$\eta'=g_{n-1}g_{n-2}⋯g_{2}g_{1}\eta$$
由于是幺正变换，模长始终不变，最终得到的$\xi'$和$\eta'$的第一个元素模长相同，由于有复数空间的自由度，实际上两个元素相差一个相位因子$e^{i\theta}$,为除掉相位因子，可以构造一个幺正矩阵，使得仅第一个元素为相位因子，其他位置为单位阵，即：
    $$ P=\left[ \begin{matrix} e^{i\theta}&0&⋯&0&0\\ 0&1&⋯&⋯&0\\ ⋮&&\ddots&&⋮\\⋮&&&1&0\\ 0&0&⋯&0&1\end{matrix} \right]$$
其中，$e^{i\theta}=\frac{\eta'}{\xi'}$最终我们可以得到$U$的表达式：
    $$U=g_{1}^†g_{2}^†⋯g_{n-2}^†g_{n-1}^†PG_{n-1}G_{n-2}⋯G_{2}G_{1}I$$
* 对生成的$U$矩阵进行验算，并且内置在函数当中，检验：
    $$U^†U=I\\ U\xi=\eta$$
实际编程当中设定了矩阵总误差不超过$10^{-10}$才算验证通过。


###### Householder矩阵

* 求相差的相位$e^{i\theta}$
* 取与轴正交的向量$w$
* 验证与Givens类似
  