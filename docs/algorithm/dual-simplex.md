# 对偶修正单纯型法

简体中文 | [English](./dual-simplex_en.md) | [项目主页](../../README.md) | [算法总览](./README.md)

## 摘要

对偶单纯型法保持检验数与非基端点的符号相容，逐步消除基本变量的界违反。本文给出有界变量形式下的定价、比值检验与翻界公式。

## 1. 基关系与初始化

对 $Ax=b,\ell\le x\le u$，非基变量固定在选定端点，基本值满足

$$
Bx_\mathcal B=b-A_\mathcal N x_\mathcal N,\qquad
B^T\pi=c_\mathcal B,\qquad r=c-A^T\pi.
$$

对偶可行性要求：下界非基变量 $r_j\ge0$，上界非基变量 $r_j\le0$，自由非基变量 $r_j=0$，固定变量无符号限制。基本检验数恒为零。若再有 $\ell_\mathcal B\le x_\mathcal B\le u_\mathcal B$，即满足最优性条件。

在行变量形式中，逻辑列组成 $-I$，可提供非奇异初始基；但该事实本身不保证对偶可行。需通过端点选择、基交换或辅助目标建立对偶可行状态。辅助目标改变过的成本在最终最优性判断前必须撤销。

## 2. 离基行与最陡边定价

令基本位置 $i$ 的界违反为

$$
v_i=\max\{\ell_{\mathcal B_i}-x_{\mathcal B_i},\,x_{\mathcal B_i}-u_{\mathcal B_i},\,0\}.
$$

对偶最陡边权重为 $\gamma_i=\|B^{-T}e_i\|_2^2$，选择 $v_i^2/\gamma_i$ 最大的位置 $p$。令

$$
\rho=B^{-T}e_p,\qquad \alpha_j=\rho^Ta_j.
$$

若 $x_{\mathcal B_p}<\ell_{\mathcal B_p}$，取修复符号 $\kappa=1$；若超过上界，取 $\kappa=-1$。下一节以此统一两种情形。

## 3. 对偶步长

取乘子路径及对应检验数

$$
\pi(\theta)=\pi-\kappa\theta\rho,\qquad
r_j(\theta)=r_j+\kappa\theta\alpha_j.
$$

对于非固定的端点非基变量，令 $s_j=1$ 表示下界，$s_j=-1$ 表示上界。对偶可行性为 $s_jr_j\ge0$。只有

$$
\beta_j=-\kappa s_j\alpha_j>0
$$

的变量会随 $\theta$ 增加失去对偶可行性，其断点为

$$
\theta_j=\frac{s_jr_j}{\beta_j}.
$$

不作翻界时，取最小断点并令相应变量进入。自由非基变量要求检验数始终为零，故当 $\alpha_j\ne0$ 时在零步长处形成阻挡；不能赋予它不存在的有限端点。

## 4. 原始步与基交换

设进入变量为 $q$，$d=B^{-1}a_q$，离基变量所到端点为 $h_p$。令

$$
\Delta x_q=\frac{x_{\mathcal B_p}-h_p}{d_p},\qquad
\Delta x_\mathcal B=-d\,\Delta x_q.
$$

由 $d_p=\alpha_q$ 可知，上述符号条件保证 $\Delta x_q$ 沿进入变量允许的端点方向。新基以 $a_q$ 替换第 $p$ 列。对偶步长使 $r_q$ 变为零，离基变量在相应端点满足所需的检验数符号。

新基本变量允许暂时违反自身界，因为此阶段的不变量是对偶可行性。要求每次换基同时获得全部基本变量可行，会排除合法的对偶单纯型步。

## 5. 有界变量翻界

若候选变量两端有限，可在检验数过零处将其从一个端点移到另一端点。翻界引起

$$
\Delta x_\mathcal B=-B^{-1}a_j\Delta x_j.
$$

非基端点符号反转后，可继续沿对偶路径推进。多个翻界对基本值的总影响为

$$
\Delta x_\mathcal B=-B^{-1}\sum_{j\in F}a_j\Delta x_j.
$$

断点应按对偶步长的单调次序处理。翻界与最终交换必须构成一致的原始—对偶更新，不能仅移动乘子而忽略基本等式 $B^T\pi=c_\mathcal B$。

若某离基变量仍违反界，且所有非基变量在允许方向上均无法改善该违反，则可由 $\rho^TAx=\rho^Tb$ 与变量盒区间导出矛盾，构成不可行分离证明。浮点候选筛选为空不足以单独确立该结论，必须验证完整定价行与分离不等式。

## 6. 成对单元素列的结构对偶化

考虑

$$
\min\ c_z^Tz+(c^+)^Tp+(c^-)^Tq,
\quad A_zz+D^+p+D^-q=b,\quad z,p,q\ge0,
$$

其中 $D^+,D^-$ 为对角矩阵，$d_i^+>0$、$d_i^-<0$，且 $c_i^+,c_i^->0$。其对偶为

$$
\max_\pi b^T\pi,\qquad
A_z^T\pi\le c_z,\qquad
\frac{c_i^-}{d_i^-}\le\pi_i\le\frac{c_i^+}{d_i^+}.
$$

成对单元素列因此转化为对偶变量的有限双界。由对偶问题相应约束的最优乘子恢复 $z$ 后，令 $r=b-A_zz$，再取

$$
p_i=\frac{\max(r_i,0)}{d_i^+},\qquad
q_i=\frac{\min(r_i,0)}{d_i^-}.
$$

这给出原问题可行解；在精确最优乘子与互补条件成立时，它亦为原问题最优解。该化简依赖上述结构条件，并非一般等式模型都具有有限对偶盒界。
