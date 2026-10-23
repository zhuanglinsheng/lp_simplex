# Primal–Dual Hybrid Gradient and Its Connection to Active-Set Methods

[简体中文](./pdhg.md) | English | [Project home](../../README_en.md) | [Algorithm overview](./README_en.md)

## Abstract

Primal–dual hybrid gradient solves the saddle-point formulation of a linear program through box projections and equality-residual updates. Iterates may be infeasible and serve as approximate solutions. Feasibility recovery and active-set descent connect these approximations to feasible optimal points at the required accuracy.

## 1. Saddle-point formulation and iterations

For $\min_{x\in C}c^Tx$ subject to $Ax=b$, consider

$$
\min_x\max_y\ c^Tx+\delta_C(x)+y^T(Ax-b),
$$

Here $\delta_C$ is the indicator function of $C$. With positive diagonal step matrices $T,\Sigma$, a primal-first PDHG iteration is

$$
\begin{aligned}
x^{k+1}&=\Pi_C\bigl(x^k-T(c+A^Ty^k)\bigr),\\
\bar x^{k+1}&=2x^{k+1}-x^k,\\
y^{k+1}&=y^k+\Sigma(A\bar x^{k+1}-b).
\end{aligned}
$$

Projection onto a box in a diagonal metric is coordinatewise: $[\Pi_C(z)]_j=\min\{u_j,\max\{\ell_j,z_j\}\}$. An objective constant does not affect the iteration.

Under standard assumptions including existence of a saddle point and fixed step sizes satisfying

$$
\|\Sigma^{1/2}AT^{1/2}\|_2<1
$$

convergence and ergodic saddle-gap estimates on bounded test sets can be established. A global gap over an unrestricted dual domain or an unbounded primal domain need not be finite, so it cannot be stated unconditionally as a finite error bound over the full domain.

## 2. Diagonal preconditioning and scale balancing

For nonzero rows and columns, one may take

$$
T_{jj}=\frac{\eta}{\sum_i|a_{ij}|},\qquad
\Sigma_{ii}=\frac{\eta}{\sum_j|a_{ij}|},\qquad 0<\eta<1.
$$

A weighted Cauchy–Schwarz inequality gives $\|\Sigma^{1/2}AT^{1/2}\|_2\le\eta$. Zero rows and columns may receive arbitrary positive step sizes; their feasibility or independent objective contributions require separate consideration.

The transformation $T\leftarrow\omega T,\Sigma\leftarrow\Sigma/\omega$ preserves this operator norm and adjusts the relative scales of primal and dual displacement. If costs are scaled to $c/s$, an iteration multiplier $\widehat y$ corresponds to the original-objective multiplier $y=s\widehat y$. Fixed-step convergence theorems do not apply unchanged to adaptive metrics or restarts without additional conditions.

## 3. Averaging, residuals, and restarts

Weighted average iterates are

$$
\widetilde x^K=\frac{\sum_{k=1}^K w_kx^k}{\sum_{k=1}^K w_k},\qquad
\widetilde y^K=\frac{\sum_{k=1}^K w_ky^k}{\sum_{k=1}^K w_k},\quad w_k>0.
$$

The fixed-point residual

$$
R(x,y)^2=\|\Sigma^{1/2}(Ax-b)\|_2^2
 +\left\|T^{-1/2}\bigl[x-\Pi_C(x-T(c+A^Ty))\bigr]\right\|_2^2
$$

can compare current and averaged pairs. The condition $R=0$ is equivalent to primal feasibility and the box KKT conditions. Averaging can reduce oscillations. A restart selects a current or averaged pair as a new starting point and resets the averaging history. The meaning of a small nonzero residual depends on variable scales and step sizes; it does not directly replace accuracy checks in original units.

## 4. Box dual bounds and duality gaps

Let $g=c+A^Ty$. A dual lower bound is

$$
q(y)=c_0-b^Ty+\sum_j\varphi_j(g_j),\qquad
\varphi_j(g_j)=\inf_{\ell_j\le t\le u_j}g_jt.
$$

When the infimum is finite, use $g_j\ell_j$ for $g_j>0$, $g_ju_j$ for $g_j<0$, and zero for $g_j=0$. If the required improving endpoint is infinite, $\varphi_j=-\infty$. Omitting this term would produce a false finite lower bound.

At a primal feasible point, the gap is

$$
c^Tx+c_0-q(y)\ge0.
$$

A zero gap proves optimality. Floating-point termination should control primal residuals, box sign residuals, and the original-objective gap together. An approximately infeasible point may have an apparently small objective difference that lacks the same optimality interpretation.

## 5. Combining feasibility recovery and BDA

The hybrid method comprises three mathematical tasks: saddle-point iterations produce an approximate primal–dual pair; bounded least squares restores primal feasibility; and a tangent-cone or deficient-basis method continues descent and constructs multipliers. A recovered point must satisfy the original constraints and the agreed objective ceiling. Multipliers passed between phases must also be restored to the original cost scale.

An approximately feasible point need not have linearly independent positive support. It cannot be assumed to define a feasible deficient basis for classical BDA without support compression. A tangent-cone method can instead express all allowed motions on the native box before using independent active columns to solve a direction subproblem.

Residual decrease, phase switching, and computational budgets in this combination do not change the definition of final optimality. Termination of one auxiliary problem does not automatically terminate the linear program. The final conclusion still rests on original-model primal feasibility, dual feasibility, and complementarity.
