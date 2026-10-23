# Linear Programming Algorithms

[简体中文](./README.md) | English | [Project home](../../README_en.md) | [Documentation](../README_en.md)

This collection studies equivalent transformations, feasibility construction, optimality conditions, and iterative methods for continuous linear programming. A minimization convention is used throughout. Except in discussions of numerical error, propositions are stated in exact arithmetic. Finite termination and convergence claims apply only under the conditions explicitly stated.

## 1. Problem and notation

The general model is

$$
\min_x\ c^Tx+c_0,\qquad \ell^r\le Ax\le u^r,\qquad \ell\le x\le u.
$$

Infinite endpoints indicate the absence of the corresponding constraint; equal finite row bounds represent an equality. After introducing row variables, the model can be expressed using equalities and box constraints:

$$
\tag{P}\min_x\ c^Tx+c_0,\qquad Ax=b,\quad x\in C=[\ell,u].
$$

Hereafter, $A\in\mathbb R^{m\times n}$, $b\in\mathbb R^m$, and $c\in\mathbb R^n$ refer to the current equivalent model. The symbol $a_j$ denotes column $j$, $A_J$ a column submatrix, and $\mathcal R(B)$ the column space. Infinite endpoints are allowed, but $C$ must be nonempty and no coordinate may be forced to take an infinite value.

A full basis $B=A_{\mathcal B}$ is a nonsingular square matrix of order $m$. A deficient basis has $k\le m$ linearly independent columns. The set $\mathcal N$ indexes nonbasic variables, and $x_\mathcal B$ denotes the basic coordinates. The two notions of basis remain distinct even when the same symbol is used.

## 2. Primal–dual conditions

Define the Lagrangian by

$$
L(x,y)=c^Tx+c_0+y^T(Ax-b),\qquad g=c+A^Ty.
$$

The dual function is

$$
q(y)=c_0-b^Ty+\sum_j\inf_{\ell_j\le t\le u_j}g_jt.
$$

If a finite primal optimum is attained, strong duality for linear programming guarantees an optimal multiplier. A primal feasible point $x$ and a multiplier $y$ form an optimal pair when

$$
\tag{KKT}
Ax=b,\quad x\in C,\qquad
\begin{cases}
g_j\ge0,&x_j=\ell_j<u_j,\\
g_j\le0,&x_j=u_j>\ell_j,\\
g_j=0,&\ell_j<x_j<u_j.
\end{cases}
$$

There is no sign restriction on $g_j$ at a fixed coordinate. For one-sided and free variables, these conditions are interpreted using only the endpoints that exist. Equivalently, $0\in g+N_C(x)$, where
$N_C(x)=\{v:v^T(z-x)\le0,\ \forall z\in C\}$ is the normal cone.

Simplex notation commonly takes $\pi=-y$, so that the reduced costs are $r=c-A^T\pi=g$. The simplex chapters use $\pi$, whereas the tangent-cone and first-order chapters use $y$; this relation connects the two conventions throughout.

## 3. Chapters and dependencies

| Chapter | Topics |
|---|---|
| [Equivalent transformations and presolve](presolve_en.md) | Standard form, activity intervals, implied bounds, substitution, scaling, and recovery |
| [Two-phase primal simplex](primal-simplex_en.md) | Artificial variables, primal feasible iterations, ratio tests, and pricing |
| [Dual revised simplex](dual-simplex_en.md) | Dual feasibility, leaving-variable pricing, dual steps, bound flipping, and structural dualization |
| [Degeneracy and anti-cycling](anti-degeneracy_en.md) | Degenerate geometry, Bland's rule, symbolic lexicographic ordering, and finite perturbations |
| [Nonnegative least squares and feasibility](nnls_en.md) | Active-set least squares, feasible deficient bases, and infeasibility separation certificates |
| [Pan/BDA dynamic deficient bases](pan-bda_en.md) | Minimum-norm multipliers, basis extension, simultaneous tied deletions, and unbounded directions |
| [Tangent-cone directions and feasibility recovery](tangent-cone_en.md) | Native box constraints, mixed NNLS, support compression, and multiplier refinement |
| [Primal–dual hybrid gradient](pdhg_en.md) | Saddle-point iterations, diagonal preconditioning, averaging, restarts, and hybrid optimization |
| [Basis linear algebra and numerical error](basis-linear-algebra_en.md) | Basis updates, block elimination, QR, rank determination, and iterative refinement |

Presolve supplies equivalent problems to these methods. Two-phase simplex and nonnegative least squares provide distinct feasibility constructions; dual simplex repairs primal feasibility from a dual feasible state. BDA and tangent-cone methods descend from a feasible point, whereas PDHG generates approximate primal–dual pairs through saddle-point iterations. All methods ultimately share the original model's feasibility and optimality conditions.
