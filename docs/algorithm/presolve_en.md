# Equivalent Model Transformations and Presolve

[简体中文](./presolve.md) | English | [Project home](../../README_en.md) | [Algorithm overview](./README_en.md)

## Abstract

Presolve reduces the dimension of a linear program through algebraic transformations that preserve the projection of its feasible set and its objective value. Correctness requires every reduced feasible point to admit recovery to an original feasible point, with objective values differing only by a known constant. The discussion begins with a row-bound model.

## 1. Equality models and standard form

For $\ell^r\le Ax\le u^r$, introduce $s=Ax$ to obtain

$$
[A\;-I]\binom{x}{s}=0,\qquad \ell^r\le s\le u^r.
$$

This formulation preserves native variable bounds. If nonnegative standard form is required, use

$$
\begin{array}{c|c}
\text{Variable type}&\text{Transformation}\\\hline
x_j\ge\ell_j&x_j=\ell_j+z_j,\ z_j\ge0\\
x_j\le u_j&x_j=u_j-z_j,\ z_j\ge0\\
\ell_j\le x_j\le u_j&x_j=\ell_j+z_j,\ z_j+w_j=u_j-\ell_j,\ z_j,w_j\ge0\\
x_j\text{ free}&x_j=z_j^+-z_j^-,\ z_j^+,z_j^-\ge0\\
x_j=\ell_j=u_j&x_j=\ell_j.
\end{array}
$$

Adding positive slack variables to less-than-or-equal constraints and negative slack variables to greater-than-or-equal constraints yields $\widetilde A z=\widetilde b,z\ge0$. If $x=h+Tz$, the transformed objective is

$$
\widetilde c=T^Tc,\qquad \widetilde c_0=c_0+c^Th.
$$

A maximization problem becomes a minimization problem by negating its objective; the original sign is restored during recovery.

## 2. Fixed variables, empty rows, and empty columns

A fixed variable $x_j=\alpha$ can be eliminated by setting $b\leftarrow b-a_j\alpha$ and $c_0\leftarrow c_0+c_j\alpha$. A zero row is feasible precisely when $\ell_i^r\le0\le u_i^r$. A zero column requires only independent minimization of $c_jx_j$: choose the finite lower bound if $c_j>0$, the finite upper bound if $c_j<0$, and any feasible value if $c_j=0$. If the improving endpoint is infinite, an unbounded ray exists provided the remaining constraints are feasible. Until that feasibility is known, the conclusion remains conditional.

## 3. Activity intervals and implied bounds

Define lower and upper activity bounds for row $i$ by

$$
L_i=\sum_j\min\{a_{ij}\ell_j,a_{ij}u_j\},\qquad
U_i=\sum_j\max\{a_{ij}\ell_j,a_{ij}u_j\}.
$$

Terms with zero coefficients are defined as zero. If $L_i>u_i^r$ or $U_i<\ell_i^r$, the problem is infeasible. If $[L_i,U_i]\subseteq[\ell_i^r,u_i^r]$, the row is redundant. If a finite activity extremum equals the row bound that forces attainment of that extremum, every nonzero term must attain its corresponding endpoint. The associated variables can therefore be fixed; such a row is called a forcing row.

For $a_{ij}\ne0$, let $L_{i,-j},U_{i,-j}$ denote activity bounds excluding term $j$. The inequalities

$$
\ell_i^r-U_{i,-j}\le a_{ij}x_j\le u_i^r-L_{i,-j}
$$

imply tighter bounds on $x_j$; division reverses the inequalities when $a_{ij}<0$. Singleton rows are a special case. Infinite endpoints require interval arithmetic: undefined subtraction of infinities cannot justify a new bound. Safe tightening and elimination are repeated until no further effective reduction is obtained.

If two rows satisfy $a_{kj}=\lambda a_{ij}$ for all $j$, with $\lambda\ne0$, divide the bounds of row $k$ by $\lambda$ and intersect them with those of row $i$, exchanging endpoints when $\lambda<0$. An empty intersection proves infeasibility; otherwise, retain the stronger common interval.

### 3.1 Dominance and preservation of the optimum

For nonnegative variables, negate less-than-or-equal rows to obtain covering inequalities. If $0\le a\le d$ componentwise and $b\ge e$, then $a^Tx\ge b$ implies $d^Tx\ge e$, so the latter row can be deleted. Coefficients need not be binary and right-hand sides need not coincide.

For minimization, suppose column $k$ costs no more than column $j$, and decreasing $x_j$ while increasing $x_k$ equally never weakens an inequality or changes an equality activity. If $k$ has no upper bound, $j$ can be fixed at its finite lower bound. This rule permits signed coefficients; coefficients on equality rows must match.

A finite upper bound on $k$ invalidates that argument in general. The implementation permits clipping the transfer only under sufficient conditions: $j$ has lower bound zero, both costs are nonnegative, every affected row is a monotone covering row on nonnegative variables, and $k$ at its upper bound can satisfy every affected row alone. Clipping an overflowing transfer then preserves coverage without increasing cost. Insufficient capacity, negative costs, or incompatible signs cannot use this bounded rule.

Column fixing preserves the optimal value and at least one optimum, not every original feasible point or every optimum. Existing postsolve restores fixed values without aggregating remaining column capacities. The implementation uses exact coefficient comparisons and checks capacity-product rounding. Its bit index is capped at 32 MiB; larger models use sparse posting lists and support fingerprints. Fingerprints only reject impossible candidates; hash collisions never authorize deletion. Compact postings retain movable columns and prioritize cost and support size, while acceptance still checks original coefficients, costs, bounds, and row directions. All fixed-point passes share a finite search budget. Unsorted or duplicate-index sparse views skip these rules, and budget exhaustion only skips reductions that have not been proved.

### 3.2 Witnessed equality dependency elimination

Currently disabled by default and enabled only through explicit experimental configuration; complete-solve performance has not passed acceptance. See [configuration boundaries](../runtime-configuration_en.md).

For an equality subsystem, suppose a sparse integer matrix $W$ satisfies $WA=0$ and $Wb=0$, and the square submatrix $W_D$ on selected columns is invertible. The original equations indexed by $D$ may be removed together: $W_D(A_Dx-b_D)=-W_K(A_Kx-b_K)$, so retained equations imply deleted equations. Variables, objective and retained matrix coefficients remain unchanged.

The implementation generates local and cross-group conservation witnesses from adjacency of exact signed-unit marginal equations. Overflow-checked integer sums verify complete coefficients and RHS before a nonzero minor modulo a prime proves independence. Approximate coefficients are excluded; modular zero alone cannot establish redundancy. Bounded search or storage refusal only loses optional reductions. After publishing the deletion set, ordinary queues and subsequent substitutions resume, and final acceptance audits the original model. Budgets are specified in the [numerical policy](../numerical-policy_en.md).

## 4. Equality substitution and projection

An equality with nonzero coefficient $a_{ij}$ gives

$$
x_j=\frac{b_i-\sum_{k\ne j}a_{ik}x_k}{a_{ij}}.
$$

An equivalent problem results from substituting this expression into the other constraints and the objective while retaining the projected constraints induced by $\ell_j\le x_j\le u_j$. If the variable occurs only in this row, no other row changes. If the row contains only two variables, the substitution is an affine relation in one variable.

For example, $a x_p+d x_q=\beta$ gives

$$
x_p=\beta/a-(d/a)x_q,\qquad
\widehat c_q=c_q-c_p d/a,\qquad
\widehat c_0=c_0+c_p\beta/a.
$$

A projected bound may be omitted only if the substitution expression satisfies the eliminated variable's bound throughout the admissible range of the remaining variables. Approximate satisfaction at one point or a small local residual does not establish redundancy.

Different retained constraints may imply the two bounds. For example, $x+a-b=1$ and $x-c-d=0$, with $a,c,d\ge0$ and $b\le8$, together imply $0\le x\le9$. Eliminating $x$ leaves $a-b+c+d=1$, which preserves both implications with the remaining bounds. Recheck the current matrix after every elimination instead of retaining a deleted proof source.

The additional cross-row proof currently accepts signed-unit equalities only. TwoSum errors and outward rounding enclose the interval endpoints; user feasibility tolerance never authorizes dropping a bound. A pivot row of length $r$ and column of degree $d$ delete $r+d-1$ entries and insert at most $(r-1)(d-1)$; accept only nonincreasing fill. The additional degree-three/four rule also requires signed-unit incident rows, conservatively limiting numerical changes in mixed-coefficient representations, independently of dataset names. Reverse journal replay uses compensated addition and FMA product remainders to retain small coordinates through cancellation. Final acceptance still requires the original-model audit.

## 5. Scaling and inverse transformations

For an equality-constrained box model, choose positive diagonal matrices $R,D$ and set $x=Dz$. Then

$$
\widehat A=RAD,\quad \widehat b=Rb,\quad
\widehat c=Dc,\quad D^{-1}\ell\le z\le D^{-1}u.
$$

If $\widehat y$ is a multiplier for the scaled problem, the original multiplier is $y=R\widehat y$. Scaling can improve numerical balance but does not change exact optimality or repair genuine rank deficiency.

After several substitutions, recovery must proceed in reverse dependency order. Since $x=h+Tz$, a reduced-variable error $\Delta z$ induces original-variable error $T\Delta z$. A small reduced residual therefore need not imply a comparable original-model residual. Constraint violations and the objective must ultimately be recomputed in the original model. The primal recovery map alone does not recover all original-model dual multipliers.
