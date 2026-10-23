# Nonnegative Least Squares and Feasible Deficient Bases

[简体中文](./nnls.md) | English | [Project home](../../README_en.md) | [Algorithm overview](./README_en.md)

## Abstract

Nonnegative least squares reduces feasibility of a standard-form linear program to a convex quadratic problem. An active-set method alternates between nonnegative points and unconstrained least-squares solutions, and can produce an independent active column set that need not be square.

## 1. Auxiliary problem and separation conditions

Consider

$$
\min_{z\ge0}\ \phi(z)=\frac12\|Az-b\|_2^2.
$$

Let $e=b-Az$ and $w=A^Te$. The KKT conditions are

$$
z\ge0,\qquad w\le0,\qquad z_jw_j=0.
$$

If the optimal residual is zero, $z$ is feasible for $Az=b,z\ge0$. If the optimal residual is nonzero, then

$$
A^Te\le0,\qquad b^Te=(Az+e)^Te=\|e\|_2^2>0.
$$

Consequently, $h=-e$ satisfies $A^Th\ge0,b^Th<0$. This is a Farkas separation certificate proving infeasibility of the nonnegative equality system. The derivation requires the auxiliary KKT conditions to hold. Exhausting an iteration budget, failing a rank test, or observing temporary residual stagnation does not establish that premise.

## 2. Active-set iterations

Let $P$ be the passive set, with all variables outside it fixed at zero. One may start with $P=\varnothing,z=0$. If $w_j>0$, that coordinate has a negative objective directional derivative and may enter $P$. Pricing by $w_j/\|a_j\|_2$ reduces the influence of column scaling; zero columns provide no useful descent.

On $P$, solve

$$
\widehat z_P\in\arg\min_v\|A_Pv-b\|_2^2,
\qquad \widehat z_{P^c}=0.
$$

Accept this solution if every passive component is strictly positive. Otherwise, move along $z(\alpha)=z+\alpha(\widehat z-z)$ with

$$
\alpha^*=\min_{i\in P:\widehat z_i\le0,\ z_i>\widehat z_i}
\frac{z_i}{z_i-\widehat z_i}.
$$

Remove components that reach zero and resolve least squares on the smaller passive set. A newly added zero component that cannot move positively must also be removed. The segment step preserves nonnegativity; least-squares optimality and convexity ensure that the objective does not increase.

## 3. Independent support and feasible deficient bases

If active columns remain linearly independent, $|P|\le m$ and the least-squares solution is unique. At zero residual, let $B=A_P$. Then

$$
Bz_P=b,\qquad z_P\ge0,\qquad \operatorname{rank}(B)=|P|\le m.
$$

This is precisely a feasible initial state for a deficient-basis method. A dependent support can be reduced through null-space directions that eliminate redundant positive coordinates, or a rank-revealing factorization can select an independent representation. Small components cannot simply be deleted without checking the equalities.

Rejection of a column by a numerical rank test does not establish exact dependence in the original matrix. An infeasibility conclusion therefore requires auxiliary KKT and separation relations verified against that matrix.

## 4. Mixed free and nonnegative variables

More generally, let $d_F$ be free and $d_N\ge0$, and solve

$$
\min_{d_F,\,d_N\ge0}\frac12\|C_Fd_F+C_Nd_N-t\|_2^2.
$$

Free coordinates always participate in least squares; only constrained coordinates enter the passive set or leave upon reaching zero. With $e=t-Cd$, the KKT conditions are

$$
C_F^Te=0,\qquad C_N^Te\le0,\qquad d_N^TC_N^Te=0.
$$

This formulation supports tangent-cone descent and feasibility recovery. Its objective measures auxiliary fitting error only. Any connection to optimality of the original linear program requires the relevant geometric derivation.
