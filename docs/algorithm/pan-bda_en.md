# Pan/BDA with Dynamic Deficient Bases

[简体中文](./pan-bda.md) | English | [Project home](../../README_en.md) | [Algorithm overview](./README_en.md)

## Abstract

A basis-deficiency-allowing generalized simplex method represents feasible points using linearly independent active columns without requiring their number to equal the number of constraints. At degenerate ratio ties, several coordinates reaching zero can leave simultaneously, expressing degeneracy as a reduction in active dimension.

## 1. Feasible states and multipliers

Consider

$$
\min_z c^Tz+c_0,\qquad Az=b,\quad z\ge0.
$$

Let $J$ be the active set and let $B=A_J\in\mathbb R^{m\times k}$ have full column rank, with

$$
Bz_J=b,\qquad z_J\ge0,\qquad z_{J^c}=0,\qquad k\le m.
$$

Such a state can be constructed by [nonnegative least squares](nnls_en.md). Compute the minimum-norm multiplier

$$
y^*=\arg\min_y\left\{\tfrac12\|y\|_2^2:B^Ty=-c_J\right\}.
$$

For $k>0$, its algebraic expression is

$$
y^*=-B(B^TB)^{-1}c_J.
$$

For an empty active set, take $y^*=0$. This expression serves the derivation and does not require evaluation through normal equations. Define $g_j=c_j+a_j^Ty^*$; active columns have $g_j=0$.

## 2. Pricing and optimality

If every inactive column has $g_j\ge0$, then $y^*$ is dual feasible and

$$
c^Tz=-y^{*T}Bz_J=-b^Ty^*.
$$

The primal and dual objectives coincide, proving optimality of the current point. If some $g_j<0$, choose a column maximizing $-g_j$, or use normalized pricing $-g_j/\|a_j\|_2$.

When $k<m$, the multiplier equation has multiple solutions. A violating column for the minimum-norm multiplier does not necessarily imply that the current feasible point is nonoptimal. A basis extension may eliminate the violation by changing the multiplier representation.

## 3. Extension by an independent column

For entering column $a_q$, form the orthogonal projection

$$
\delta=\arg\min_v\|Bv-a_q\|_2,\qquad h=a_q-B\delta.
$$

If $h\ne0$, the augmented column set $[B,a_q]$ retains full column rank. Extend the basis with $z_q=0$. The feasible point and objective remain unchanged, and the active column count becomes $k+1$. The next multiplier satisfies the additional equation $a_q^Ty=-c_q$, allowing the previous pricing violation to be corrected.

## 4. Dependent-column exchanges and simultaneous deletions

If $h=0$, then $a_q=B\delta$. Movement along

$$
z_q(\theta)=\theta,\qquad z_J(\theta)=z_J-\theta\delta
$$

preserves the equalities and has objective slope

$$
c_q-c_J^T\delta=c_q+a_q^Ty^*=g_q<0.
$$

If some $\delta_i>0$, take

$$
\theta^*=\min_{\delta_i>0}\frac{z_{J_i}}{\delta_i},\qquad
Q=\left\{i:\delta_i>0,\ \frac{z_{J_i}}{\delta_i}=\theta^*\right\}.
$$

After updating, delete all zero components indexed by $Q$ and add $q$. The new active column count is

$$
k'=k-|Q|+1.
$$

The new column set remains independent. Indeed, if $\lambda a_q+\sum_{i\notin Q}\mu_iB_i=0$, substitute $a_q=B\delta$. At any $i\in Q$, independence of the original columns gives $\lambda\delta_i=0$. Since $\delta_i>0$, it follows that $\lambda=0$, and then all $\mu_i=0$.

Simultaneous deletion thus preserves both feasibility and column independence. When $|Q|>1$, a smaller basis is obtained without filling it with zero-valued columns to make it square. In finite precision, approximate ratio ties require verification of the equality residual after deletion. Near-zero coordinates cannot be treated unconditionally as exact zeros.

## 5. Unboundedness and termination conditions

If all $\delta_i\le0$ in the dependent case, the path remains feasible for every $\theta\ge0$ and decreases the objective with negative slope $g_q$, proving unboundedness. A zero column with $c_q<0$ is a special case.

Optimal termination requires primal feasibility, the active multiplier equation, and reduced-cost signs for all columns. Primal feasibility may be checked by resolving the full-column-rank least-squares problem for active values. For nearly dependent columns, the reliability of projection residuals and rank decisions directly affects basis extension and ray identification; both require checks against the original matrix.

Allowing basis deficiency improves the representation of degeneracy, but these local invariants alone do not prove finite termination for arbitrary floating-point problems. Combinations with [multi-column tangent-cone directions](tangent-cone_en.md) and [PDHG hybrid methods](pdhg_en.md) must also satisfy the direction, feasibility, and optimality conditions of each constituent phase.
