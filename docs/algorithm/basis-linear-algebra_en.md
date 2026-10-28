# Basis Linear Algebra, Orthogonal Factorization, and Numerical Error

[简体中文](./basis-linear-algebra.md) | English | [Project home](../../README_en.md) | [Algorithm overview](./README_en.md)

## Abstract

Simplex and deficient-basis methods rely on linear equations, least squares, minimum-norm multipliers, and rank decisions. This chapter studies only their algebraic relations and error properties.

## 1. Square-basis solves and single-column updates

Two basic operations for a full basis are

$$
Bd=a_q,\qquad B^T\rho=e_p.
$$

The first determines the entering direction; the second determines the leaving pricing row. They satisfy

$$
d_p=e_p^TB^{-1}a_q=\rho^Ta_q.
$$

This identity permits an independent consistency check between the two calculations. An LU factorization can be written as $PBQ=LU$, where $P,Q$ are permutation matrices, reducing a solve to two triangular systems.

When $a_q$ replaces column $p$ of $B$, let $d=B^{-1}a_q$ and let $E$ be the identity with column $p$ replaced by $d$. Then

$$
B'=BE,\qquad (B')^{-1}=E^{-1}B^{-1},\qquad \det E=d_p.
$$

If $d_p\ne0$, the solution $z=E^{-1}v$ is

$$
z_p=v_p/d_p,\qquad z_i=v_i-d_i z_p\quad(i\ne p).
$$

Successive exchanges correspond to a product of such transformations. Small $|d_p|$ or large direction components amplify rounding error. Refactorization after accumulated updates can control representation error, but does not remove ill-conditioning of the basis itself.

## 2. Exact block elimination for logical columns

After suitable row and column permutations, a basis containing logical columns has the form

$$
B=\begin{bmatrix}K&0\\H&-I\end{bmatrix}.
$$

For $Bx=f$, first solve $Kx_R=f_R$, then set $x_D=Hx_R-f_D$. For $B^Ty=g$, first set $y_D=-g_D$, then solve

$$
K^Ty_R=g_R+H^Tg_D.
$$

Only the structural block $K$ therefore needs factorization. This is exact block elimination rather than an approximation of the original constraints.

Solution lifting also uses sparsity. The product $Hx_R$ skips only exactly zero core components. For $H^Tg_D$, an optional structural-column inverse map and the existing CSR view visit basic columns adjacent to nonzero logical components. Comparing the candidate row entries with the core CSC nonzero count selects the cheaper traversal and avoids expensive row scans for dense right-hand sides. Missing CSR or refused optional storage retains CSC evaluation; residual audits and refinement are unchanged.

KLU factors can additionally expose an optional support-driven solve graph through the public extraction API. Separate L, U and off-block stages preserve the block traversal; transpose solves collect contributions in stored-column order. Sparse event propagation falls back to the existing solve when support grows or bounded storage is refused. Packed RHS support also passes through logical-block recovery and mixed eta chains. No small-value dropping is used; full residual checks and refinement remain authoritative. Integer support sorting borrows the existing slot array and rebuilds its mappings.

KLU eta chains use the existing deterministic amortization estimate: accumulated update-application work reaching estimated factor work requests refactorization, alongside update-length and weak-pivot protection. Large internal-LU bases retain their size guard against expensive refactorization; original row count no longer disables KLU's cost rule as well. This changes factor rebuild scheduling without relaxing residual or success conditions.

## 3. Rectangular-basis QR and minimum-norm solutions

Let $B\in\mathbb R^{m\times k}$ have full column rank. Its thin QR factorization is

$$
B=QR,\qquad Q^TQ=I_k,
$$

where $R$ is nonsingular and upper triangular. The least-squares projection and residual are

$$
\delta=R^{-1}Q^Ta,\qquad h=(I-QQ^T)a.
$$

The minimum-norm solution of $B^Ty=v$ is

$$
y=QR^{-T}v.
$$

Every other solution adds a vector in $\mathcal R(B)^\perp$, and therefore cannot have smaller norm. Although the normal equations $B^TB\delta=B^Ta$ are algebraically equivalent, $\kappa_2(B^TB)=\kappa_2(B)^2$ in the Euclidean norm, increasing sensitivity to errors in ill-conditioned problems.

## 4. Orthogonal extensions and downdates

When adding column $a$, first compute $q=Q^Ta$ and $h=a-Qq$. If $\|h\|_2>0$, then

$$
[B,a]=[Q,h/\|h\|_2]
\begin{bmatrix}R&q\\0&\|h\|_2\end{bmatrix}.
$$

If $h=0$, the new column is dependent and cannot be treated as an independent extension. After deleting an active column, orthogonal transformations can restore triangular form to the remaining factor. Householder transformations $H=I-2vv^T/(v^Tv)$ and Givens plane rotations preserve the Euclidean norm and express these extensions and downdates.

Floating-point rank decisions must be interpreted relative to column scales and error bounds. A very small projection residual indicates numerical near-dependence only. Rank-revealing factorizations distinguish reliably independent columns from uncertain directions; inability to determine rank reliably does not prove LP infeasibility.

## 5. Residuals, backward error, and refinement

For an approximate solution $\widehat z$, define $r=f-B\widehat z$. A componentwise backward error is

$$
\eta=\max_i\frac{|r_i|}{(|B||\widehat z|+|f|)_i},
$$

A zero numerator and denominator contribute zero; a nonzero numerator with a zero denominator contributes infinity. Iterative refinement solves $B\Delta z=r$ and sets $\widehat z\leftarrow\widehat z+\Delta z$. Its effectiveness depends on residual evaluation, correction accuracy, and the condition number; it does not guarantee repair of every ill-conditioned system.

Final primal residuals must be computed using the original matrix and variable bounds. Residuals from scaling, substitution, and objective perturbation belong to different models and cannot be compared directly without inverse transformations. Primal feasibility alone also cannot establish optimality when reduced-cost signs are incorrect.

If every variable has finite two-sided bounds, the gap at a primal feasible point for a given multiplier is

$$
c^Tx+c_0-q(y)=\sum_j\left[g_jx_j-\inf_{\ell_j\le t\le u_j}g_jt\right].
$$

For a nonbasic variable at an endpoint, an incorrect reduced-cost sign contributes its bound width times the corresponding sign violation to the gap. Small coordinatewise reduced costs may therefore accumulate into a substantial objective error on wide-bound problems. With infinite bounds and a nonfinite dual function, this width-based estimate is not a strict duality-gap certificate.
