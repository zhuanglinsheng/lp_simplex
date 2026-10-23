# Tangent-Cone Directions and Feasibility Recovery with Native Box Constraints

[简体中文](./tangent-cone.md) | English | [Project home](../../README_en.md) | [Algorithm overview](./README_en.md)

## Abstract

Mixed free–nonnegative least squares can construct multi-column descent directions directly on native box constraints. This formulation expresses feasible motion on the current face and interprets auxiliary least-squares residuals through the original problem's KKT conditions.

## 1. Parameterizing the tangent cone

At a feasible point $x$, the tangent cone of the box and equality constraints is

$$
T(x)=\{d:Ad=0,\ d_j\ge0\text{ at a lower bound},\ d_j\le0\text{ at an upper bound},\ d_j=0\text{ at a fixed coordinate}\}.
$$

Directions of interior and free coordinates are unrestricted. Introduce oriented coordinates $d=Sv$: use column $a_j$ at a lower bound and $-a_j$ at an upper bound, with the corresponding $v_j\ge0$. Interior directions are free, and fixed directions are removed. Define

$$
M=AS,\qquad h=S^Tc.
$$

If $Mv=0$ and $h^Tv<0$, then $d=Sv$ is a feasible descent direction. By convexity, absence of such a direction implies optimality of the current feasible point.

## 2. Augmented mixed least squares

For $\gamma>0$, define

$$
C=\begin{bmatrix}M\\h^T/\gamma\end{bmatrix},\qquad
 t=\binom{0}{-1},\qquad
\min_{v_F\in\mathbb R,\ v_N\ge0}\frac12\|Cv-t\|_2^2.
$$

A zero residual gives $Mv=0,h^Tv=-\gamma$, hence a normalized descent direction. The objective can also be written as

$$
\frac12\|Mv\|_2^2+rac12(h^Tv/\gamma+1)^2.
$$

Explicit objective-row augmentation and an implicit scalar objective term are mathematically equivalent; the latter leaves the auxiliary problem unchanged.

## 3. Dual interpretation of a nonzero optimal residual

Suppose the auxiliary KKT conditions hold, and let $e=t-Cv=(p,\eta)$. Mixed least squares yields

$$
C_F^Te=0,\qquad C_N^Te\le0,\qquad v^TC^Te=0.
$$

Thus $t^Te=\|e\|_2^2$, while $t^Te=-\eta$. If $e\ne0$, then $\eta<0$; define

$$
y=\frac{\gamma}{\eta}p.
$$

Using $M^Tp+(h/\gamma)\eta=C^Te$ gives

$$
h_F+M_F^Ty=0,\qquad h_N+M_N^Ty\ge0.
$$

These are exactly the reduced-cost conditions in native direction coordinates, so $x,y$ satisfy the original LP's KKT conditions. If $e=0$, the preceding section instead yields a descent direction. This alternative requires true auxiliary optimality and primal feasibility of the current point. An approximate residual cannot serve as a certificate without verification.

## 4. Two-sided step limits and unbounded rays

For a verified descent direction $d$, take

$$
\alpha_{\max}=\min\left\{
\frac{u_j-x_j}{d_j}:d_j>0,u_j<\infty;
\frac{\ell_j-x_j}{d_j}:d_j<0,\ell_j>-\infty
\right\}.
$$

If there is a finite blocker, every $0\le\alpha\le\alpha_{\max}$ preserves box feasibility. Without a finite blocker, $Ad=0$ and $c^Td<0$ make $d$ an unbounded ray. Numerical calculations must independently verify $A(x+\alpha d)=b$, all variable bounds, and actual objective decrease. In particular, a small auxiliary residual obtained through cancellation of original constraint errors is insufficient.

## 5. Feasibility recovery and support compression

Let the residual of an approximate point $\widetilde x$ be $r=b-A\widetilde x$. After fixing reliable active bounds, first solve on the remaining free face $F$:

$$
\min_{\Delta x_F}\|A_F\Delta x_F-r\|_2^2.
$$

If this correction violates bounds, allow active-bound coordinates to move inward and repair the point subject to $\ell-\widetilde x\le\Delta x\le u-\widetilde x$. The associated bound directions admit a mixed-NNLS representation, with finite steps preserving bound feasibility. Recovery is valid only if the new point satisfies the original constraints, original variable bounds, and the required objective ceiling.

If the free support $F$ is dependent, choose a nonzero direction satisfying $A_Fv=0$ and move along $v$ or $-v$ until at least one free coordinate reaches a bound. Both signs permit a small step on the current free face, so a sign with $c^Tv\le0$ may be selected. With a finite blocker, this compresses the free support without increasing the objective. A descent direction without a finite blocker must instead be analyzed as an unbounded ray. Compression must preserve the equalities; small coordinates cannot simply be deleted.

## 6. Multiplier refinement

In oriented coordinates, let $g(y)=h+M^Ty$. Consider minimizing the convex piecewise-quadratic loss

$$
\Psi(y)=\frac12\sum_{j\in F}\frac{g_j(y)^2}{\omega_j^2}
 +\frac12\sum_{j\in N}\frac{\min\{g_j(y),0\}^2}{\omega_j^2},
\qquad \omega_j>0.
$$

For a fixed violated set $V$, its Hessian is

$$
H=\sum_{j\in F\cup V}\frac{M_jM_j^T}{\omega_j^2}\succeq0.
$$

Gradient descent or approximate Newton corrections with a line search can seek a multiplier. Positive semidefiniteness does not imply invertibility. A decrease in this loss establishes only a reduction in auxiliary sign errors; original-model KKT and duality-gap checks remain necessary.
