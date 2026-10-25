# The Two-Phase Primal Simplex Method

[简体中文](./primal-simplex.md) | English | [Project home](../../README_en.md) | [Algorithm overview](./README_en.md)

## Abstract

Primal simplex moves between basic feasible solutions while preserving primal feasibility and improving the objective. The two-phase construction separates feasibility from optimization of the original objective. Full tableaux and revised simplex are two representations of the same basis algebra.

## 1. Basis dictionaries and reduced costs

Consider $\min c^Tx+c_0$ subject to $Ax=b,x\ge0$. Assume that $A$ has full row rank and choose a nonsingular basis $B=A_\mathcal B$. Then

$$
x_\mathcal B=B^{-1}b-B^{-1}A_\mathcal N x_\mathcal N.
$$

Let $B^T\pi=c_\mathcal B$ and $r_j=c_j-a_j^T\pi$. The objective satisfies

$$
c^Tx+c_0=c_\mathcal B^TB^{-1}b+c_0+\sum_{j\in\mathcal N}r_jx_j.
$$

If $B^{-1}b\ge0$ and all nonbasic reduced costs are nonnegative, $\pi$ is dual feasible and the current basic solution is optimal.

## 2. Phase I: artificial variables

Negate equations with negative right-hand sides so that $b\ge0$, and construct

$$
\min_{x,a}\ \mathbf1^Ta,\qquad Ax+Ia=b,\qquad x,a\ge0.
$$

The artificial columns provide the feasible initial solution $x=0,a=b$. Since the auxiliary objective is nonnegative, its optimum is zero if and only if the original problem is feasible. A strictly positive auxiliary optimum proves infeasibility.

At a zero auxiliary optimum, artificial variables may remain basic at zero. If the corresponding dictionary row admits a replacement by an original-variable column, a zero-length pivot removes the artificial column. Otherwise, that equality is redundant in the remaining system and may be deleted. Once all artificial variables have been removed, redefine multipliers and reduced costs using the original objective and begin Phase II.

## 3. Phase II: descent and the ratio test

Choose a nonbasic variable with $r_q<0$ and let $d=B^{-1}a_q$. Along

$$
x_q(\theta)=\theta,\qquad
x_\mathcal B(\theta)=x_\mathcal B-\theta d
$$

the objective changes by $\theta r_q$. If some $d_i>0$, the largest feasible step is

$$
\theta^*=\min_{d_i>0}\frac{x_{\mathcal B_i}}{d_i}.
$$

A variable attaining this minimum leaves the basis. If no $d_i>0$ exists, the direction remains feasible for every $\theta\ge0$ and strictly decreases the objective, yielding an unbounded ray. A step with $\theta^*=0$ is degenerate; anti-cycling rules are discussed in [Degeneracy and anti-cycling](anti-degeneracy_en.md).

## 4. Extension to native variable bounds

For $Ax=b,\ell\le x\le u$, assign movement sign $s_q=1$ to a nonbasic variable at its lower bound and $s_q=-1$ at its upper bound. For a free variable, choose a direction with $s_qr_q<0$. Define

$$
v_q=s_q,\qquad v_\mathcal B=-s_qB^{-1}a_q.
$$

All finite endpoints determine the step limit:

$$
\theta^*=\min\left\{
\frac{u_j-x_j}{v_j}:v_j>0,u_j<\infty;
\frac{\ell_j-x_j}{v_j}:v_j<0,\ell_j>-\infty
\right\}.
$$

The minimum includes the entering variable's opposite endpoint. If that endpoint is reached first, only a bound flip is required. If a basic variable reaches a bound first, exchange basis columns. With no finite blocker and $c^Tv<0$, the vector $v$ is an unbounded direction of the native box model.

## 5. Pricing and stable ratio tests

Dantzig pricing chooses the most negative directional reduced cost $s_jr_j$. Steepest-edge pricing divides its magnitude by the Euclidean norm of the corresponding full-space edge direction. Devex approximates edge length in reference coordinates. For a reference set $\mathcal R$, the quantity

$$
\omega_q=\mathbf1_{q\in\mathcal R}
 +\sum_i\mathbf1_{\mathcal B_i\in\mathcal R}(B^{-1}a_q)_i^2
$$

represents a reference squared edge norm, whose estimate can be updated through basis exchanges. The reference set is rebuilt when the estimate departs substantially from the reference norm. These rules choose the trajectory; they do not replace optimality conditions.

The Harris two-pass ratio test first determines a step window within the permitted feasibility error, then favors a larger pivot within that window to reduce error amplification. This is a numerical rule. A strict anti-cycling conclusion still requires an exact minimum ratio and suitable tie-breaking; it does not follow from a tolerance window alone.

If the working objective contains reversible perturbations, restore the original costs once a primal feasible point is obtained and refine optimality from that point. Optimality for the working objective does not imply optimality for the original objective.

Tableau Phase I, artificial-basis cleanup and Phase II share one cumulative budget. Each actual cleanup pivot consumes a work unit; exhaustion stops before the next pivot and is distinct from numerical failure. Existing pricing/terminal-check work accounting is retained.

The private report records input preparation, Phase I, artificial-basis cleanup and Phase II separately; phase work sums to invocation work. A complete point is exported only after successful Phase II; writing it does not imply independent acceptance. This report exports no multipliers, basis, factorization or lower bound. Actual `tableau` invocation events and AUTO final reports retain algorithm facts and controller conclusions separately; see the [AUTO handoff contract](../auto-policy_en.md#candidates-checkpoints-and-termination).

Bound conversion is a model representation: its standard form has a distinct runtime identity, `tableau-model` links the input, and `tableau-standard` records the numerical kernel on that representation. On success, the outer call recovers input coordinates and retains the actual kernel source; failure exports no recovered point. Standard form is not another algorithm.

The Tableau implementation combines Bland/Dantzig entering selection with a shared two-pass stable leaving window and explicit unit-column maintenance. Bland entering indices do not establish the strict minimum-ratio anti-cycling theorem. Tableau-coordinate tolerance still requires input-model audit; see [numerical policy](../numerical-policy_en.md#tableau).
