# The Dual Revised Simplex Method

[简体中文](./dual-simplex.md) | English | [Project home](../../README_en.md) | [Algorithm overview](./README_en.md)

## Abstract

Dual simplex preserves compatibility between reduced-cost signs and nonbasic bound states while progressively eliminating bound violations of basic variables. This chapter derives pricing, ratio tests, and bound flipping for bounded variables.

## 1. Basis relations and initialization

For $Ax=b,\ell\le x\le u$, fix nonbasic variables at selected endpoints. The basic values satisfy

$$
Bx_\mathcal B=b-A_\mathcal N x_\mathcal N,\qquad
B^T\pi=c_\mathcal B,\qquad r=c-A^T\pi.
$$

Dual feasibility requires $r_j\ge0$ for a nonbasic variable at its lower bound, $r_j\le0$ at its upper bound, and $r_j=0$ for a free nonbasic variable; fixed variables have no sign restriction. Basic reduced costs are zero. If also $\ell_\mathcal B\le x_\mathcal B\le u_\mathcal B$, optimality follows.

In the row-variable formulation, logical columns form $-I$ and supply a nonsingular initial basis. This alone does not guarantee dual feasibility. Bound choices, basis exchanges, or an auxiliary objective must establish a dual feasible state. Any auxiliary changes to costs must be undone before final optimality is assessed.

## 2. Leaving rows and steepest-edge pricing

Define the bound violation at basic position $i$ as

$$
v_i=\max\{\ell_{\mathcal B_i}-x_{\mathcal B_i},\,x_{\mathcal B_i}-u_{\mathcal B_i},\,0\}.
$$

The dual steepest-edge weight is $\gamma_i=\|B^{-T}e_i\|_2^2$. Choose position $p$ maximizing $v_i^2/\gamma_i$, and set

$$
\rho=B^{-T}e_p,\qquad \alpha_j=\rho^Ta_j.
$$

If $x_{\mathcal B_p}<\ell_{\mathcal B_p}$, take repair sign $\kappa=1$; if it exceeds its upper bound, take $\kappa=-1$. This convention unifies the two cases below.

Cancellation can make recursively updated weights underestimate edge lengths. After choosing a leaving row, reuse its BTRAN result to recompute its squared norm. If the old weight is less than one quarter of that norm, update the feasibility heap and reselect. At a fixed basis, at most all rows are corrected before one final selection; repricing does not advance pivot work or change the point or objective. A nonfinite or nonpositive norm causes numerical failure. This ratio bounds optimistic pricing error without changing feasibility or optimality tolerances.

Dense CSC pricing and sparse CSR scatter share the same pricing-row storage. Dense writes must register their coordinates too, so the next indexed clear removes every old value; untouched coordinates must be exactly zero. Periodic full CSC audits do not replace this storage invariant.

## 3. Dual step lengths

Consider the multiplier path and reduced costs

$$
\pi(\theta)=\pi-\kappa\theta\rho,\qquad
r_j(\theta)=r_j+\kappa\theta\alpha_j.
$$

For a nonfixed nonbasic variable at an endpoint, let $s_j=1$ denote the lower bound and $s_j=-1$ the upper bound. Dual feasibility is $s_jr_j\ge0$. Only variables satisfying

$$
\beta_j=-\kappa s_j\alpha_j>0
$$

lose dual feasibility as $\theta$ increases. Their breakpoints are

$$
\theta_j=\frac{s_jr_j}{\beta_j}.
$$

Without bound flipping, take the smallest breakpoint and let its variable enter. A free nonbasic variable must retain zero reduced cost, so $\alpha_j\ne0$ gives a blocker at zero dual step. Such a variable cannot be assigned a finite endpoint that does not exist.

With feasibility tolerance $t$, Harris' first-pass bound is $(s_j r_j+t)/\beta_j$. If $s_jr_j$ is already negative within tolerance, only the primary breakpoint is clipped to zero; the first-pass bound must retain that negative value. Replacing it by zero would replenish spent feasibility allowance at every pivot. The second pass compares primary ratios inside this bound without a ratio-space epsilon expansion, whose reduced-cost effect depends on the column coefficient.

## 4. Primal steps and basis exchanges

Let $q$ enter, $d=B^{-1}a_q$, and $h_p$ be the endpoint reached by the leaving variable. Set

$$
\Delta x_q=\frac{x_{\mathcal B_p}-h_p}{d_p},\qquad
\Delta x_\mathcal B=-d\,\Delta x_q.
$$

Since $d_p=\alpha_q$, the sign conditions ensure that $\Delta x_q$ moves into the allowed direction from its endpoint. Replace column $p$ of the basis by $a_q$. The dual step makes $r_q$ zero, and the leaving variable satisfies the required reduced-cost sign at its endpoint.

The new basic variable may temporarily violate its own bounds because dual feasibility is the invariant of this phase. Requiring every exchange to make all basic variables feasible would exclude valid dual-simplex steps.

## 5. Bound flipping

A candidate with two finite bounds may move from one endpoint to the other when its reduced cost crosses zero. The change in basic values is

$$
\Delta x_\mathcal B=-B^{-1}a_j\Delta x_j.
$$

After reversing the nonbasic bound sign, the dual path may continue. Multiple flips have total effect

$$
\Delta x_\mathcal B=-B^{-1}\sum_{j\in F}a_j\Delta x_j.
$$

Breakpoints must be processed in monotone dual-step order. Flips and the eventual exchange must form a consistent primal–dual update; moving the multiplier alone cannot disregard the basic equation $B^T\pi=c_\mathcal B$.

If the leaving variable still violates its bound and no nonbasic variable can improve that violation in an allowed direction, the identity $\rho^TAx=\rho^Tb$ and the variable box may yield a contradiction, providing an infeasibility separation proof. An empty floating-point candidate set is insufficient by itself: the full pricing row and separating inequality must be verified.

## 6. Structural dualization of paired singleton columns

Consider

$$
\min\ c_z^Tz+(c^+)^Tp+(c^-)^Tq,
\quad A_zz+D^+p+D^-q=b,\quad z,p,q\ge0,
$$

where $D^+,D^-$ are diagonal, $d_i^+>0$, $d_i^-<0$, and $c_i^+,c_i^->0$. The dual is

$$
\max_\pi b^T\pi,\qquad
A_z^T\pi\le c_z,\qquad
\frac{c_i^-}{d_i^-}\le\pi_i\le\frac{c_i^+}{d_i^+}.
$$

The paired singleton columns thus become finite two-sided bounds on dual variables. After recovering $z$ from optimal multipliers of the corresponding dual constraints, set $r=b-A_zz$ and take

$$
p_i=\frac{\max(r_i,0)}{d_i^+},\qquad
q_i=\frac{\min(r_i,0)}{d_i^-}.
$$

This gives a primal feasible solution. With exact optimal multipliers and complementary slackness, it is also primal optimal. The reduction depends on the stated structure; a general equality model need not have finite dual box bounds.

The auxiliary classification is not an input classification. The signed-pair structure establishes input feasibility; certified auxiliary infeasibility inside a conservative finite box establishes input unboundedness by LP duality. Auxiliary unboundedness contradicts the box and is a numerical failure. Auxiliary optimality still requires recovered input point/multiplier KKT audit; see the [Singleton handoff contract](../auto-policy_en.md#candidates-checkpoints-and-termination).


## 7. Implementation boundaries and call path

Dual separates its algorithm entry, numerical-backend recovery and individual numerical attempt:

```mermaid
flowchart LR
    C["Caller: AUTO / BDA service / explicit Dual"] --> E["simplex_dual.c: input adaptation and recovery"]
    E --> B["simplex_dual_backend.c: attempts and shared grant"]
    B --> K["simplex_dual_kernel.c: iteration and certification"]
```

| Layer | Inputs and outputs | Responsibility and exits |
|---|---|---|
| [Entry](../../src/dual/simplex_dual.c) | Immutable invocation input, options, optional model-bound hint → point, multipliers, status and diagnostics in that input space | Check hint identity before scaling; create scaling/progress bridge, falling back to the original representation for irreversible scaling; recover coordinates and recompute input diagnostics. Reject bad hints or cold-start under existing rules; exit on resource/bridge failure. |
| [Backend](../../src/dual/simplex_dual_backend.c) | Fixed-representation input, complete Dual grant, retry key, optional hint → selected candidate and combined work | After default-LU numerical failure/singularity, try internal LU only when compiled, eligible under the exact-request ledger and within remaining work. Save, audit and restore candidates and producers; do not rerun mathematical conclusions or resource/budget stops. |
| [Kernel](../../src/dual/simplex_dual_kernel.c) | Borrowed kernel-space input, explicit factorization backend, attempt work grant → point/multipliers/optional basis, status and actual work | Own the attempt's DualState and factors; initialize, iterate, perform local numerical recovery and mathematical certification, then destroy state. Do not create scaling models or choose another backend. |

See the private [backend](../../src/dual/simplex_dual_backend.h) and [kernel](../../src/dual/simplex_dual_kernel.h) contracts. An invocation input may already be an outer reduced model; it is not always the public original LP. A basis hint transfers model-bound indices/status only and still requires refactorization and classification. No live factor transfers between backends. Diagnostic parents and the numerical failure ledger retain distinct roles and do not replace model identity.

Backend point selection only assembles output. Actual attempts retain local work; the selection summary carries their inclusive total. Both attempts and selection belong to the invocation owner, with distinct triggering and point-producing sources. Inclusive work is neither reported as one backend's local work nor nested under the last attempt. See the [cross-layer recovery chain](../auto-policy_en.md#cross-layer-recovery-chain).

These files implement one Dual algorithm. LU replacement is its numerical-backend recovery; AUTO still owns cross-algorithm selection, grants and candidate commits among Dual, Tableau and BDA. Checkpoints, request comparison, bounds, pricing, degeneracy and state-lifecycle modules retain their roles. The kernel still contains coupled iteration/cleanup logic; this separation does not close all large-file, transformed-model or resource-exit acceptance.
