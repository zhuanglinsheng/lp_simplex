# Degeneracy, Perturbations, and Anti-Cycling

[简体中文](./anti-degeneracy.md) | English | [Project home](../../README_en.md) | [Algorithm overview](./README_en.md)

## Abstract

Degeneracy allows distinct bases to represent the same feasible point, so an algebraically valid basis exchange need not improve the objective. This chapter distinguishes combinatorial rules with exact-arithmetic termination guarantees from finite perturbations used to recover from floating-point stagnation.

## 1. Degeneracy and zero progress

A basic feasible solution in nonnegative standard form is degenerate when a basic coordinate is zero. For bounded variables, the corresponding definition is that a basic variable lies at a finite endpoint. A zero-ratio primal-simplex step may change the basis while leaving the point unchanged. Dual simplex may likewise exchange bases repeatedly on the same dual-objective level set.

Basis indices together with nonbasic bound states determine a basis state. Repeating that state under the same working objective indicates cycling along the iteration path; it does not establish optimality or infeasibility.

## 2. Bland's rule

For standard-form primal simplex, choose the smallest-index entering variable among those with negative reduced cost. Among variables tied at the exact minimum ratio, choose the smallest-index leaving variable. Bland's rule excludes cycling in exact arithmetic. Since there are finitely many bases, the feasible phase eventually reaches optimality or detects unboundedness.

Extending this principle to bounded variables requires consistent treatment of bound states and ties involving bound flips. Using smallest-index pricing for only some steps, or allowing candidates beyond the exact minimum ratio, does not directly justify the same finite-termination theorem.

## 3. Symbolic lexicographic ordering

Let $\varepsilon$ be a formal infinitesimal and write the right-hand side as

$$
b(\varepsilon)=b+\varepsilon v_1+\varepsilon^2v_2+\cdots.
$$

The primal ratio test then compares tuples

$$
\left(\frac{(B^{-1}b)_i}{d_i},
\frac{(B^{-1}v_1)_i}{d_i},
\frac{(B^{-1}v_2)_i}{d_i},\ldots\right)
$$

lexicographically, considering only $d_i>0$. Suitable independent perturbations remove ties in the formal problem, and a lexicographic potential strictly advances along valid iterations.

A dual rule can be constructed from symbolic costs

$$
c(\varepsilon)=c+\sum_{h\ge1}\varepsilon^h w_h
$$

For the current basis, the reduced cost at level $h$ is

$$
r_j^{(h)}=(w_h)_j-a_j^TB^{-T}(w_h)_\mathcal B.
$$

Dual breakpoints may therefore be compared first by their real ratios $s_jr_j/\beta_j$, and then by successive symbolic ratios $s_jr_j^{(h)}/\beta_j$. A formal infinitesimal cannot simply be replaced by a fixed small floating-point number. Floating-point evaluation of symbolic coefficients likewise does not guarantee exact comparisons for arbitrarily ill-conditioned problems.

## 4. Finite objective perturbations

Increase $c_j$ by $\delta_j\ge0$ for a nonbasic variable at its lower bound, and decrease it by $\delta_j$ at its upper bound. If basic costs are unchanged, the multiplier is unchanged while the directional reduced cost $s_jr_j$ increases by $\delta_j$, separating near-zero breakpoints. Fixed and free variables require separate treatment; endpoint signs do not apply to them unchanged.

Finite perturbations can reduce ties, but they optimize $\widetilde c=c+\Delta c$. Even a small $\|\Delta c\|$ may have an amplified objective effect on wide or unbounded coordinates. Cumulative perturbations must therefore be bounded. Once a primal feasible point is obtained, restore the original costs and perform primal-simplex refinement or feasibility repair under the original objective.

Bound perturbations similarly change the feasible set, so feasibility must be checked again after restoring the original bounds. Undoing perturbations is a necessary return to the original problem.

## 5. Stagnation measures and numerical recovery

The quantity

$$
M(x_\mathcal B)=\sum_i\operatorname{dist}
 (x_{\mathcal B_i},[\ell_{\mathcal B_i},u_{\mathcal B_i}])
$$

measures primal infeasibility but is not a globally monotone potential for the dual phase. A local decrease in $M$ on an objective level set does not exclude a later return to a previous state. Repetition history must be interpreted together with the current working objective.

The scale of dual-objective progress can be measured by $|\theta v_p|$. Changes below the scale of rounding and feasibility errors should be treated as numerically indistinguishable progress. This assessment selects recovery strategies; it does not replace KKT checks.

When the pricing row, entering direction, or basis equations become inconsistent, recompute and refine the basis equations before deciding on perturbations. Relaxing safe-pivot requirements or directly setting anomalous reduced costs to zero changes primal–dual relations and cannot establish optimality.
