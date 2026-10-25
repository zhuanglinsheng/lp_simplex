# State handoff for restricted-column models

[简体中文](./column-handoff.md) | English | [Documentation](./README_en.md) | [AUTO Scheduling](./auto-policy_en.md)

A cover master restricts columns of its parent LP: row order/types/RHS and selected-column coefficients/costs/bounds are unchanged; omitted variables are fixed at zero. It is not a public algorithm or permission to reuse a basis merely because dimensions match. `simplex_column_handoff.c/.h` expresses this exact relation. Arbitrary row changes, scaling, presolve and other related LPs require different contracts.

## Identity and relation

`simplex_ProblemIdentity` records the live input object, runtime context/model, revision, dimensions, nonzero count, and storage pointers for CSC/CSR, RHS, costs, bounds and row types. AUTO candidates, Dual seeds and checkpoints use it. Action parents and diagnostic journals are excluded. In-place edits still require `simplex_problem_mark_changed`; pointer checks are not whole-array checksums. Saturated revisions cannot authorize reuse; current-invocation audits remain separate from reusability.

`simplex_column_projection_bind` checks valid identities in the same runtime context, identical rows, an injective bidirectional column map (permutations allowed), exact selected coefficients/costs/effective bounds/variable types, identical row types/RHS, and zero admissibility for omitted variables. It also checks the child's CSR against CSC using caller-owned row cursors. This adapter requires CSR assembled in child-column order.

Binding checks the mathematical relation. Transfers recheck identities, revisions, storage and maps. Models/maps remain immutable while borrowed; both models must outlive the descriptor, which is retired before child destruction. Rejection preserves any previous descriptor; callers must inspect the return code rather than treating an earlier binding as a new success.

## Output permissions

| Handoff | Checks and output | Further audit |
| --- | --- | --- |
| Child basis to parent-coordinate carrier | Seed bound to this child; index ranges/uniqueness, BASIC set and status enums; structural ids mapped, logical ids become `parent.columns + row` | Index/status hint only, without factors, points, multipliers or certificates |
| Parent carrier to new master seed | Current parent identity; each old basic structural column still selected; reverse mapping and explicit target binding | Dual refactors, recomputes reduced costs/basic values and reclassifies feasibility |
| Child point to parent staging point | Child feasibility, zero padding and full-parent feasibility recomputation | Commit staging only on success; feasibility is not optimality |
| Child row multipliers to parent row candidate | Same row coordinates/sign convention; finite components | Full-parent pricing and certification; changing the column set changes dual feasibility |

Callers provide correctly sized, disjoint source/target basis arrays. Layout checks use caller-owned `ceil((parent.columns + parent.rows)/8)` bytes. Rejection preserves existing seed arrays, binding and validity. Point expansion uses separate staging storage, committed only after audit. A parent carrier borrows only its parent and survives child destruction; it cannot survive destruction of the parent.

`simplex_dual_seed_init` explicitly binds seed storage. Passing a seed to another model, even with equal dimensions, rejects import and solves cold without exporting into foreign-bound arrays. Related-model reuse requires checked projection/rebinding; an input revision also invalidates the old seed.

## Cover calls and resources

Each master receives a distinct shared-context `column-model` identity. `cover-master` records the actual Dual call and relative work grant, with scaling/backend calls nested beneath it. Column admission and tighter-tolerance refinement report different mechanism reasons. `column-handoff` records role/rejection without fabricating an LP result; a new master's basis cause references the actual basis-lift event that produced the parent carrier. The carrier retains a producer reference with history, context, token, kind/role and model. Restriction validates the same live context/event and completion, not only an index/history address. Reusing the same slot/kind/model cannot impersonate the old producer. Missing/incomplete/exhausted-token or different-history references become zero without changing identities, maps, numerical validity or work.

The fixed 256 MiB cover workspace limit now includes parent-point staging and the layout bitmap; it is not a total peak-RSS limit. Row cursors reuse master construction workspace, and old masters need not stay allocated. Full-column pricing, parent primal/dual checks, gap and rounding protection still certify the LP. A restricted master's infeasibility, budget or numerical failure is not a full-LP infeasibility proof.

### P06 complete flow exits

Admission permits nonnegative costs/matrix, nonnegative GE RHS and zero-lower LO/finite-BS variables. One individually sufficient-capacity column per row is a sufficient startup condition; its absence returns CondUnsatisfied, not infeasibility, because combinations may still suffice. Reset the invocation report at entry: unsupported input is not-applicable; refusal of four root workspaces is resource with zero work and no old output flags.

Each master uses the remaining grant and charges actual work; full-column pricing charges one more unit. Invalid work exhausts that remaining grant and loses classification permission. A nonoptimal master cannot classify the parent infeasible/unbounded; budget/resource retain their categories. Audit point, row-multiplier and basis mappings before committing. Later allocation failure may preserve an already-audited parent feasible point and consumed work; it exports no uncertified row multipliers.

A new master must add selected columns or strictly tighten dual tolerance after unsuccessful full-model certification, with at most three consecutive refinements. A failed master exits immediately; new objects/event numbers are not new consumed conditions. Success still requires full-column lower bounds, gap/roundoff and parent primal/dual audits; restriction alone proves none of these. Retire masters/projections before the next round and free all roots at exit. The 256 MiB bound covers fixed cover-owned workspace, not factories, factorization or total RSS.

## Verification

[test_column_handoff](../test/test_column_handoff.c) checks mathematical relations, invalid maps/identity, basis/point permissions, stale provenance and output preservation on rejection. [test_cover](../test/test_cover.c) exercises real column admission, accounting, owned allocation refusal and history-capacity independence. Delegated storage is tested at its owner; the local cover cap is not a global RSS guarantee. See [development history](./development/README.md).
