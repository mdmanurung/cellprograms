<div class="page-header">

# CPUB local candidate: support and migration

</div>

<div id="cpub-local-candidate-support-and-migration" class="section level1">

This is a local implementation candidate, not a scientific release certificate. Ordinary donor-level EBMF retains the default normal activity / point-Laplace loading estimator inherited from flashier/ebnm.

<div class="section level2">

## Intentional fail-closed changes<a href="#intentional-fail-closed-changes" class="anchor" aria-label="anchor"></a>

- Residualized projection requires trained design state. Refit old serialized residualized objects; do not synthesize missing state by rerunning the query formula. Custom covariate functions require an explicitly trained preprocessing recipe outside the API. Supported formula transforms are linear/factor terms, interactions, `poly`, `scale(variable)`, <a href="https://rdrr.io/r/splines/ns.html" class="external-link"><code>splines::ns</code></a>, <a href="https://rdrr.io/r/splines/bs.html" class="external-link"><code>splines::bs</code></a>, `I`, arithmetic, `log`, `log1p`, and `sqrt`.
- Preprocessing specs now carry units and a verified digest. Intact 1.0 specs retain unknown units. A digest mismatch requires deliberate refreezing.
- Numeric `var_floor` is **fixed known residual SD**, not a lower bound on an estimated SD. `var_floor="vignette"` retains its historical deterministic small-noise recipe; using it does not establish calibration.
- Optional torch uses a different estimator and requires the explicit settings in NEWS. Its real backend integration is not established by pure loss-trace tests or mocked calls.
- Score/private-first block fitting is undefined and rejected. Block membership, EB-CP memberships and uncalibrated loading sharing are descriptive.
- CMP missing-support guards preserve NA and expose callability; permutation validity with selected models or arbitrary missingness is not established by fixing singular QR. Node-family adjustment does not control false edges or merged CMP-set membership.

</div>

<div class="section level2">

## Release prerequisites still open<a href="#release-prerequisites-still-open" class="anchor" aria-label="anchor"></a>

DESCRIPTION reuses the maintainer’s existing public repository contact under delegated repository-completion approval. Optional heavy dependencies are not installed merely to remove skips. Cross-platform CI must be judged from actual run receipts, not from configuration alone. Data permissions, final manifests, independent scientific campaigns and release publication remain separate gates.

</div>

</div>
