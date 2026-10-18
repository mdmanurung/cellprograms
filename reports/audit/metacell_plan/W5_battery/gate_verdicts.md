# W5 SEACells-builder battery: gate verdicts

- G1 rate 0.0: 3/10 seeds with rare_abs >= 0.7x ceiling (mean ratio 0.658, range 0.573-0.760) -> FAIL (registered bar: >=8/10)
- G1 rate 0.1: 0/10 seeds with rare_abs >= 0.7x ceiling (mean ratio 0.092, range 0.031-0.233) -> FAIL (registered bar: >=8/10)
- G1 rate 0.3: 0/10 seeds with rare_abs >= 0.7x ceiling (mean ratio 0.053, range 0.024-0.191) -> FAIL (registered bar: >=8/10)
- G2 rate 0.0: mean matched cor 0.427 (range 0.393-0.475), 0/10 >= 0.8 -> FAIL (registered bar: >=0.8)
- G2 rate 0.1: mean matched cor 0.411 (range 0.345-0.445), 0/10 >= 0.8 -> FAIL (registered bar: >=0.8)
- G2 rate 0.3: mean matched cor 0.330 (range 0.283-0.371), 0/10 >= 0.8 -> FAIL (registered bar: >=0.8)
- G4 rate 0.0: strong-factor min projection cor 0.0172 (worst 0.0005), 0/10 >= 0.99 -> FAIL (registered bar: >=0.99 on strong factors)
- G4 rate 0.1: strong-factor min projection cor 0.0070 (worst 0.0009), 0/10 >= 0.99 -> FAIL (registered bar: >=0.99 on strong factors)
- G4 rate 0.3: strong-factor min projection cor 0.0077 (worst 0.0023), 0/10 >= 0.99 -> FAIL (registered bar: >=0.99 on strong factors)
- G3 (dubious <10% + ZeroRate <= k-means): mcRigor dubious rates are computed by the W1/W2 smoke path, not in this battery's JSONs; the battery records ZeroRate proxies only (n_mc_zero_counts, see per_combo_diagnostics.csv). Reported as NOT DIRECTLY MEASURED here; the smoke test covers the mcRigor QC path.
