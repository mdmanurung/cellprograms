<div class="page-header">

# Coordination programs from an EB-CP fit

<div class="d-none name">

`eb_cp_programs.Rd`

</div>

</div>

<div class="ref-description section level2">

Extracts per-component cell-type memberships (a cell type is a member of component r iff its \|W\| reaches \`member_thresh\` of the component maximum) and donor activities. Components with fewer than \`min_members\` member cell types are private programs and are flagged, not dropped.

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
eb_cp_programs(cp, member_thresh = cp$member_thresh, min_members = 2L)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

cp<a href="#arg-cp" class="anchor" aria-label="anchor"></a><br>
An \`eb_cp_fit\` object from \[fit_eb_cp()\] or \[eb_cp_from_fit()\].

member_thresh<a href="#arg-member-thresh" class="anchor" aria-label="anchor"></a><br>
Membership threshold on max-normalized \|W\|.

min_members<a href="#arg-min-members" class="anchor" aria-label="anchor"></a><br>
Minimum member cell types for a "coordination" program.

</div>

<div class="section level2">

## Value<a href="#value" class="anchor" aria-label="anchor"></a>

A list with \`programs\` (one row per component: program id, member cell types, n_members, is_coordination, held-out z) and \`activity\` (donors x components, only for coordination components).

</div>
