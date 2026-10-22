# Coordination programs from an EB-CP fit

Extracts per-component cell-type memberships (a cell type is a member of
component r iff its \|W\| reaches \`member_thresh\` of the component
maximum) and donor activities. Components with fewer than
\`min_members\` member cell types are private programs and are flagged,
not dropped.

## Usage

``` r
eb_cp_programs(cp, member_thresh = cp$member_thresh, min_members = 2L)
```

## Arguments

- cp:

  An \`eb_cp_fit\` object from \[fit_eb_cp()\] or \[eb_cp_from_fit()\].

- member_thresh:

  Membership threshold on max-normalized \|W\|.

- min_members:

  Minimum member cell types for a "coordination" program.

## Value

A list with \`programs\` (one row per component: program id, member cell
types, n_members, is_coordination, held-out z) and \`activity\` (donors
x components, only for coordination components).
