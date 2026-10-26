# Empty rows, unresolved cells, and independent summaries

All three features are opt-in and work with Plotly and Matplotlib, including
expanded mutation-type/numeric tracks. Run the complete synthetic example:

	source .venv/bin/activate
	PYTHONPATH=src python docs/examples/cell_states.py

It writes .venv/issue-1-smoke/cell-states.png and cell-states.html. No image-export
extra is required. The PNG uses Matplotlib; the HTML uses Plotly.

## Preserve explicit rows

Set `preserve_empty_genes=True` together with a nonempty, unique
`include_genes` sequence. Requested genes remain in that exact order,
even if every selected gene has zero events; `ignore_genes` still removes
explicit exclusions. This order overrides pathway sorting; pathway decorations
are split into contiguous runs. An empty mutation DataFrame with the required
columns is accepted in this mode. Supply samples using a state table or metadata
(with `metadata_require_mutations=False, show_all_samples=True`) for a
cohort with no events. Ordinary sample-selection rules still apply: use
`show_all_samples=True` to retain samples without selected mutations.
Without this option the original row-selection behavior is unchanged.

## Non-mutation states

`cell_states` is a separate DataFrame with the same gene and sample column
names as the mutation data plus a fixed `state` column. Example:

| gene | sample | state |
| --- | --- | --- |
| TP53 | S2 | unresolved |
| EMPTY | S3 | not_assessed |

Both standard states have distinct default colors and explicit "Cell state: ..."
legend entries. Matplotlib additionally hatches them. Override their colors or
add named states with `cell_state_palette={"unresolved": "#E69F00", ...}`.
All such states are non-events, not mutation categories. There must be exactly
one state per pair; missing identifiers, unknown states without a palette,
duplicates, and mutation/state collisions are errors. Inputs are not modified.

States overlay every track for their gene/sample, including numeric tracks;
they never become mutation/variant tiles. State-only samples for selected genes
are included automatically. Sample filters restrict state samples, too; mutation
filters apply only to mutation events (not the state table). Conflicts are checked
against retained events after filtering. Rows not selected for display are ignored.
States never affect event ranking, mutation counts, mutation-type summaries,
variant aggregation, or computed TMB. Computed TMB remains a count of recorded
mutation rows, not a claim that unrecorded/unassessed variants are negative;
state-only samples therefore have zero recorded events. Custom TMB is caller data.

Unmarked background cells retain the existing no-recorded-event appearance.
They are not evidence of assay negativity or completeness. If that distinction
matters, supply the relevant unresolved/not-assessed observations explicitly.

## Denominators are explicit, never inferred from unresolved cells

Supplying `cell_states` disables implicit cohort denominators for row
percentages (even for an empty state table). Count bars work without a denominator,
and hover then reports "percentage unavailable" rather than treating unknowns
as negatives. To request percentage bars or percentage labels, supply a positive,
finite denominator for every displayed gene via `gene_denominators`:

	oncoplot(events, ..., cell_states=states,
	         gene_denominators={"TP53": 12, "KRAS": 15, "EMPTY": 8},
	         draw_gene_bar=True,
	         options={"gene_bar_mode": "percent", "show_gene_bar_labels": True})

Declare the scientifically appropriate assessed/eligible population per gene;
unresolved/not-assessed subjects should be excluded from an assessed-only
denominator. Nothing automatically subtracts subjects or assumes an assay panel.
Explicit denominators may also be used without states. Counts are distinct
mutated displayed samples per gene (multi-hit cells count once). Percentage bars
are 100 * each mutation-type count / declared gene denominator; total percentage
labels are 100 * total mutated samples / that denominator. These are prevalence,
not composition bars. Declared denominators may refer to a wider population than
the visible subjects; filtering does not rescale them. They must not be smaller
than the numerator. Without these opt-ins, existing total_samples behavior and
mutation-type composition percentages are unchanged.

## Independently supplied aggregates

`row_summaries` is a DataFrame with the configured gene column,
`count`, and optionally `denominator`. Identifiers align by name,
not position. Supply each displayed gene exactly once; extra non-displayed genes
are ignored. Missing rows are errors, not assumed zero. Counts must be finite,
nonnegative integers. A supplied denominator must be finite, positive, and at
least the count. Every displayed row needs a denominator for percentage bars or
percentage labels; count-only displays do not need one. Supply external
denominators in this table, not simultaneously through gene_denominators.

External bars and their subplot title are labeled "External summary", including
when axis labels are hidden. Their percentages use each supplied denominator.
External summaries replace only row-bar totals; they are not stacked into
fabricated mutation types and do not create subject-level cells, reorder genes,
or affect TMB. The caller owns the meaning/population of these independent counts.

Both tables support the standard public parameter-dictionary/YAML loading forms:
a DataFrame, a CSV/TSV path, or a mapping with path and read_csv arguments. YAML
relative paths resolve against the config file. All five options are available
through oncoplot keyword arguments or params; the transformation options are also
accepted by prepare_oncoplot_data. See the runnable Python example above.
