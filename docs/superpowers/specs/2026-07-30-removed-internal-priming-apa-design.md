# Removed Internal Priming APA Output Design

## Goal

Preserve every APA site removed by the internal priming filter and write those
sites to a separate tab-delimited output with the A-content value that caused
the removal.

## Data Model

Add a `RemovedAPA` dataclass in `trek/internal_priming_filter.py` with parallel
lists:

- `site`: removed 1-based genomic positions.
- `count`: original supporting-read counts.
- `abundance`: original abundances from the unfiltered `TranscriptAPA`.
- `a_content`: downstream-window A-content values used by the filter.

The pipeline-level `removed_apa` object is
`Dict[str, RemovedAPA]`, keyed by transcript ID. It contains only transcripts
with at least one removed site. Removed-site abundances are not renormalized.

## Filter Interface

`InternalPrimingFilter.filter_apa_results()` returns a tuple:

```python
(filtered_apa, removed_apa)
```

`filtered_apa` retains its existing `Dict[str, TranscriptAPA]` shape.
`removed_apa` uses the new data model. The per-transcript filter returns both
the filtered `TranscriptAPA` and any removed site data so the A-content value is
captured during the same calculation that makes the filtering decision.

Single-site transcripts, transcripts missing from the annotation, dominant
sites, and non-dominant sites at or below the threshold do not produce removed
records. Existing filtering decisions remain unchanged.

## Pipeline Flow

When filtering is enabled, `ApaFinderPipeline` receives both dictionaries from
the filter. When `--no-filter-priming` is used, the filtered result remains the
raw APA result and `removed_apa` is an empty dictionary.

The pipeline always invokes the removed-result writer. This guarantees a
predictable file set whether filtering is enabled, disabled, or removes no
sites.

## Output File

Write removed sites to:

```text
{prefix}.internal_priming_removed.txt
```

The file uses one row per removed site and these columns:

```text
transcript_id
gene_id
gene_name
chromosome
strand
ID
site_position
site_count
site_abundance
transcript_biotype
a_content
```

`ID` uses the existing `chromosome:position:strand` format. Abundance and
A-content are formatted to four decimal places. If no sites were removed, the
file contains only the header. Missing transcript metadata is skipped in the
same way as the main result writer.

## Error Handling

The new result collection introduces no recoverable error state. File-system
errors continue to propagate through the pipeline's existing exception handler.
Parallel-list construction occurs from the same removed indices, keeping all
fields aligned.

## Testing

- Verify removed sites, counts, original abundances, and A-content values are
  returned together.
- Verify retained abundances are still renormalized while removed abundances
  remain original.
- Verify both genomic strands record their calculated A-content correctly.
- Verify the removed output header and row formatting.
- Verify filtering-disabled and no-removal runs create a header-only file.
- Verify the main APA result and summary outputs remain unchanged.
