# Internal Priming Filter Toggle Design

## Goal

Allow CLI and programmatic users to disable internal priming filtering while
keeping the existing default behavior enabled and preserving the 20 nt
downstream window.

## Interface

- Add `filter_priming: bool = True` to `ApaFinderPipeline`.
- Add the CLI flag `--no-filter-priming` using `argparse`'s `store_true`
  action.
- Pass `filter_priming=not args.no_filter_priming` when constructing the
  pipeline.

The positive pipeline property keeps the internal API readable even though the
CLI uses the conventional negative flag.

## Pipeline Behavior

`ApaFinderPipeline.run()` performs Step 5 only when `filter_priming` is true.
When it is false, the pipeline logs that internal priming filtering was skipped
and passes the unmodified APA results to the output stage.

The downstream sequence calculation, A-content threshold, dominant-site
behavior, and 20 nt default window remain unchanged.

## Documentation

Keep the documented `--no-filter-priming` option and correct the remaining
reference to an "upstream window" so all descriptions consistently say
downstream in transcript direction.

## Testing

- Verify filtering is enabled by default.
- Verify `--no-filter-priming` disables Step 5.
- Verify the CLI flag maps to the positive pipeline setting.
- Retain strand-aware sequence tests and add exact 20 nt boundary and threshold
  checks where needed.

Tests must run without invoking alignment or writing output files. Pipeline
stage tests will use a small subclass or patched stage methods to observe
whether the filtering stage is called.

## Error Handling

The toggle introduces no new error states. Existing filter errors continue to
propagate when filtering is enabled; disabling filtering avoids constructing or
running the filter.
