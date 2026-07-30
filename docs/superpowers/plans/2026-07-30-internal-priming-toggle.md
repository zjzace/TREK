# Internal Priming Filter Toggle Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a working `--no-filter-priming` option while retaining the default 20 nt downstream internal-priming filter.

**Architecture:** The pipeline owns a positive `filter_priming` setting and conditionally executes Step 5. The CLI exposes a conventional negative flag and maps it to that positive setting, so programmatic and CLI callers can both disable filtering.

**Tech Stack:** Python 3.12, `argparse`, `unittest`, `unittest.mock`, Biopython

---

## File Structure

- Modify `trek/run_trek.py`: own the toggle, conditionally execute Step 5, define the CLI flag, and map it into pipeline configuration.
- Create `trek/test_run_trek.py`: cover enabled/disabled pipeline behavior and CLI mapping without alignment or output writes.
- Modify `trek/test_internal_priming_filter.py`: characterize exact 20 nt boundaries and threshold behavior.
- Modify `README.md`: make downstream-window wording consistent with the implementation.

### Task 1: Make Pipeline Filtering Configurable

**Files:**
- Create: `trek/test_run_trek.py`
- Modify: `trek/run_trek.py:37-96`

- [ ] **Step 1: Write the failing disabled-pipeline test**

Create `trek/test_run_trek.py`:

```python
import tempfile
import unittest
from unittest.mock import patch

import run_trek


class RunTrekTest(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.temp_dir.cleanup()

    def _run_pipeline_with_stubbed_stages(self, pipeline):
        transcripts = {"tx1": object()}
        raw_results = {"tx1": object()}
        filtered_results = {"tx1": object()}

        with patch.object(
            pipeline, "_process_gtf", return_value=(transcripts, {}, {})
        ), patch.object(
            pipeline, "_load_assignments_if_valid", return_value={}
        ), patch.object(
            pipeline, "_find_apa_sites", return_value=raw_results
        ), patch.object(
            pipeline, "_filter_internal_priming", return_value=filtered_results
        ) as filter_stage, patch.object(
            pipeline, "_write_results"
        ) as writer:
            pipeline.run()

        return transcripts, raw_results, filtered_results, filter_stage, writer

    def test_pipeline_skips_internal_priming_when_disabled(self):
        pipeline = run_trek.ApaFinderPipeline(
            gtf_file="annotation.gtf",
            genome_fasta="genome.fa",
            fastq_files=[],
            output_dir=self.temp_dir.name,
            filter_priming=False,
        )

        transcripts, raw_results, _, filter_stage, writer = (
            self._run_pipeline_with_stubbed_stages(pipeline)
        )

        filter_stage.assert_not_called()
        writer.assert_called_once_with(raw_results, transcripts)
```

- [ ] **Step 2: Run the test and verify RED**

```bash
uv run --isolated --no-project \
  --with biopython --with interlap --with numpy --with scikit-learn \
  --with joblib --with tqdm --with pysam \
  -- python -m unittest discover -s trek -p 'test_run_trek.py' -v
```

Expected: ERROR with `unexpected keyword argument 'filter_priming'`.

- [ ] **Step 3: Add the pipeline setting and conditional stage**

Add to `ApaFinderPipeline.__init__()`:

```python
filter_priming=True,
```

Store it with the other filter settings:

```python
self.filter_priming = filter_priming
```

Replace the unconditional Step 5 block in `run()` with:

```python
if self.filter_priming:
    logger.info("STEP 5: Filtering internal priming artifacts")
    apa_results = self._filter_internal_priming(apa_results, transcripts)
else:
    logger.info("STEP 5: Internal priming filtering disabled; skipping")
```

- [ ] **Step 4: Run the focused test and verify GREEN**

Run the command from Step 2.

Expected: PASS.

- [ ] **Step 5: Add the default-enabled regression test**

Append to `RunTrekTest`:

```python
def test_pipeline_filters_internal_priming_by_default(self):
    pipeline = run_trek.ApaFinderPipeline(
        gtf_file="annotation.gtf",
        genome_fasta="genome.fa",
        fastq_files=[],
        output_dir=self.temp_dir.name,
    )

    transcripts, _, filtered_results, filter_stage, writer = (
        self._run_pipeline_with_stubbed_stages(pipeline)
    )

    filter_stage.assert_called_once()
    writer.assert_called_once_with(filtered_results, transcripts)
```

- [ ] **Step 6: Run both pipeline tests**

Run the command from Step 2.

Expected: both tests PASS.

- [ ] **Step 7: Commit configurable pipeline behavior**

```bash
git add trek/test_run_trek.py trek/run_trek.py
git commit -m "feat: make internal priming filtering configurable"
```

### Task 2: Expose the CLI Flag

**Files:**
- Modify: `trek/test_run_trek.py`
- Modify: `trek/run_trek.py:299-345`

- [ ] **Step 1: Write the failing CLI mapping test**

Add `import sys` and append to `RunTrekTest`:

```python
def test_no_filter_priming_flag_disables_pipeline_filter(self):
    argv = [
        "trek", "-g", "annotation.gtf", "-f", "genome.fa",
        "-q", "reads.fastq", "--no-filter-priming",
    ]

    with patch.object(sys, "argv", argv), \
         patch.object(run_trek.Path, "exists", return_value=True), \
         patch.object(run_trek, "ApaFinderPipeline") as pipeline_class:
        result = run_trek.main()

    self.assertEqual(result, 0)
    self.assertFalse(pipeline_class.call_args.kwargs["filter_priming"])
    pipeline_class.return_value.run.assert_called_once_with()
```

- [ ] **Step 2: Run the test and verify RED**

Run the Task 1 test command.

Expected: ERROR with `unrecognized arguments: --no-filter-priming`.

- [ ] **Step 3: Add the CLI flag and mapping**

Add to `parse_arguments()`:

```python
parser.add_argument(
    '--no-filter-priming',
    action='store_true',
    help='Disable internal priming filtering (enabled by default)'
)
```

Add to the `ApaFinderPipeline(...)` call in `main()`:

```python
filter_priming=not args.no_filter_priming,
```

- [ ] **Step 4: Run the focused tests and verify GREEN**

Run the Task 1 test command.

Expected: all three tests PASS.

- [ ] **Step 5: Commit the CLI flag**

```bash
git add trek/test_run_trek.py trek/run_trek.py
git commit -m "feat: add internal priming filter CLI toggle"
```

### Task 3: Lock Down the 20 nt Filter Contract

**Files:**
- Modify: `trek/test_internal_priming_filter.py`

- [ ] **Step 1: Add exact-boundary characterization tests**

For `position = 61` and `window_size = 20`, add tests with A/T only at the
excluded site base and excluded 21st downstream base. Assert the result is
`0.0` on both strands. Retain the existing tests that assert `1.0` when all 20
included downstream bases are A/T.

Use these zero-based indices:

```python
# Plus included window [61:81], excluded site 60, excluded next base 81.
# Minus included window [40:60], excluded site 60, excluded next base 39.
```

- [ ] **Step 2: Add threshold characterization tests**

Import `TranscriptAPA` from `apa_finder`. Build this secondary-site case:

```python
apa = TranscriptAPA(
    site=[31, 61],
    count=[20, 10],
    abundance=[2 / 3, 1 / 3],
)
transcript = SimpleNamespace(
    chromosome=self.chromosome,
    strand="+",
    transcript_id="tx1",
)
```

Set `filter_obj.a_content_threshold = 0.5`. With exactly 10 As in genomic slice
`[61:81]`, assert site 61 remains; with 11 As, assert site 61 is removed.

- [ ] **Step 3: Run the characterization tests**

```bash
uv run --isolated --no-project \
  --with biopython --with interlap --with numpy --with scikit-learn \
  --with joblib --with tqdm \
  -- python -m unittest discover -s trek \
  -p 'test_internal_priming_filter.py' -v
```

Expected: all tests PASS because they characterize behavior already present in
the corrected working tree.

- [ ] **Step 4: Commit filter contract tests**

```bash
git add trek/test_internal_priming_filter.py
git commit -m "test: cover internal priming filter boundaries"
```

### Task 4: Correct Documentation and Run Full Verification

**Files:**
- Modify: `README.md:66-81`

- [ ] **Step 1: Correct the remaining directional wording**

Change:

```markdown
3. Calculate A proportion in the upstream window:
```

to:

```markdown
3. Calculate A proportion in the downstream transcript-direction window:
```

Keep the documented 20 nt default and `--no-filter-priming` option.

- [ ] **Step 2: Run the complete unit suite**

```bash
uv run --isolated --no-project \
  --with biopython --with interlap --with numpy --with scikit-learn \
  --with joblib --with tqdm --with pysam \
  -- python -m unittest discover -s trek -p 'test*.py' -v
```

Expected: all tests PASS.

- [ ] **Step 3: Run static checks**

```bash
python3 -m compileall -q trek
git diff --check
```

Expected: both commands exit 0 with no errors.

- [ ] **Step 4: Inspect the final scope**

```bash
git status --short
git diff -- README.md trek/run_trek.py trek/internal_priming_filter.py \
  trek/test_internal_priming_filter.py trek/test_run_trek.py
```

Expected: the 20 nt downstream filter is preserved, the toggle is implemented,
and unrelated untracked files remain untouched.

- [ ] **Step 5: Commit documentation**

```bash
git add README.md
git commit -m "docs: clarify downstream priming window"
```
