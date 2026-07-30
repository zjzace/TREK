# Removed Internal Priming APA Output Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Preserve internally primed APA sites in a structured `removed_apa` object and always write them to a separate tab-delimited result file.

**Architecture:** A `RemovedAPA` dataclass carries removed site coordinates, original counts and abundances, and calculated A-content. The filter returns retained and removed dictionaries together; the pipeline propagates both and uses a dedicated writer that also creates a header-only file when no sites were removed.

**Tech Stack:** Python 3.12, dataclasses, type hints, `unittest`, Biopython

---

## File Structure

- Modify `trek/internal_priming_filter.py`: define `RemovedAPA` and collect removed sites during filtering.
- Modify `trek/run_trek.py`: propagate `removed_apa` and write `{prefix}.internal_priming_removed.txt`.
- Modify `trek/test_internal_priming_filter.py`: verify removed data, original abundance, A-content, and strand handling.
- Modify `trek/test_run_trek.py`: verify pipeline propagation, removed-file rows, and header-only behavior.
- Modify `README.md`: document the new output file.

### Task 1: Add the RemovedAPA Model and Output Writer

**Files:**
- Modify: `trek/internal_priming_filter.py:7-16`
- Modify: `trek/run_trek.py:207-254`
- Modify: `trek/test_run_trek.py`

- [ ] **Step 1: Write the failing removed-output row test**

Add imports:

```python
from types import SimpleNamespace

from internal_priming_filter import RemovedAPA
```

Add this test to `ApaFinderPipelineTest`:

```python
def test_write_removed_internal_priming_results(self):
    transcript = SimpleNamespace(
        ncbi_gene_id="gene1",
        gene_name="GENE1",
        chromosome="chr1",
        strand="+",
        transcript_biotype="protein_coding",
    )
    removed_apa = {
        "tx1": RemovedAPA(
            site=[61],
            count=[10],
            abundance=[1 / 3],
            a_content=[0.55],
        )
    }

    with tempfile.TemporaryDirectory() as output_dir:
        pipeline = ApaFinderPipeline(
            gtf_file="annotation.gtf",
            genome_fasta="genome.fa",
            fastq_files=[],
            output_dir=output_dir,
            prefix="sample",
        )
        pipeline._write_removed_apa_results(removed_apa, {"tx1": transcript})

        output = (
            run_trek.Path(output_dir) / "sample.internal_priming_removed.txt"
        ).read_text()

    self.assertEqual(
        output,
        "transcript_id\tgene_id\tgene_name\tchromosome\tstrand\tID\t"
        "site_position\tsite_count\tsite_abundance\ttranscript_biotype\t"
        "a_content\n"
        "tx1\tgene1\tGENE1\tchr1\t+\tchr1:61:+\t61\t10\t0.3333\t"
        "protein_coding\t0.5500\n",
    )
```

- [ ] **Step 2: Run the focused test and verify RED**

```bash
uv run --python 3.12 --isolated --no-project \
  --with biopython --with interlap --with numpy --with scikit-learn \
  --with joblib --with tqdm --with pysam \
  -- python -m unittest discover -s trek -p 'test_run_trek.py' -v
```

Expected: import error because `RemovedAPA` does not exist.

- [ ] **Step 3: Add the data model and dedicated writer**

In `trek/internal_priming_filter.py`, import `dataclass` and `List`, then add:

```python
@dataclass
class RemovedAPA:
    """APA sites removed as internal priming artifacts."""

    site: List[int]
    count: List[int]
    abundance: List[float]
    a_content: List[float]
```

In `ApaFinderPipeline`, add:

```python
def _write_removed_apa_results(self, removed_apa, transcripts):
    """Write APA sites removed by the internal priming filter."""
    output_file = self.output_dir / f"{self.prefix}.internal_priming_removed.txt"

    with open(output_file, 'w') as f:
        f.write(
            "transcript_id\tgene_id\tgene_name\tchromosome\tstrand\t"
            "ID\tsite_position\tsite_count\tsite_abundance\t"
            "transcript_biotype\ta_content\n"
        )

        for transcript_id, apa in removed_apa.items():
            transcript = transcripts.get(transcript_id)
            if not transcript:
                continue

            for position, count, abundance, a_content in zip(
                apa.site, apa.count, apa.abundance, apa.a_content
            ):
                locus_id = (
                    f"{transcript.chromosome}:{position}:{transcript.strand}"
                )
                f.write(
                    f"{transcript_id}\t{transcript.ncbi_gene_id}\t"
                    f"{transcript.gene_name}\t{transcript.chromosome}\t"
                    f"{transcript.strand}\t{locus_id}\t{position}\t{count}\t"
                    f"{abundance:.4f}\t{transcript.transcript_biotype}\t"
                    f"{a_content:.4f}\n"
                )

    logger.info(f"Saved removed internal priming sites: {output_file}")
```

- [ ] **Step 4: Run the focused test and verify GREEN**

Run the command from Step 2.

Expected: the new row test and all existing pipeline tests PASS.

- [ ] **Step 5: Add and run the header-only writer test**

```python
def test_write_empty_removed_internal_priming_results(self):
    with tempfile.TemporaryDirectory() as output_dir:
        pipeline = ApaFinderPipeline(
            gtf_file="annotation.gtf",
            genome_fasta="genome.fa",
            fastq_files=[],
            output_dir=output_dir,
            prefix="sample",
        )
        pipeline._write_removed_apa_results({}, {})

        lines = (
            run_trek.Path(output_dir) / "sample.internal_priming_removed.txt"
        ).read_text().splitlines()

    self.assertEqual(len(lines), 1)
    self.assertTrue(lines[0].endswith("\ta_content"))
```

Run the command from Step 2.

Expected: all tests PASS.

- [ ] **Step 6: Commit the model and writer**

```bash
git add trek/internal_priming_filter.py trek/run_trek.py trek/test_run_trek.py
git commit -m "feat: add removed internal priming output writer"
```

### Task 2: Collect and Propagate Removed APA Sites

**Files:**
- Modify: `trek/internal_priming_filter.py:40-147`
- Modify: `trek/run_trek.py:70-205`
- Modify: `trek/test_internal_priming_filter.py`
- Modify: `trek/test_run_trek.py`

- [ ] **Step 1: Write the failing filter result test**

Add this test to `InternalPrimingFilterTest`:

```python
def test_filter_returns_removed_site_with_original_values(self):
    sequence = list("C" * 120)
    sequence[self.position:self.position + 11] = "A" * 11
    filter_obj = self._filter_for_sequence(sequence)
    filter_obj.a_content_threshold = 0.5
    apa = TranscriptAPA(
        site=[31, self.position],
        count=[20, 10],
        abundance=[2 / 3, 1 / 3],
    )
    transcript = SimpleNamespace(
        chromosome=self.chromosome,
        strand="+",
        transcript_id="tx1",
    )

    filtered, removed = filter_obj.filter_apa_results(
        {"tx1": apa}, {"tx1": transcript}
    )

    self.assertEqual(filtered["tx1"].site, [31])
    self.assertEqual(filtered["tx1"].abundance, [1.0])
    self.assertEqual(removed["tx1"].site, [self.position])
    self.assertEqual(removed["tx1"].count, [10])
    self.assertEqual(removed["tx1"].abundance, [1 / 3])
    self.assertEqual(removed["tx1"].a_content, [0.55])
```

In the two existing threshold tests, replace the single return assignment with:

```python
filtered, removed = filter_obj._filter_transcript_apa(
    TranscriptAPA(
        site=[31, self.position],
        count=[20, 10],
        abundance=[2 / 3, 1 / 3],
    ),
    SimpleNamespace(
        chromosome=self.chromosome,
        strand="+",
        transcript_id="tx1",
    ),
)
```

In the half-A test, assert:

```python
self.assertEqual(filtered.site, [31, self.position])
self.assertEqual(removed.site, [])
```

In the above-half test, assert:

```python
self.assertEqual(filtered.site, [31])
self.assertEqual(removed.site, [self.position])
```

- [ ] **Step 2: Run the filter tests and verify RED**

```bash
uv run --python 3.12 --isolated --no-project \
  --with biopython --with interlap --with numpy --with scikit-learn \
  --with joblib --with tqdm \
  -- python -m unittest discover -s trek \
  -p 'test_internal_priming_filter.py' -v
```

Expected: unpacking fails because the filter currently returns only retained
results.

- [ ] **Step 3: Return retained and removed data from the filter**

Import `Tuple` and update the public and per-transcript return annotations:

```python
def filter_apa_results(
    self,
    apa_results: Dict[str, TranscriptAPA],
    transcripts: Dict[str, Transcript],
) -> Tuple[Dict[str, TranscriptAPA], Dict[str, RemovedAPA]]:
```

```python
def _filter_transcript_apa(
    self, apa: TranscriptAPA, transcript: Transcript
) -> Tuple[TranscriptAPA, RemovedAPA]:
```

In `_filter_transcript_apa`, initialize and populate the two removed-data
collections in the existing decision loop:

```python
keep_indices = [0]
removed_indices = []
removed_a_content = []

for idx in range(1, len(apa.site)):
    position = apa.site[idx]
    a_proportion = self._calculate_a_content(
        chromosome=transcript.chromosome,
        position=position,
        strand=transcript.strand,
    )

    if a_proportion <= self.a_content_threshold:
        keep_indices.append(idx)
    else:
        removed_indices.append(idx)
        removed_a_content.append(a_proportion)
```

Keep the existing retained-site construction and abundance normalization, then
return:

```python
return (
    TranscriptAPA(
        site=filtered_sites,
        count=filtered_counts,
        abundance=filtered_abundances,
    ),
    RemovedAPA(
        site=[apa.site[i] for i in removed_indices],
        count=[apa.count[i] for i in removed_indices],
        abundance=[apa.abundance[i] for i in removed_indices],
        a_content=removed_a_content,
    ),
)
```

In `filter_apa_results`, initialize `removed_results` beside
`filtered_results`, then replace the multiple-site filtering assignment with:

```python
filtered_apa, removed_apa = self._filter_transcript_apa(apa, transcript)
filtered_results[transcript_id] = filtered_apa
if removed_apa.site:
    removed_results[transcript_id] = removed_apa
```

Keep the existing statistics based on `filtered_apa`, then return:

```python
return filtered_results, removed_results
```

- [ ] **Step 4: Add the negative-strand removed A-content test**

```python
def test_minus_strand_removed_site_records_t_as_a_content(self):
    sequence = list("C" * 120)
    sequence[self.position - self.window_size - 1:self.position - 10] = "T" * 11
    filter_obj = self._filter_for_sequence(sequence)
    filter_obj.a_content_threshold = 0.5

    _, removed = filter_obj._filter_transcript_apa(
        TranscriptAPA(
            site=[91, self.position],
            count=[20, 10],
            abundance=[2 / 3, 1 / 3],
        ),
        SimpleNamespace(
            chromosome=self.chromosome,
            strand="-",
            transcript_id="tx1",
        ),
    )

    self.assertEqual(removed.site, [self.position])
    self.assertEqual(removed.a_content, [0.55])
```

Run the command from Step 2.

Expected: all filter tests PASS.

- [ ] **Step 5: Write failing pipeline propagation assertions**

In the enabled pipeline test, define and configure:

```python
removed_results = {"transcript": object()}
pipeline._filter_internal_priming = Mock(
    return_value=(filtered_results, removed_results)
)
pipeline._write_removed_apa_results = Mock()
```

Enabled assertion:

```python
pipeline._write_removed_apa_results.assert_called_once_with(
    removed_results, transcripts
)
```

Disabled assertion:

```python
pipeline._write_removed_apa_results.assert_called_once_with({}, transcripts)
```

In the disabled test, add only
`pipeline._write_removed_apa_results = Mock()` before `pipeline.run()`; its
filter mock remains unused and its existing raw-result assertion remains.

- [ ] **Step 6: Run pipeline tests and verify RED**

Run the Task 1 test command.

Expected: failures because `run()` neither unpacks the filter tuple nor calls
the removed-result writer.

- [ ] **Step 7: Propagate removed data through the pipeline**

In `run()`, initialize and propagate removed data:

```python
removed_apa = {}
if self.filter_priming:
    logger.info("STEP 5: Filtering internal priming artifacts")
    apa_results, removed_apa = self._filter_internal_priming(
        apa_results, transcripts
    )
else:
    logger.info("STEP 5: Internal priming filtering disabled; skipping")

logger.info("STEP 6: Writing results")
self._write_results(apa_results, transcripts)
self._write_removed_apa_results(removed_apa, transcripts)
```

`_filter_internal_priming()` already forwards the filter return value, so no
additional adapter is needed.

- [ ] **Step 8: Run filter and pipeline tests and verify GREEN**

Run both focused commands from Steps 2 and 6.

Expected: all focused tests PASS.

- [ ] **Step 9: Commit collection and propagation**

```bash
git add trek/internal_priming_filter.py trek/run_trek.py \
  trek/test_internal_priming_filter.py trek/test_run_trek.py
git commit -m "feat: retain removed internal priming APA sites"
```

### Task 3: Document and Verify the New Output

**Files:**
- Modify: `README.md:83-88`
- Modify: `trek/test_run_trek.py`

- [ ] **Step 1: Add a main-output regression test**

Import `TranscriptAPA` from `apa_finder`, then add:

```python
def test_main_output_format_remains_unchanged(self):
    transcript = SimpleNamespace(
        ncbi_gene_id="gene1",
        gene_name="GENE1",
        chromosome="chr1",
        strand="+",
        transcript_biotype="protein_coding",
    )
    apa_results = {
        "tx1": TranscriptAPA(site=[31], count=[20], abundance=[1.0])
    }

    with tempfile.TemporaryDirectory() as output_dir:
        pipeline = ApaFinderPipeline(
            gtf_file="annotation.gtf",
            genome_fasta="genome.fa",
            fastq_files=[],
            output_dir=output_dir,
            prefix="sample",
        )
        pipeline._write_results(apa_results, {"tx1": transcript})

        main_output = (
            run_trek.Path(output_dir) / "sample.apa_sites.txt"
        ).read_text()
        summary = (
            run_trek.Path(output_dir) / "sample.summary.txt"
        ).read_text()

    self.assertEqual(
        main_output,
        "transcript_id\tgene_id\tgene_name\tchromosome\tstrand\tID\t"
        "site_position\tsite_count\tsite_abundance\ttranscript_biotype\n"
        "tx1\tgene1\tGENE1\tchr1\t+\tchr1:31:+\t31\t20\t1.0000\t"
        "protein_coding\n",
    )
    self.assertIn("Total transcripts analyzed: 1", summary)
    self.assertIn("Transcripts with alternative TES: 0", summary)
```

- [ ] **Step 2: Run pipeline tests**

Run the Task 1 test command.

Expected: all tests PASS because the main output implementation remains
unchanged.

- [ ] **Step 3: Document the output file**

Add to the root README output list:

```markdown
- `{prefix}.internal_priming_removed.txt`: APA sites removed as internal
  priming artifacts, including downstream-window A-content
```

State that this file contains only its header when filtering is disabled or no
sites exceed the threshold.

- [ ] **Step 4: Run the complete Python 3.12 test suite**

```bash
uv run --python 3.12 --isolated --no-project \
  --with biopython --with interlap --with numpy --with scikit-learn \
  --with joblib --with tqdm --with pysam \
  -- python -m unittest discover -s trek -p 'test*.py' -v
```

Expected: all tests PASS.

- [ ] **Step 5: Run static and scope checks**

```bash
python3 -m compileall -q trek
git diff --check
git status --short
```

Expected: compile and diff checks exit 0; only intended tracked files and the
user's pre-existing unrelated untracked files appear.

- [ ] **Step 6: Commit documentation and regression coverage**

```bash
git add README.md trek/test_run_trek.py
git commit -m "docs: describe removed internal priming output"
```
