import sys
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import Mock
from unittest.mock import patch

import run_trek
from apa_finder import TranscriptAPA
from internal_priming_filter import RemovedAPA
from run_trek import ApaFinderPipeline


class ApaFinderPipelineTest(unittest.TestCase):
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
                run_trek.Path(output_dir)
                / "sample.internal_priming_removed.txt"
            ).read_text()

        self.assertEqual(
            output,
            "transcript_id\tgene_id\tgene_name\tchromosome\tstrand\tID\t"
            "site_position\tsite_count\tsite_abundance\ttranscript_biotype\t"
            "a_content\n"
            "tx1\tgene1\tGENE1\tchr1\t+\tchr1:61:+\t61\t10\t0.3333\t"
            "protein_coding\t0.5500\n",
        )

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
                run_trek.Path(output_dir)
                / "sample.internal_priming_removed.txt"
            ).read_text().splitlines()

        self.assertEqual(len(lines), 1)
        self.assertTrue(lines[0].endswith("\ta_content"))

    def test_no_filter_priming_flag_disables_pipeline_filter(self):
        argv = [
            "trek",
            "-g", "annotation.gtf",
            "-f", "genome.fa",
            "-q", "reads.fastq",
            "--no-filter-priming",
        ]

        with patch.object(sys, "argv", argv), \
             patch.object(run_trek.Path, "exists", return_value=True), \
             patch.object(run_trek, "ApaFinderPipeline") as pipeline_class:
            result = run_trek.main()

        self.assertEqual(result, 0)
        self.assertFalse(pipeline_class.call_args.kwargs["filter_priming"])
        pipeline_class.return_value.run.assert_called_once_with()

    def test_run_skips_internal_priming_filter_when_disabled(self):
        transcripts = {"transcript": object()}
        raw_results = {"transcript": object()}
        filtered_results = {"transcript": object()}

        with tempfile.TemporaryDirectory() as output_dir:
            pipeline = ApaFinderPipeline(
                gtf_file="annotation.gtf",
                genome_fasta="genome.fa",
                fastq_files=["reads.fastq"],
                output_dir=output_dir,
                filter_priming=False,
            )
            pipeline._process_gtf = Mock(return_value=(transcripts, {}, {}))
            pipeline._load_assignments_if_valid = Mock(return_value={})
            pipeline._find_apa_sites = Mock(return_value=raw_results)
            pipeline._filter_internal_priming = Mock(return_value=filtered_results)
            pipeline._write_results = Mock()
            pipeline._write_removed_apa_results = Mock()

            pipeline.run()

            pipeline._filter_internal_priming.assert_not_called()
            pipeline._write_results.assert_called_once_with(raw_results, transcripts)
            pipeline._write_removed_apa_results.assert_called_once_with(
                {}, transcripts
            )

    def test_run_filters_internal_priming_by_default(self):
        transcripts = {"transcript": object()}
        raw_results = {"transcript": object()}
        filtered_results = {"transcript": object()}
        removed_results = {"transcript": object()}

        with tempfile.TemporaryDirectory() as output_dir:
            pipeline = ApaFinderPipeline(
                gtf_file="annotation.gtf",
                genome_fasta="genome.fa",
                fastq_files=["reads.fastq"],
                output_dir=output_dir,
            )
            pipeline._process_gtf = Mock(return_value=(transcripts, {}, {}))
            pipeline._load_assignments_if_valid = Mock(return_value={})
            pipeline._find_apa_sites = Mock(return_value=raw_results)
            pipeline._filter_internal_priming = Mock(
                return_value=(filtered_results, removed_results)
            )
            pipeline._write_results = Mock()
            pipeline._write_removed_apa_results = Mock()

            pipeline.run()

            pipeline._filter_internal_priming.assert_called_once_with(
                raw_results, transcripts
            )
            pipeline._write_results.assert_called_once_with(filtered_results, transcripts)
            pipeline._write_removed_apa_results.assert_called_once_with(
                removed_results, transcripts
            )


if __name__ == "__main__":
    unittest.main()
