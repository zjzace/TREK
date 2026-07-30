import tempfile
import unittest
from unittest.mock import Mock

from run_trek import ApaFinderPipeline


class ApaFinderPipelineTest(unittest.TestCase):
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

            pipeline.run()

            pipeline._filter_internal_priming.assert_not_called()
            pipeline._write_results.assert_called_once_with(raw_results, transcripts)

    def test_run_filters_internal_priming_by_default(self):
        transcripts = {"transcript": object()}
        raw_results = {"transcript": object()}
        filtered_results = {"transcript": object()}

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
            pipeline._filter_internal_priming = Mock(return_value=filtered_results)
            pipeline._write_results = Mock()

            pipeline.run()

            pipeline._filter_internal_priming.assert_called_once_with(
                raw_results, transcripts
            )
            pipeline._write_results.assert_called_once_with(filtered_results, transcripts)


if __name__ == "__main__":
    unittest.main()
