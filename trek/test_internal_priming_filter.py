import unittest
from types import SimpleNamespace

from Bio.Seq import Seq

from apa_finder import TranscriptAPA
from internal_priming_filter import InternalPrimingFilter


class InternalPrimingFilterTest(unittest.TestCase):
    chromosome = "chr1"
    position = 61
    window_size = 20

    def _filter_for_sequence(self, sequence):
        filter_obj = InternalPrimingFilter.__new__(InternalPrimingFilter)
        filter_obj.window_size = self.window_size
        filter_obj.genome_seqs = {
            self.chromosome: SimpleNamespace(seq=Seq("".join(sequence)))
        }
        return filter_obj

    def test_plus_strand_uses_downstream_sequence(self):
        sequence = list("C" * 120)
        sequence[self.position - self.window_size - 1:self.position - 1] = "A" * self.window_size
        filter_obj = self._filter_for_sequence(sequence)

        a_content = filter_obj._calculate_a_content(
            self.chromosome, self.position, "+"
        )

        self.assertEqual(a_content, 0.0)

    def test_plus_strand_counts_downstream_a_content(self):
        sequence = list("C" * 120)
        sequence[self.position:self.position + self.window_size] = "A" * self.window_size
        filter_obj = self._filter_for_sequence(sequence)

        a_content = filter_obj._calculate_a_content(
            self.chromosome, self.position, "+"
        )

        self.assertEqual(a_content, 1.0)

    def test_minus_strand_uses_downstream_sequence(self):
        sequence = list("C" * 120)
        sequence[self.position:self.position + self.window_size] = "T" * self.window_size
        filter_obj = self._filter_for_sequence(sequence)

        a_content = filter_obj._calculate_a_content(
            self.chromosome, self.position, "-"
        )

        self.assertEqual(a_content, 0.0)

    def test_minus_strand_counts_downstream_t_content(self):
        sequence = list("C" * 120)
        sequence[self.position - self.window_size - 1:self.position - 1] = "T" * self.window_size
        filter_obj = self._filter_for_sequence(sequence)

        a_content = filter_obj._calculate_a_content(
            self.chromosome, self.position, "-"
        )

        self.assertEqual(a_content, 1.0)

    def test_plus_strand_excludes_site_and_21st_downstream_base(self):
        sequence = list("C" * 120)
        sequence[self.position - 1] = "A"
        sequence[self.position + self.window_size] = "A"
        filter_obj = self._filter_for_sequence(sequence)

        a_content = filter_obj._calculate_a_content(
            self.chromosome, self.position, "+"
        )

        self.assertEqual(a_content, 0.0)

    def test_minus_strand_excludes_site_and_21st_downstream_base(self):
        sequence = list("C" * 120)
        sequence[self.position - 1] = "T"
        sequence[self.position - self.window_size - 2] = "T"
        filter_obj = self._filter_for_sequence(sequence)

        a_content = filter_obj._calculate_a_content(
            self.chromosome, self.position, "-"
        )

        self.assertEqual(a_content, 0.0)

    def test_secondary_site_with_half_a_content_is_kept(self):
        sequence = list("C" * 120)
        sequence[self.position:self.position + 10] = "A" * 10
        filter_obj = self._filter_for_sequence(sequence)
        filter_obj.a_content_threshold = 0.5

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

        self.assertEqual(filtered.site, [31, self.position])
        self.assertEqual(removed.site, [])

    def test_secondary_site_above_half_a_content_is_removed(self):
        sequence = list("C" * 120)
        sequence[self.position:self.position + 11] = "A" * 11
        filter_obj = self._filter_for_sequence(sequence)
        filter_obj.a_content_threshold = 0.5

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

        self.assertEqual(filtered.site, [31])
        self.assertEqual(removed.site, [self.position])

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

    def test_minus_strand_removed_site_records_t_as_a_content(self):
        sequence = list("C" * 120)
        sequence[
            self.position - self.window_size - 1:self.position - 10
        ] = "T" * 11
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


if __name__ == "__main__":
    unittest.main()
