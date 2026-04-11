import csv
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).parent.parent))
import dnds

FILES_DIR = Path(__file__).parent / "files"
ALN_FILE = FILES_DIR / "aln.fas"
CDS_FILE = FILES_DIR / "cds.fa"

ALN_IDS = [
    "EAL73549", "EAL71978", "EAL71980", "EAL71985", "EAL72046",
    "EAL71292", "EDR41074", "EAL67447", "EAL67458", "EAL64800",
]


class TestParseArgs(unittest.TestCase):
    def test_required_args_and_defaults(self):
        with patch("sys.argv", ["dnds.py", "-m", "aln.fas", "-f", "cds.fa", "-o", "out.tsv"]):
            args = dnds.parse_args()
        self.assertEqual(args.msa, "aln.fas")
        self.assertEqual(args.fasta, "cds.fa")
        self.assertEqual(args.output, "out.tsv")
        self.assertEqual(args.msa_format, "fasta")
        self.assertEqual(args.yn00_command, "yn00")
        self.assertFalse(args.yn00_verbose)

    def test_custom_options(self):
        with patch("sys.argv", [
            "dnds.py", "-m", "a.fas", "-f", "b.fa", "-o", "c.tsv",
            "--msa-format", "clustal",
            "--yn00-command", "/usr/bin/yn00",
            "--yn00-verbose",
        ]):
            args = dnds.parse_args()
        self.assertEqual(args.msa_format, "clustal")
        self.assertEqual(args.yn00_command, "/usr/bin/yn00")
        self.assertTrue(args.yn00_verbose)


class TestReadProteinAlignment(unittest.TestCase):
    def test_reads_aln_file_ids(self):
        _, ids = dnds.read_protein_alignment(ALN_FILE, "fasta")
        self.assertEqual(ids, ALN_IDS)

    def test_returns_alignment_object(self):
        alignment, ids = dnds.read_protein_alignment(ALN_FILE, "fasta")
        self.assertIsNotNone(alignment)
        self.assertEqual(len(ids), len(ALN_IDS))

    def test_duplicate_ids_raise_value_error(self):
        content = ">SEQ_A\nMTK\n>SEQ_A\nMTK\n"
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fas", delete=False) as fh:
            fh.write(content)
            tmp = Path(fh.name)
        try:
            with self.assertRaises(ValueError) as ctx:
                dnds.read_protein_alignment(tmp, "fasta")
            self.assertIn("SEQ_A", str(ctx.exception))
        finally:
            tmp.unlink()


class TestBuildNucleotideIndex(unittest.TestCase):
    def test_builds_index_containing_all_aln_ids(self):
        index = dnds.build_nucleotide_index(CDS_FILE)
        for seq_id in ALN_IDS:
            self.assertIn(seq_id, index)

    def test_values_are_seqrecords(self):
        index = dnds.build_nucleotide_index(CDS_FILE)
        for record in index.values():
            self.assertIsInstance(record, SeqRecord)

    def test_duplicate_ids_raise_value_error(self):
        content = ">SEQ_A\nATGCCC\n>SEQ_A\nATGCCC\n"
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fh:
            fh.write(content)
            tmp = Path(fh.name)
        try:
            with self.assertRaises(ValueError) as ctx:
                dnds.build_nucleotide_index(tmp)
            self.assertIn("SEQ_A", str(ctx.exception))
        finally:
            tmp.unlink()


class TestFindMatchingNucleotides(unittest.TestCase):
    def setUp(self):
        self.index = dnds.build_nucleotide_index(CDS_FILE)

    def test_all_present_returns_no_missing(self):
        matched, missing = dnds.find_matching_nucleotides(ALN_IDS, self.index)
        self.assertEqual(len(matched), len(ALN_IDS))
        self.assertEqual(missing, [])

    def test_missing_id_reported(self):
        ids = ["EAL73549", "MISSING_ID"]
        matched, missing = dnds.find_matching_nucleotides(ids, self.index)
        self.assertEqual(len(matched), 1)
        self.assertEqual(missing, ["MISSING_ID"])

    def test_all_missing(self):
        matched, missing = dnds.find_matching_nucleotides(["X", "Y"], self.index)
        self.assertEqual(matched, [])
        self.assertEqual(missing, ["X", "Y"])

    def test_empty_ids(self):
        matched, missing = dnds.find_matching_nucleotides([], self.index)
        self.assertEqual(matched, [])
        self.assertEqual(missing, [])


class TestMakeShortIdMap(unittest.TestCase):
    def test_maps_to_s_prefixed_padded_ids(self):
        result = dnds.make_short_id_map(["alpha", "beta", "gamma"])
        self.assertEqual(result["alpha"], "S00001")
        self.assertEqual(result["beta"], "S00002")
        self.assertEqual(result["gamma"], "S00003")

    def test_length_matches_input(self):
        result = dnds.make_short_id_map(["a", "b", "c", "d"])
        self.assertEqual(len(result), 4)

    def test_empty_input(self):
        self.assertEqual(dnds.make_short_id_map([]), {})

    def test_index_exceeding_padding_width(self):
        ids = [f"SEQ_{i}" for i in range(100001)]
        result = dnds.make_short_id_map(ids)
        self.assertEqual(result["SEQ_0"], "S00001")
        self.assertEqual(result["SEQ_100000"], "S100001")


class TestRestoreOriginalIds(unittest.TestCase):
    def setUp(self):
        self.reverse_map = {"S00001": "alpha", "S00002": "beta"}

    def test_flat_dict_keys_restored(self):
        result = dnds._restore_original_ids({"S00001": 1, "S00002": 2}, self.reverse_map)
        self.assertEqual(result, {"alpha": 1, "beta": 2})

    def test_nested_dict_keys_restored(self):
        value = {"S00001": {"S00002": {"omega": 0.5}}}
        result = dnds._restore_original_ids(value, self.reverse_map)
        self.assertIn("alpha", result)
        self.assertIn("beta", result["alpha"])
        self.assertEqual(result["alpha"]["beta"]["omega"], 0.5)

    def test_list_items_not_substituted(self):
        result = dnds._restore_original_ids(["S00001", "S00002"], self.reverse_map)
        self.assertEqual(result, ["S00001", "S00002"])

    def test_scalar_returned_unchanged(self):
        self.assertEqual(dnds._restore_original_ids(42, self.reverse_map), 42)

    def test_unmapped_key_preserved(self):
        result = dnds._restore_original_ids({"S00001": 1, "UNKNOWN": 2}, self.reverse_map)
        self.assertIn("alpha", result)
        self.assertIn("UNKNOWN", result)


class TestResolveYn00Command(unittest.TestCase):
    def test_resolves_command_from_path(self):
        result = dnds.resolve_yn00_command("python3")
        self.assertTrue(Path(result).is_file())

    def test_resolves_absolute_path(self):
        import shutil
        python_path = shutil.which("python3")
        result = dnds.resolve_yn00_command(python_path)
        self.assertEqual(result, python_path)

    def test_raises_file_not_found_for_missing_command(self):
        with self.assertRaises(FileNotFoundError):
            dnds.resolve_yn00_command("no_such_command_xyzzy_12345")


class TestRunYn00(unittest.TestCase):
    def _make_records(self, names, seq="ATGCCCATG"):
        return [SeqRecord(Seq(seq), id=name, name=name, description="") for name in names]

    @patch("dnds.yn00.Yn00")
    def test_restores_original_ids_in_result(self, MockYn00):
        names = ["SEQ_A", "SEQ_B"]
        records = self._make_records(names)
        mock_runner = MagicMock()
        MockYn00.return_value = mock_runner
        mock_runner.run.return_value = {
            "S00001": {"S00002": {"NG86": {"omega": 0.5}}},
            "S00002": {"S00001": {"NG86": {"omega": 0.5}}},
        }

        result = dnds.run_yn00(records, names, "yn00", False)

        self.assertIn("SEQ_A", result)
        self.assertIn("SEQ_B", result["SEQ_A"])

    @patch("dnds.yn00.Yn00")
    def test_raises_runtime_error_when_result_is_none(self, MockYn00):
        names = ["SEQ_A", "SEQ_B"]
        records = self._make_records(names)
        mock_runner = MagicMock()
        MockYn00.return_value = mock_runner
        mock_runner.run.return_value = None

        with self.assertRaises(RuntimeError):
            dnds.run_yn00(records, names, "yn00", False)

    @patch("dnds.yn00.Yn00")
    def test_raises_runtime_error_when_runner_raises(self, MockYn00):
        names = ["SEQ_A", "SEQ_B"]
        records = self._make_records(names)
        mock_runner = MagicMock()
        MockYn00.return_value = mock_runner
        mock_runner.run.side_effect = RuntimeError("yn00 failed")

        with self.assertRaises(RuntimeError):
            dnds.run_yn00(records, names, "yn00", False)

    @patch("dnds.yn00.Yn00")
    def test_index_error_triggers_pairwise_fallback(self, MockYn00):
        """When the full-alignment parse raises IndexError, pairwise results are returned."""
        names = ["SEQ_A", "SEQ_B"]
        records = self._make_records(names)
        mock_runner = MagicMock()
        MockYn00.return_value = mock_runner
        pairwise_result = {
            "S00001": {"S00002": {"NG86": {"omega": 0.5}}},
            "S00002": {"S00001": {"NG86": {"omega": 0.5}}},
        }
        # First call (full alignment) raises IndexError; second (pairwise pair) returns result.
        mock_runner.run.side_effect = [IndexError("parse error"), pairwise_result]

        result = dnds.run_yn00(records, names, "yn00", False)

        self.assertIn("SEQ_A", result)
        self.assertIn("SEQ_B", result["SEQ_A"])
        self.assertEqual(result["SEQ_A"]["SEQ_B"]["NG86"]["omega"], 0.5)

    @patch("dnds.yn00.Yn00")
    def test_pairwise_fallback_failed_pair_yields_empty_dict(self, MockYn00):
        """When a pairwise run fails, that pair is represented as an empty dict (→ NA)."""
        names = ["SEQ_A", "SEQ_B"]
        records = self._make_records(names)
        mock_runner = MagicMock()
        MockYn00.return_value = mock_runner
        # Full alignment fails with IndexError; the single pairwise run also fails.
        mock_runner.run.side_effect = [IndexError("full parse"), Exception("pair failed")]

        result = dnds.run_yn00(records, names, "yn00", False)

        self.assertEqual(result.get("SEQ_A", {}).get("SEQ_B"), {})
        self.assertEqual(result.get("SEQ_B", {}).get("SEQ_A"), {})


class TestFmt(unittest.TestCase):
    def test_none_returns_na(self):
        self.assertEqual(dnds._fmt(None), "NA")

    def test_nan_returns_na(self):
        self.assertEqual(dnds._fmt(float("nan")), "NA")

    def test_inf_returns_na(self):
        self.assertEqual(dnds._fmt(float("inf")), "NA")

    def test_neg_inf_returns_na(self):
        self.assertEqual(dnds._fmt(float("-inf")), "NA")

    def test_regular_float_returns_string(self):
        self.assertEqual(dnds._fmt(1.23), "1.23")

    def test_zero_returns_string(self):
        self.assertEqual(dnds._fmt(0), "0")

    def test_string_returned_unchanged(self):
        self.assertEqual(dnds._fmt("hello"), "hello")


class TestCalculatePctIdentity(unittest.TestCase):
    def test_identical_sequences(self):
        pct_a, pct_b = dnds.calculate_pct_identity("ATGATG", "ATGATG")
        self.assertAlmostEqual(pct_a, 1.0)
        self.assertAlmostEqual(pct_b, 1.0)

    def test_completely_different_sequences(self):
        pct_a, pct_b = dnds.calculate_pct_identity("AAAA", "TTTT")
        self.assertAlmostEqual(pct_a, 0.0)
        self.assertAlmostEqual(pct_b, 0.0)

    def test_gap_in_b_reduces_pct_a_but_not_pct_b(self):
        # A = "ATGATG", B = "ATG---"
        # matches = 3 (positions 0-2), non_gap_A = 6, non_gap_B = 3
        # pct_A_in_B = 3/6 = 0.5, pct_B_in_A = 3/3 = 1.0
        pct_a, pct_b = dnds.calculate_pct_identity("ATGATG", "ATG---")
        self.assertAlmostEqual(pct_a, 0.5)
        self.assertAlmostEqual(pct_b, 1.0)

    def test_all_gaps_in_a_returns_zero_for_both(self):
        # A is all gaps: no non-gap positions so 0 matches
        pct_a, pct_b = dnds.calculate_pct_identity("---", "ATG")
        self.assertAlmostEqual(pct_a, 0.0)
        self.assertAlmostEqual(pct_b, 0.0)

    def test_all_gaps_in_b_returns_zero_for_b(self):
        pct_a, pct_b = dnds.calculate_pct_identity("ATG", "---")
        self.assertAlmostEqual(pct_a, 0.0)
        self.assertAlmostEqual(pct_b, 0.0)

    def test_partial_identity(self):
        # matches at positions 0,1,2; mismatch at position 3
        pct_a, pct_b = dnds.calculate_pct_identity("AATT", "AATG")
        self.assertAlmostEqual(pct_a, 0.75)
        self.assertAlmostEqual(pct_b, 0.75)

    def test_custom_gap_char(self):
        # A = "ATG..", B = "ATGTT"; gap_char="."
        # matches = 3 (positions 0-2), non_gap_A = 3, non_gap_B = 5
        pct_a, pct_b = dnds.calculate_pct_identity("ATG..", "ATGTT", gap_char=".")
        self.assertAlmostEqual(pct_a, 1.0)
        self.assertAlmostEqual(pct_b, 0.6)


class TestComputePairwisePctIdentities(unittest.TestCase):
    def _make_alignment(self, seqs: dict):
        records = [SeqRecord(Seq(seq), id=name, name=name) for name, seq in seqs.items()]
        return MultipleSeqAlignment(records)

    def test_both_orderings_stored(self):
        aln = self._make_alignment({"A": "ATGATG", "B": "ATGATG"})
        pct_ids = dnds.compute_pairwise_pct_identities(aln)
        self.assertIn(("A", "B"), pct_ids)
        self.assertIn(("B", "A"), pct_ids)

    def test_identical_sequences_give_100_pct(self):
        aln = self._make_alignment({"A": "ATGATG", "B": "ATGATG"})
        pct_ids = dnds.compute_pairwise_pct_identities(aln)
        self.assertAlmostEqual(pct_ids[("A", "B")], 1.0)
        self.assertAlmostEqual(pct_ids[("B", "A")], 1.0)

    def test_asymmetric_gaps_give_asymmetric_pct(self):
        # A = "ATGATG", B = "ATG---"
        # pct_A_in_B = 3/6 = 0.5, pct_B_in_A = 3/3 = 1.0
        aln = self._make_alignment({"A": "ATGATG", "B": "ATG---"})
        pct_ids = dnds.compute_pairwise_pct_identities(aln)
        self.assertAlmostEqual(pct_ids[("A", "B")], 0.5)
        self.assertAlmostEqual(pct_ids[("B", "A")], 1.0)

    def test_three_sequences_has_six_entries(self):
        aln = self._make_alignment({"A": "ATGATG", "B": "ATGATG", "C": "ATGATG"})
        pct_ids = dnds.compute_pairwise_pct_identities(aln)
        # 3 pairs × 2 orderings = 6 entries
        self.assertEqual(len(pct_ids), 6)

    def test_single_sequence_returns_empty(self):
        aln = self._make_alignment({"A": "ATGATG"})
        pct_ids = dnds.compute_pairwise_pct_identities(aln)
        self.assertEqual(pct_ids, {})


class TestBuildTsvHeader(unittest.TestCase):
    def test_first_four_columns(self):
        header = dnds.build_tsv_header()
        self.assertEqual(header[0], "seq1")
        self.assertEqual(header[1], "seq2")
        self.assertEqual(header[2], "pct_id_seq1_in_seq2")
        self.assertEqual(header[3], "pct_id_seq2_in_seq1")

    def test_pct_id_columns_present(self):
        header = dnds.build_tsv_header()
        self.assertIn("pct_id_seq1_in_seq2", header)
        self.assertIn("pct_id_seq2_in_seq1", header)

    def test_all_method_fields_present(self):
        header = dnds.build_tsv_header()
        for method in dnds.METHODS:
            for field in dnds.METHOD_FIELDS[method]:
                self.assertIn(f"{method}_{field}", header)

    def test_no_duplicate_columns(self):
        header = dnds.build_tsv_header()
        self.assertEqual(len(header), len(set(header)))


class TestFlattenResults(unittest.TestCase):
    def _empty_methods(self):
        return {m: {} for m in dnds.METHODS}

    def test_basic_pair(self):
        results = {"A": {"B": self._empty_methods()}}
        rows = dnds.flatten_results(results)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["seq1"], "A")
        self.assertEqual(rows[0]["seq2"], "B")

    def test_deduplicates_symmetric_pairs(self):
        methods = self._empty_methods()
        results = {"A": {"B": methods}, "B": {"A": methods}}
        rows = dnds.flatten_results(results)
        self.assertEqual(len(rows), 1)

    def test_missing_method_data_yields_na(self):
        results = {"A": {"B": {}}}
        rows = dnds.flatten_results(results)
        self.assertEqual(len(rows), 1)
        for method in dnds.METHODS:
            for field in dnds.METHOD_FIELDS[method]:
                self.assertEqual(rows[0][f"{method}_{field}"], "NA")

    def test_empty_results(self):
        self.assertEqual(dnds.flatten_results({}), [])

    def test_skips_non_dict_partners(self):
        rows = dnds.flatten_results({"A": "not_a_dict"})
        self.assertEqual(rows, [])

    def test_populates_ng86_fields(self):
        methods = {
            "NG86": {"omega": 0.3, "dN": 0.1, "dS": 0.5},
            **{m: {} for m in dnds.METHODS if m != "NG86"},
        }
        rows = dnds.flatten_results({"A": {"B": methods}})
        self.assertEqual(rows[0]["NG86_omega"], "0.3")
        self.assertEqual(rows[0]["NG86_dN"], "0.1")
        self.assertEqual(rows[0]["NG86_dS"], "0.5")

    def test_pct_id_values_included_in_output(self):
        pct_ids = {("A", "B"): 0.8, ("B", "A"): 0.9}
        rows = dnds.flatten_results({"A": {"B": self._empty_methods()}}, pct_ids)
        self.assertEqual(rows[0]["pct_id_seq1_in_seq2"], "80.0000")
        self.assertEqual(rows[0]["pct_id_seq2_in_seq1"], "90.0000")

    def test_pct_id_is_na_when_not_provided(self):
        rows = dnds.flatten_results({"A": {"B": self._empty_methods()}}, None)
        self.assertEqual(rows[0]["pct_id_seq1_in_seq2"], "NA")
        self.assertEqual(rows[0]["pct_id_seq2_in_seq1"], "NA")

    def test_zero_pct_id_seq1_in_seq2_suppresses_yn00(self):
        methods = {"NG86": {"omega": 0.3, "dN": 0.1, "dS": 0.5}, **{m: {} for m in dnds.METHODS if m != "NG86"}}
        pct_ids = {("A", "B"): 0.0, ("B", "A"): 0.8}
        rows = dnds.flatten_results({"A": {"B": methods}}, pct_ids)
        self.assertEqual(rows[0]["NG86_omega"], "NA")
        self.assertEqual(rows[0]["NG86_dN"], "NA")
        # pct_id column should still be present and set
        self.assertEqual(rows[0]["pct_id_seq1_in_seq2"], "0.0000")

    def test_zero_pct_id_seq2_in_seq1_suppresses_yn00(self):
        methods = {"NG86": {"omega": 0.3, "dN": 0.1, "dS": 0.5}, **{m: {} for m in dnds.METHODS if m != "NG86"}}
        pct_ids = {("A", "B"): 0.8, ("B", "A"): 0.0}
        rows = dnds.flatten_results({"A": {"B": methods}}, pct_ids)
        self.assertEqual(rows[0]["NG86_omega"], "NA")
        self.assertEqual(rows[0]["pct_id_seq2_in_seq1"], "0.0000")

    def test_nonzero_pct_id_preserves_yn00_results(self):
        methods = {"NG86": {"omega": 0.3, "dN": 0.1, "dS": 0.5}, **{m: {} for m in dnds.METHODS if m != "NG86"}}
        pct_ids = {("A", "B"): 0.8, ("B", "A"): 0.9}
        rows = dnds.flatten_results({"A": {"B": methods}}, pct_ids)
        self.assertEqual(rows[0]["NG86_omega"], "0.3")

    def test_none_pct_ids_does_not_suppress_yn00(self):
        methods = {"NG86": {"omega": 0.3}, **{m: {} for m in dnds.METHODS if m != "NG86"}}
        rows = dnds.flatten_results({"A": {"B": methods}}, None)
        self.assertEqual(rows[0]["NG86_omega"], "0.3")


class TestWriteTsv(unittest.TestCase):
    def _make_row(self):
        row = {"seq1": "A", "seq2": "B", "pct_id_seq1_in_seq2": "100.0000", "pct_id_seq2_in_seq1": "100.0000"}
        for method in dnds.METHODS:
            for field in dnds.METHOD_FIELDS[method]:
                row[f"{method}_{field}"] = "NA"
        return row

    def test_writes_header_and_data_row(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "results.tsv"
            dnds.write_tsv([self._make_row()], output)
            self.assertTrue(output.exists())
            with output.open() as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["seq1"], "A")
            self.assertEqual(rows[0]["seq2"], "B")

    def test_creates_parent_directories(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "nested" / "dir" / "results.tsv"
            dnds.write_tsv([], output)
            self.assertTrue(output.exists())

    def test_empty_rows_writes_only_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "results.tsv"
            dnds.write_tsv([], output)
            lines = output.read_text().strip().splitlines()
            self.assertEqual(len(lines), 1)
            self.assertTrue(lines[0].startswith("seq1\tseq2"))


if __name__ == "__main__":
    unittest.main()

