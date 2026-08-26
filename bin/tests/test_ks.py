"""Unit tests for bin/ks.py.

The PAML fixtures below are copied verbatim from real codeml/yn00 runs, since
the output format is the thing most likely to break silently on a version bump.
"""

import itertools
import math
import random
import sys
import tempfile
import unittest
from io import StringIO
from pathlib import Path

from Bio import Phylo

sys.path.insert(0, str(Path(__file__).parent.parent))
import ks


# A real M0 (`runmode = 0`) codeml result for a six-sequence alignment.
CODEML_TREE_OUTPUT = """
lnL(ntime:  9  np: 11):  -2577.851885      +0.000000

tree length =  0.499467

Detailed output identifying parameters

kappa (ts/tv) =  1.42208

omega (dN/dS) =  0.04619

dN & dS for each branch

 branch          t       N       S   dN/dS      dN      dS  N*dN  S*dS

   7..8      0.021   880.2   319.8  0.0462  0.0011  0.0228   0.9   7.3
   8..9      0.049   880.2   319.8  0.0462  0.0025  0.0546   2.2  17.4
   9..1      0.041   880.2   319.8  0.0462  0.0021  0.0451   1.8  14.4
   9..2      0.030   880.2   319.8  0.0462  0.0015  0.0329   1.3  10.5
   8..10     0.046   880.2   319.8  0.0462  0.0024  0.0515   2.1  16.5
  10..3      0.033   880.2   319.8  0.0462  0.0017  0.0364   1.5  11.6
  10..4      0.032   880.2   319.8  0.0462  0.0016  0.0351   1.4  11.2
   7..5      0.119   880.2   319.8  0.0462  0.0061  0.1322   5.4  42.3
   7..6      0.130   880.2   319.8  0.0462  0.0066  0.1438   5.8  46.0

tree length for dN:       0.0256
tree length for dS:       0.5543
"""

CODEML_PAIRWISE_OUTPUT = """
2 (S00002) ... 1 (S00001)
lnL = -1234.567890
   0.06930  1.92340  0.14850

t= 0.0693  S=   311.4  N=   888.6  dN/dS=  0.1485  dN = 0.0093  dS = 0.0625
"""

# yn00 reports several methods; only the Yang & Nielsen block should be read.
YN00_OUTPUT = """
(A) Nei-Gojobori (1986) method

Nei & Gojobori 1986. dN/dS (dN, dS)
S00001
S00002  0.1234 (0.0111 0.0900)

(B) Yang & Nielsen (2000) method

Yang Z, Nielsen R (2000) Estimating synonymous and nonsynonymous substitution
rates under realistic evolutionary models. Mol. Biol. Evol. 17:32-43

(equal weighting of pathways)

seq. seq.     S       N        t   kappa   omega     dN +- SE    dS +- SE

   2    1   370.8   829.2   0.0694  4.6000  0.2265 0.0112 +- 0.0037  0.0497 +- 0.0119
"""


class TestBranchTableParsing(unittest.TestCase):
    def test_parses_every_branch(self):
        edges = ks.parse_branch_table(CODEML_TREE_OUTPUT)
        self.assertEqual(len(edges), 9)
        self.assertIn((9, 1, 0.0021, 0.0451), edges)
        self.assertIn((7, 6, 0.0066, 0.1438), edges)

    def test_returns_nothing_without_the_section(self):
        self.assertEqual(ks.parse_branch_table("no table here"), [])

    def test_branch_dS_sums_to_reported_tree_length(self):
        # Confirms the dS column is the one being read. Only 3 places: codeml
        # prints per-branch dS to 4 decimals, so summing nine of them drifts
        # from the reported total by rounding alone.
        edges = ks.parse_branch_table(CODEML_TREE_OUTPUT)
        self.assertAlmostEqual(sum(e[3] for e in edges), 0.5543, places=3)


class TestTreePathDistances(unittest.TestCase):
    def setUp(self):
        self.edges = ks.parse_branch_table(CODEML_TREE_OUTPUT)
        self.nodes = {f"S{i:05d}": i for i in range(1, 7)}

    def test_sibling_path_sums_two_branches(self):
        # S00001 and S00002 are joined through node 9: 0.0451 + 0.0329
        distances = ks.tree_path_distances(self.edges, self.nodes, ["S00001"])
        d_s, d_n = distances[("S00001", "S00002")]
        self.assertAlmostEqual(d_s, 0.0780, places=4)
        self.assertAlmostEqual(d_n, 0.0036, places=4)

    def test_deeper_path_traverses_internal_branches(self):
        # S00001 -> 9 -> 8 -> 10 -> S00003
        distances = ks.tree_path_distances(self.edges, self.nodes, ["S00001"])
        d_s, _ = distances[("S00001", "S00003")]
        self.assertAlmostEqual(d_s, 0.0451 + 0.0546 + 0.0515 + 0.0364, places=4)

    def test_only_requested_tips_are_computed(self):
        distances = ks.tree_path_distances(self.edges, self.nodes, ["S00001"])
        self.assertTrue(all(pair[0] == "S00001" for pair in distances))

    def test_distances_are_symmetric(self):
        distances = ks.tree_path_distances(self.edges, self.nodes, ["S00001", "S00003"])
        self.assertAlmostEqual(
            distances[("S00001", "S00003")][0], distances[("S00003", "S00001")][0], places=6
        )


class TestResultLineParsing(unittest.TestCase):
    def test_codeml_pairwise_line(self):
        match = ks._RE_PAIRWISE.search(CODEML_PAIRWISE_OUTPUT)
        self.assertIsNotNone(match)
        values = {k: float(v) for k, v in match.groupdict().items()}
        self.assertAlmostEqual(values["dS"], 0.0625)
        self.assertAlmostEqual(values["dN"], 0.0093)
        self.assertAlmostEqual(values["omega"], 0.1485)
        self.assertAlmostEqual(values["S"], 311.4)

    def test_yn00_line_read_from_the_right_section(self):
        marker = YN00_OUTPUT.find(ks._YN00_SECTION)
        self.assertNotEqual(marker, -1)
        match = ks._RE_YN00.search(YN00_OUTPUT[marker:])
        self.assertIsNotNone(match)
        self.assertAlmostEqual(float(match.group("dS")), 0.0497)
        self.assertAlmostEqual(float(match.group("dN")), 0.0112)
        self.assertAlmostEqual(float(match.group("omega")), 0.2265)

    def test_model_level_statistics(self):
        self.assertAlmostEqual(
            float(ks._RE_OMEGA.search(CODEML_TREE_OUTPUT).group(1)), 0.04619
        )
        self.assertAlmostEqual(
            float(ks._RE_KAPPA.search(CODEML_TREE_OUTPUT).group(1)), 1.42208
        )
        self.assertAlmostEqual(
            float(ks._RE_LNL.search(CODEML_TREE_OUTPUT).group(1)), -2577.851885
        )


class TestPairwiseGapStripping(unittest.TestCase):
    def test_drops_codons_gapped_in_either_sequence(self):
        a = "AAA" "---" "CCC" "GGG"
        b = "AAA" "TTT" "---" "GGG"
        out_a, out_b = ks.strip_pair_gaps(a, b)
        self.assertEqual(out_a, "AAAGGG")
        self.assertEqual(out_b, "AAAGGG")

    def test_keeps_everything_when_ungapped(self):
        a = b = "AAACCCGGG"
        self.assertEqual(ks.strip_pair_gaps(a, b), (a, b))

    def test_result_length_is_a_whole_number_of_codons(self):
        a = "AAA" "-CC" "GGG"
        b = "AAA" "TTT" "GGG"
        out_a, _ = ks.strip_pair_gaps(a, b)
        self.assertEqual(len(out_a) % 3, 0)

    def test_no_shared_columns_gives_empty(self):
        self.assertEqual(ks.strip_pair_gaps("---AAA", "AAA---"), ("", ""))


class TestPctIdentity(unittest.TestCase):
    def test_identical_sequences(self):
        self.assertEqual(ks.pct_identity("AAACCC", "AAACCC"), (1.0, 1.0))

    def test_gaps_count_as_mismatches_asymmetrically(self):
        # 6 non-gap in a, 3 in b, 3 matching positions
        a_in_b, b_in_a = ks.pct_identity("AAACCC", "AAA---")
        self.assertAlmostEqual(a_in_b, 0.5)
        self.assertAlmostEqual(b_in_a, 1.0)

    def test_all_gap_sequence_is_zero(self):
        self.assertEqual(ks.pct_identity("------", "AAACCC"), (0.0, 0.0))


class TestShortIds(unittest.TestCase):
    def test_ids_are_paml_safe_and_ordered(self):
        mapping = ks.make_short_id_map(["long|gene|name|1", "b", "c"])
        self.assertEqual(mapping["long|gene|name|1"], "S00001")
        self.assertEqual(mapping["c"], "S00003")
        self.assertTrue(all(len(v) <= 10 for v in mapping.values()))


class TestPhylipWriting(unittest.TestCase):
    def test_strict_sequential_layout(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "a.phy"
            ks.write_phylip([("S00001", "AAACCC"), ("S00002", "AAAGGG")], path)
            lines = path.read_text().splitlines()
        self.assertEqual(lines[0].strip(), "2 6")
        self.assertEqual(lines[1], "S00001")
        self.assertEqual(lines[2], "AAACCC")

    def test_rejects_ragged_alignment(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ValueError):
                ks.write_phylip([("a", "AAA"), ("b", "AA")], Path(tmp) / "a.phy")

    def test_rejects_empty_alignment(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ValueError):
                ks.write_phylip([], Path(tmp) / "a.phy")


class TestPamlTimeout(unittest.TestCase):
    def test_timeout_raises_a_clear_error(self):
        # Saturated sequences make CODEML crawl toward an unbounded dS, so a
        # single orthogroup can otherwise eat a whole task's time budget.
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(RuntimeError) as caught:
                ks._run_paml(["sleep", "30"], Path(tmp), timeout=1)
        self.assertIn("1s limit", str(caught.exception))

    def test_zero_means_no_limit(self):
        with tempfile.TemporaryDirectory() as tmp:
            proc = ks._run_paml(["true"], Path(tmp), timeout=0)
        self.assertEqual(proc.returncode, 0)

    def test_process_within_the_limit_is_untouched(self):
        with tempfile.TemporaryDirectory() as tmp:
            proc = ks._run_paml(["echo", "hi"], Path(tmp), timeout=30)
        self.assertEqual(proc.stdout.strip(), "hi")


class TestColumns(unittest.TestCase):
    def test_species_is_reported_per_pair(self):
        # Pairs are formed within a species, so each row names the species it
        # belongs to rather than leaving the caller to infer it.
        self.assertIn("species", ks.COLUMNS)
        self.assertEqual(ks.COLUMNS[:4], ["orthogroup", "species", "gene_a", "gene_b"])


class TestTsvWriting(unittest.TestCase):
    def test_uses_unix_line_endings(self):
        # csv.DictWriter defaults to \r\n, which leaves a stray \r on the last
        # field of every row and breaks awk/cut and any naive column split.
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "out.tsv"
            ks.write_tsv([{"orthogroup": "OG1", "gene_a": "a", "gene_b": "b"}], path)
            raw = path.read_bytes()
        self.assertNotIn(b"\r", raw)

    def test_writes_the_full_column_set_with_NA_for_missing(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "out.tsv"
            ks.write_tsv([{"orthogroup": "OG1", "gene_a": "a", "gene_b": "b"}], path)
            lines = path.read_text().splitlines()
        self.assertEqual(lines[0].split("\t"), ks.COLUMNS)
        row = dict(zip(ks.COLUMNS, lines[1].split("\t")))
        self.assertEqual(row["orthogroup"], "OG1")
        self.assertEqual(row["tree_dS"], "NA")


class TestTreePreparation(unittest.TestCase):
    def _tree_file(self, tmp, newick):
        path = Path(tmp) / "t.nwk"
        path.write_text(newick)
        return path

    def test_prunes_relabels_and_unroots(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self._tree_file(tmp, "(((a,b),(c,d)),e);")
            id_map = {n: f"S{i:05d}" for i, n in enumerate("abce", start=1)}
            tree = ks.load_and_prepare_tree(path, {"a", "b", "c", "e"}, id_map)
            names = sorted(t.name for t in tree.get_terminals())
            self.assertEqual(names, ["S00001", "S00002", "S00003", "S00004"])
            # 'd' is gone, and the root is no longer bifurcating
            self.assertGreaterEqual(len(tree.root.clades), 3)

    def test_rejects_tree_missing_alignment_sequences(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self._tree_file(tmp, "((a,b),c);")
            with self.assertRaises(ValueError):
                ks.load_and_prepare_tree(path, {"a", "b", "z"}, {})

    def test_rejects_tree_too_small_to_help(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self._tree_file(tmp, "((a,b),c);")
            with self.assertRaises(ValueError):
                ks.load_and_prepare_tree(path, {"a", "b"}, {"a": "S00001", "b": "S00002"})

    def test_written_tree_has_no_branch_lengths(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self._tree_file(tmp, "(((a:0.1,b:0.2):0.3,c:0.4),d:0.5);")
            id_map = {n: f"S{i:05d}" for i, n in enumerate("abcd", start=1)}
            tree = ks.load_and_prepare_tree(path, set("abcd"), id_map)
            out = Path(tmp) / "out.nwk"
            ks.write_tree(tree, out)
            self.assertNotIn(":", out.read_text())


class TestPolytomyResolution(unittest.TestCase):
    """CODEML's MAXNSONS is 3; OrthoFinder emits wider nodes freely.

    A single polytomy anywhere in the tree made codeml abort with `too many
    daughter nodes`, costing the whole orthogroup its tree-based dS -- the
    primary estimate -- and doing it silently, because run_codeml discarded the
    console log that carried the reason.
    """

    def _tree_file(self, tmp, newick):
        path = Path(tmp) / "t.nwk"
        path.write_text(newick)
        return path

    def _max_degree(self, tree):
        return max(len(c.clades) for c in tree.find_clades() if c.clades)

    def test_internal_polytomy_is_split(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self._tree_file(tmp, "((a,b,c,d,e),f,g);")
            names = "abcdefg"
            id_map = {n: f"S{i:05d}" for i, n in enumerate(names, start=1)}
            tree = ks.load_and_prepare_tree(path, set(names), id_map)
            self.assertLessEqual(self._max_degree(tree), ks.MAX_PAML_DAUGHTERS)

    def test_wide_root_polytomy_is_split(self):
        with tempfile.TemporaryDirectory() as tmp:
            names = [chr(ord("a") + i) for i in range(9)]
            path = self._tree_file(tmp, "(" + ",".join(names) + ");")
            id_map = {n: f"S{i:05d}" for i, n in enumerate(names, start=1)}
            tree = ks.load_and_prepare_tree(path, set(names), id_map)
            self.assertLessEqual(self._max_degree(tree), ks.MAX_PAML_DAUGHTERS)

    def test_every_tip_survives_resolution(self):
        with tempfile.TemporaryDirectory() as tmp:
            names = [chr(ord("a") + i) for i in range(12)]
            path = self._tree_file(tmp, "((" + ",".join(names[:6]) + "),(" + ",".join(names[6:]) + "));")
            id_map = {n: f"S{i:05d}" for i, n in enumerate(names, start=1)}
            tree = ks.load_and_prepare_tree(path, set(names), id_map)
            self.assertEqual(
                sorted(t.name for t in tree.get_terminals()), sorted(id_map.values())
            )

    def test_bifurcating_tree_is_left_alone(self):
        # Trees that already worked must not be perturbed by the fix.
        tree = Phylo.read(StringIO("(((a,b),(c,d)),e);"), "newick")
        self.assertEqual(ks._resolve_polytomies(tree), 0)

    def test_very_wide_polytomy_stays_shallow(self):
        # Resolved as a ladder this would be 2000 deep, and Biopython walks
        # trees recursively, so it has to come out balanced instead.
        names = [f"t{i}" for i in range(2000)]
        tree = Phylo.read(StringIO("(" + ",".join(names) + ");"), "newick")
        ks._resolve_polytomies(tree)
        self.assertLessEqual(self._max_degree(tree), ks.MAX_PAML_DAUGHTERS)
        self.assertEqual(len(tree.get_terminals()), len(names))


class TestPamlDiagnostics(unittest.TestCase):
    """PAML writes a partial output file *and* exits non-zero on a fatal error,
    so the console log is the only place the reason appears."""

    def test_error_line_is_picked_out_of_the_console_log(self):
        run = ks.PamlRun(
            text="partial",
            stdout=(
                "CODONML in paml version 4.10.7\n"
                "error: too many daughter nodes, raise MAXNSONS\n"
            ),
            returncode=1,
        )
        self.assertIn("too many daughter nodes", run.diagnostic())

    def test_falls_back_to_the_tail_when_nothing_says_error(self):
        run = ks.PamlRun(text="", stdout="some chatter\nmore chatter\n", returncode=1)
        self.assertIn("chatter", run.diagnostic())

    def test_reports_exit_status_when_silent(self):
        run = ks.PamlRun(text="", stdout="", returncode=139)
        self.assertIn("139", run.diagnostic())


class TestCrosscheckSampling(unittest.TestCase):
    """The M0 fit covers every pair from one optimisation; the pairwise and
    YN00 cross-checks cost two processes each and pairs grow quadratically, so
    above a cap they are sampled. The sample must be reproducible, or a rerun
    would silently cross-check a different set."""

    def _pairs(self, n):
        genes = [f"g{i}" for i in range(n)]
        return [("sp", a, b) for a, b in itertools.combinations(genes, 2)]

    def test_sampling_is_seeded_on_the_orthogroup(self):
        pairs = self._pairs(60)
        first = random.Random("OG0000002").sample(pairs, 100)
        again = random.Random("OG0000002").sample(pairs, 100)
        self.assertEqual(first, again)

    def test_different_orthogroups_sample_differently(self):
        pairs = self._pairs(60)
        self.assertNotEqual(
            random.Random("OG0000002").sample(pairs, 100),
            random.Random("OG0000003").sample(pairs, 100),
        )

    def test_sample_is_a_subset_of_the_pairs(self):
        pairs = self._pairs(40)
        sample = random.Random("OG").sample(pairs, 50)
        self.assertEqual(len(sample), 50)
        self.assertTrue(set(sample) <= set(pairs))


class TestCrosscheckColumn(unittest.TestCase):
    def test_crosschecked_sits_next_to_has_tree(self):
        # Downstream needs to tell an NA caused by sampling from one caused by
        # a failure, so the flag travels with the row.
        self.assertIn("crosschecked", ks.COLUMNS)
        self.assertEqual(
            ks.COLUMNS[ks.COLUMNS.index("has_tree") + 1], "crosschecked"
        )

    def test_pair_columns_still_follow(self):
        for name in ("pair_dS", "yn00_dS", "tree_dS", "dS_tree_over_pair"):
            self.assertIn(name, ks.COLUMNS)


class TestReportableSelection(unittest.TestCase):
    def test_species_below_two_survivors_is_dropped(self):
        paralogs = {"sp1": ["a", "b", "c"], "sp2": ["d", "e"]}
        self.assertEqual(
            ks.select_reportable(paralogs, {"a", "b", "d"}), {"sp1": ["a", "b"]}
        )

    def test_empty_when_nothing_survives(self):
        self.assertEqual(ks.select_reportable({"sp1": ["a", "b"]}, {"a"}), {})


class TestCdsIndexing(unittest.TestCase):
    def _fasta(self, tmp, name, records):
        path = Path(tmp) / name
        path.write_text("".join(f">{i}\n{s}\n" for i, s in records))
        return path

    def test_indexes_across_multiple_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            a = self._fasta(tmp, "a.fa", [("g1", "ATG")])
            b = self._fasta(tmp, "b.fa", [("g2", "AAA")])
            index = ks.build_cds_index([a, b])
            self.assertEqual(sorted(index), ["g1", "g2"])

    def test_duplicate_id_across_species_is_an_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            a = self._fasta(tmp, "a.fa", [("g1", "ATG")])
            b = self._fasta(tmp, "b.fa", [("g1", "AAA")])
            with self.assertRaises(ValueError):
                ks.build_cds_index([a, b])

    def test_empty_input_is_an_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            empty = self._fasta(tmp, "a.fa", [])
            with self.assertRaises(ValueError):
                ks.build_cds_index([empty])

    def test_matching_is_exact_and_reports_misses(self):
        index = {"g1": object(), "g2": object()}
        matched, missing = ks.match_cds(["g1", "g3", "g2"], index)
        self.assertEqual(len(matched), 2)
        self.assertEqual(missing, ["g3"])

    def test_no_fuzzy_version_suffix_matching(self):
        # A versioned id must not silently match its unversioned counterpart.
        _, missing = ks.match_cds(["g1.2"], {"g1": object()})
        self.assertEqual(missing, ["g1.2"])


class TestManifest(unittest.TestCase):
    def _manifest(self, tmp, text):
        path = Path(tmp) / "m.txt"
        path.write_text(text)
        return path

    def test_groups_genes_by_species(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self._manifest(tmp, "species\tgene\nspA\tg1\nspA\tg2\nspB\tg3\n")
            self.assertEqual(ks.read_manifest(path), {"spA": ["g1", "g2"], "spB": ["g3"]})

    def test_header_is_optional(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self._manifest(tmp, "spA\tg1\nspA\tg2\n")
            self.assertEqual(ks.read_manifest(path), {"spA": ["g1", "g2"]})

    def test_blank_lines_and_padding_tolerated(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self._manifest(tmp, "\n spA \t g1 \n\n")
            self.assertEqual(ks.read_manifest(path), {"spA": ["g1"]})

    def test_single_column_is_rejected(self):
        # The old one-gene-per-line format would silently lose the species.
        with tempfile.TemporaryDirectory() as tmp:
            path = self._manifest(tmp, "g1\ng2\n")
            with self.assertRaises(ValueError):
                ks.read_manifest(path)

    def test_empty_manifest_is_an_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self._manifest(tmp, "\n\n")
            with self.assertRaises(ValueError):
                ks.read_manifest(path)


class TestNewickReconstruction(unittest.TestCase):
    def setUp(self):
        self.edges = ks.parse_branch_table(CODEML_TREE_OUTPUT)
        self.names = {i: f"gene{i}" for i in range(1, 7)}

    def test_rebuilds_topology_and_dS_branch_lengths(self):
        newick = ks.newick_from_branch_table(self.edges, self.names, "dS")
        self.assertTrue(newick.endswith(";"))
        # 9..1 has dS 0.0451 and 7..5 has dS 0.1322 in the branch table.
        self.assertIn("gene1:0.045100", newick)
        self.assertIn("gene5:0.132200", newick)
        # Root is a trifurcation: CODEML treats the tree as unrooted.
        from io import StringIO

        from Bio import Phylo

        tree = Phylo.read(StringIO(newick), "newick")
        self.assertEqual(len(tree.root.clades), 3)
        self.assertEqual(len(tree.get_terminals()), 6)

    def test_dN_uses_the_other_column(self):
        newick = ks.newick_from_branch_table(self.edges, self.names, "dN")
        self.assertIn("gene1:0.002100", newick)

    def test_every_tip_appears(self):
        newick = ks.newick_from_branch_table(self.edges, self.names, "dS")
        for name in self.names.values():
            self.assertIn(name, newick)

    def test_returns_none_without_a_single_root(self):
        self.assertIsNone(ks.newick_from_branch_table([], self.names, "dS"))

    def test_ids_needing_quoting_are_quoted(self):
        self.assertEqual(ks._quote_newick("has space"), "'has space'")
        self.assertEqual(ks._quote_newick("a,b"), "'a,b'")
        # Pipes are legal in newick labels and must be left alone.
        self.assertEqual(ks._quote_newick("gi|123|emb|X.1|"), "gi|123|emb|X.1|")

    def test_deep_ladder_does_not_hit_the_recursion_limit(self):
        # A ladder-shaped gene tree is the worst case for a recursive writer.
        depth = 2000
        edges, names = [], {}
        for i in range(1, depth):
            edges.append((i + depth, i, 0.01, 0.02))
            edges.append((i + depth, i + depth + 1, 0.01, 0.02))
            names[i] = f"t{i}"
        edges.append((2 * depth, depth, 0.01, 0.02))
        names[depth] = f"t{depth}"
        newick = ks.newick_from_branch_table(edges, names, "dS")
        self.assertIsNotNone(newick)
        self.assertIn("t1:", newick)


class TestGeneticCodeMapping(unittest.TestCase):
    def test_mycoplasma_code_maps_to_ncbi_table_4(self):
        self.assertEqual(ks.PAML_ICODE_TO_NCBI[3], 4)

    def test_universal_code_maps_to_ncbi_table_1(self):
        self.assertEqual(ks.PAML_ICODE_TO_NCBI[0], 1)

    def test_every_paml_code_resolves_to_a_real_table(self):
        from Bio.Data import CodonTable

        for ncbi_id in ks.PAML_ICODE_TO_NCBI.values():
            self.assertIn(ncbi_id, CodonTable.unambiguous_dna_by_id)


class TestFormatting(unittest.TestCase):
    def test_missing_values_render_as_NA(self):
        self.assertEqual(ks.fmt(None), "NA")
        self.assertEqual(ks.fmt(math.nan), "NA")
        self.assertEqual(ks.fmt(math.inf), "NA")

    def test_numbers_keep_significant_digits(self):
        self.assertEqual(ks.fmt(0.0625), "0.0625")
        self.assertEqual(ks.fmt("yes"), "yes")


if __name__ == "__main__":
    unittest.main()
