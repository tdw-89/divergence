"""Unit tests for bin/longest_isoform.py (both schemes and --all)."""

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
import longest_isoform as li

SCRIPT = Path(__file__).parent.parent / "longest_isoform.py"

REAL_DEFLINE = (
    "lcl|NC_133200.1_cds_NP_001001398.2_93141 [gene=fgf6a] "
    "[db_xref=GeneID:30560] [protein=fibroblast growth factor 6a precursor] "
    "[protein_id=NP_001001398.2] "
    "[location=join(17881715..17882051,17882586..17882689,17887110..17887286)] [gbkey=CDS]"
)


class TestRefseqDeflineParsing(unittest.TestCase):
    def test_real_defline(self):
        gene, protein, symbol = li.parse_refseq_cds_defline(
            REAL_DEFLINE, "lcl|NC_133200.1_cds_NP_001001398.2_93141"
        )
        self.assertEqual((gene, protein, symbol), ("GeneID:30560", "NP_001001398.2", "fgf6a"))

    def test_geneid_from_multi_source_db_xref(self):
        line = "lcl|x_cds_NP_1.1_1 [gene=abc] [db_xref=ZFIN:ZDB-1,GeneID:12345] [protein_id=NP_1.1]"
        self.assertEqual(li.parse_refseq_cds_defline(line, "x")[0], "GeneID:12345")

    def test_repeated_db_xref_fields_accumulate(self):
        line = "lcl|x_cds_NP_1.1_1 [db_xref=ZFIN:Z-1] [db_xref=GeneID:999] [protein_id=NP_1.1]"
        self.assertEqual(li.parse_refseq_cds_defline(line, "x")[0], "GeneID:999")

    def test_location_commas_do_not_swallow_later_fields(self):
        # [location=join(1..2,3..4)] must not consume [protein_id=...] after it.
        self.assertEqual(li.parse_refseq_cds_defline(REAL_DEFLINE, "x")[1], "NP_001001398.2")

    def test_protein_id_falls_back_to_lcl_record_id(self):
        rid = "lcl|NC_133200.1_cds_NP_001001398.2_93141"
        line = f"{rid} [gene=fgf6a] [db_xref=GeneID:30560]"
        self.assertEqual(li.parse_refseq_cds_defline(line, rid)[1], "NP_001001398.2")

    def test_gene_key_falls_back_to_locus_tag_then_symbol(self):
        self.assertEqual(
            li.parse_refseq_cds_defline("x [locus_tag=ABC_001] [protein_id=NP_1.1]", "x")[0],
            "locus_tag:ABC_001")
        self.assertEqual(
            li.parse_refseq_cds_defline("x [gene=sym] [protein_id=NP_1.1]", "x")[0], "gene:sym")

    def test_pseudogene_has_no_protein_id(self):
        line = "lcl|x_cds_1 [gene=psg] [db_xref=GeneID:77] [pseudo=true]"
        gene, protein, _ = li.parse_refseq_cds_defline(line, "lcl|x_cds_1")
        self.assertEqual(gene, "GeneID:77")
        self.assertIsNone(protein)


class TestRegexParsing(unittest.TestCase):
    def test_default_pattern_matches_ensembl_style(self):
        import re
        pattern = re.compile(li.DEFAULT_REGEX)
        self.assertEqual(
            li.gene_key_from_regex("ENSDARP1 pep gene:ENSDARG00000086269.4 x", pattern),
            "gene:ENSDARG00000086269.4")

    def test_no_match_returns_none(self):
        import re
        self.assertIsNone(li.gene_key_from_regex("no identifier here", re.compile(li.DEFAULT_REGEX)))


class TestAddSuffix(unittest.TestCase):
    def test_only_the_final_extension_is_separated(self):
        # RefSeq accessions carry their own dots; splitting on the first would mangle them.
        for name, want in [
            ("GCF_964276395.1_fGasAcu3.hap1.1_protein.faa",
             "GCF_964276395.1_fGasAcu3.hap1.1_protein_longest.faa"),
            ("GCF_013347855.1_fAngAng1.pri_cds_from_genomic.fna",
             "GCF_013347855.1_fAngAng1.pri_cds_from_genomic_longest.fna"),
        ]:
            self.assertEqual(li.add_suffix(Path("d") / name, "longest").name, want)


def write(path: Path, records):
    path.write_text("".join(f">{h}\n{s}\n" for h, s in records))


def cds_defline(pid, geneid, symbol):
    return (f"lcl|NC_1.1_cds_{pid}_1 [gene={symbol}] [db_xref=GeneID:{geneid}] "
            f"[protein_id={pid}] [location=1..9] [gbkey=CDS]")


def ids(path: Path):
    return [l[1:].strip() for l in path.read_text().splitlines() if l.startswith(">")]


class TestRefseqMode(unittest.TestCase):
    def run_script(self, cds, faa, *extra):
        tmp = Path(tempfile.mkdtemp())
        write(tmp / "x_cds_from_genomic.fna", cds)
        write(tmp / "x_protein.faa", faa)
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--refseq",
             "--input", str(tmp / "x_protein.faa"), "--cds", str(tmp / "x_cds_from_genomic.fna"),
             "--report", str(tmp / "r.tsv"), *extra],
            capture_output=True, text=True)
        return result, tmp

    def test_longest_cds_wins_and_default_names_get_the_suffix(self):
        cds = [(cds_defline("NP_1.1", "100", "geneA"), "ATGAAATAA"),
               (cds_defline("NP_2.1", "100", "geneA"), "ATGAAAAAAAAATAA"),
               (cds_defline("NP_3.1", "200", "geneB"), "ATGTGGTAA")]
        faa = [("NP_1.1 s", "MK"), ("NP_2.1 l", "MKKK"), ("NP_3.1 o", "MW")]
        result, tmp = self.run_script(cds, faa)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(ids(tmp / "x_protein_longest.faa"), ["NP_2.1", "NP_3.1"])
        self.assertEqual(ids(tmp / "x_cds_from_genomic_longest.fna"),
                         ids(tmp / "x_protein_longest.faa"))
        self.assertNotIn(" ", (tmp / "x_protein_longest.faa").read_text().splitlines()[0])

    def test_cds_keeps_its_stop_codon(self):
        _, tmp = self.run_script([(cds_defline("NP_1.1", "100", "g"), "ATGAAATAA")], [("NP_1.1 x", "MK")])
        self.assertTrue((tmp / "x_cds_from_genomic_longest.fna").read_text().strip().endswith("TAA"))

    def test_gene_survives_when_longest_isoform_has_no_protein(self):
        cds = [(cds_defline("NP_1.1", "100", "g"), "ATGAAATAA"),
               (cds_defline("NP_2.1", "100", "g"), "ATGAAAAAAAAATAA")]
        result, tmp = self.run_script(cds, [("NP_1.1 s", "MK")])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(ids(tmp / "x_protein_longest.faa"), ["NP_1.1"])

    def test_pseudogene_skipped(self):
        cds = [("lcl|NC_1.1_cds_1 [gene=p] [db_xref=GeneID:900] [pseudo=true]", "ATGTAA"),
               (cds_defline("NP_1.1", "100", "g"), "ATGAAATAA")]
        result, tmp = self.run_script(cds, [("NP_1.1 x", "MK")])
        self.assertEqual(ids(tmp / "x_protein_longest.faa"), ["NP_1.1"])
        self.assertIn("no protein_id", result.stderr)

    def test_duplicate_accession_emitted_once(self):
        cds = [(cds_defline("NP_1.1", "100", "a"), "ATGAAATAA"),
               (cds_defline("NP_1.1", "101", "b"), "ATGAAATAA")]
        result, tmp = self.run_script(cds, [("NP_1.1 x", "MK")])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(ids(tmp / "x_protein_longest.faa"), ["NP_1.1"])

    def test_translation_mismatch_kept_by_default(self):
        result, tmp = self.run_script([(cds_defline("NP_1.1", "100", "g"), "ATGAAATAA")],
                                      [("NP_1.1 wrong", "MWWW")])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("did not translate", result.stderr)
        self.assertEqual(ids(tmp / "x_protein_longest.faa"), ["NP_1.1"])

    def test_require_translation_fails(self):
        result, _ = self.run_script([(cds_defline("NP_1.1", "100", "g"), "ATGAAATAA")],
                                    [("NP_1.1 wrong", "MWWW")], "--require-translation")
        self.assertEqual(result.returncode, 1)

    def test_drop_untranslatable(self):
        cds = [(cds_defline("NP_1.1", "100", "a"), "ATGAAATAA"),
               (cds_defline("NP_2.1", "200", "b"), "ATGTGGTAA")]
        faa = [("NP_1.1 wrong", "MWWW"), ("NP_2.1 right", "MW")]
        result, tmp = self.run_script(cds, faa, "--drop-untranslatable")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(ids(tmp / "x_protein_longest.faa"), ["NP_2.1"])

    def test_refseq_without_cds_is_rejected(self):
        tmp = Path(tempfile.mkdtemp())
        write(tmp / "p.faa", [("NP_1.1 x", "MK")])
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--refseq", "--input", str(tmp / "p.faa")],
            capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("--refseq needs --cds", result.stderr)

    def test_report_columns(self):
        cds = [(cds_defline("NP_1.1", "100", "geneA"), "ATGAAATAA"),
               (cds_defline("NP_2.1", "100", "geneA"), "ATGAAAAAAAAATAA")]
        _, tmp = self.run_script(cds, [("NP_1.1 s", "MK"), ("NP_2.1 l", "MKKK")])
        lines = (tmp / "r.tsv").read_text().splitlines()
        self.assertEqual(lines[0].split("\t"),
                         ["gene_key", "gene_symbol", "protein_id", "length", "n_isoforms",
                          "translation_ok"])
        row = lines[1].split("\t")
        self.assertEqual([row[0], row[2], row[4], row[5]], ["GeneID:100", "NP_2.1", "2", "yes"])
        self.assertNotIn("\r", (tmp / "r.tsv").read_text())


class TestRegexMode(unittest.TestCase):
    """The original Ensembl-style behaviour: longest protein per gene."""

    def run_script(self, faa, *extra, cds=None):
        tmp = Path(tempfile.mkdtemp())
        write(tmp / "prot.faa", faa)
        cmd = [sys.executable, str(SCRIPT), "--input", str(tmp / "prot.faa")]
        if cds is not None:
            write(tmp / "cds.fna", cds)
            cmd += ["--cds", str(tmp / "cds.fna")]
        result = subprocess.run(cmd + list(extra), capture_output=True, text=True)
        return result, tmp

    def test_longest_protein_per_gene(self):
        faa = [("P1 gene:G1", "MKK"), ("P2 gene:G1", "MKKKKK"), ("P3 gene:G2", "MW")]
        result, tmp = self.run_script(faa)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(ids(tmp / "prot_longest.faa"), ["P2", "P3"])

    def test_records_without_the_pattern_are_skipped(self):
        faa = [("P1 gene:G1", "MKK"), ("P2 no-identifier", "MKKKKK")]
        result, tmp = self.run_script(faa)
        self.assertEqual(ids(tmp / "prot_longest.faa"), ["P1"])
        self.assertIn("no gene identifier", result.stderr)

    def test_custom_regex(self):
        faa = [("P1 locus=ABC", "MK"), ("P2 locus=ABC", "MKKK")]
        result, tmp = self.run_script(faa, "--regex", r"locus=\S+")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(ids(tmp / "prot_longest.faa"), ["P2"])

    def test_optional_cds_pairing_by_exact_id(self):
        faa = [("P1 gene:G1", "MK"), ("P2 gene:G1", "MKKK")]
        cds = [("P1", "ATGAAATAA"), ("P2", "ATGAAAAAAAAATAA")]
        result, tmp = self.run_script(faa, cds=cds)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(ids(tmp / "prot_longest.faa"), ["P2"])
        self.assertEqual(ids(tmp / "cds_longest.fna"), ["P2"])

    def test_no_output_written_when_nothing_resolves(self):
        result, _ = self.run_script([("P1 no-identifier", "MK")])
        self.assertEqual(result.returncode, 1)
        self.assertIn("No gene could be resolved", result.stderr)


class TestAllMode(unittest.TestCase):
    def build_tree(self):
        root = Path(tempfile.mkdtemp())
        for species, acc in [("zebrafish", "GCF_049306965.2_GRCz12tu"),
                             ("threespine", "GCF_964276395.1_fGasAcu3.hap1.1")]:
            d = root / species
            d.mkdir()
            write(d / f"{acc}_cds_from_genomic.fna",
                  [(cds_defline("NP_1.1", "100", "a"), "ATGAAATAA"),
                   (cds_defline("NP_2.1", "100", "a"), "ATGAAAAAAAAATAA")])
            write(d / f"{acc}_protein.faa", [("NP_1.1 s", "MK"), ("NP_2.1 l", "MKKK")])
        return root

    def test_processes_each_subdirectory(self):
        root = self.build_tree()
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--refseq", "--all", str(root)],
            capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue((root / "zebrafish" /
                         "GCF_049306965.2_GRCz12tu_protein_longest.faa").exists())
        # The hap1.1 name is the one a first-dot split would mangle.
        self.assertTrue((root / "threespine" /
                         "GCF_964276395.1_fGasAcu3.hap1.1_cds_from_genomic_longest.fna").exists())
        self.assertTrue((root / "zebrafish" / "zebrafish_isoforms.tsv").exists())
        self.assertIn("Processed 2 species", result.stderr)

    def test_species_missing_a_file_is_skipped_not_fatal(self):
        root = self.build_tree()
        (root / "empty").mkdir()
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--refseq", "--all", str(root)],
            capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("SKIP empty", result.stderr)
        self.assertIn("Processed 2 species (1 skipped, 0 failed)", result.stderr)

    def test_one_bad_species_does_not_abandon_the_rest(self):
        root = self.build_tree()
        bad = root / "broken"
        bad.mkdir()
        (bad / "a_cds_from_genomic.fna").write_text("")
        (bad / "a_protein.faa").write_text("")
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--refseq", "--all", str(root)],
            capture_output=True, text=True)
        self.assertEqual(result.returncode, 1)          # status still reflects it
        self.assertIn("FAILED broken", result.stderr)
        self.assertIn("Processed 2 species (0 skipped, 1 failed)", result.stderr)
        self.assertTrue((root / "zebrafish" /
                         "GCF_049306965.2_GRCz12tu_protein_longest.faa").exists())


if __name__ == "__main__":
    unittest.main()
