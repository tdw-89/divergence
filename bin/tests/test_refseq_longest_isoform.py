"""Unit tests for bin/refseq_longest_isoform.py."""

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
import refseq_longest_isoform as rli

SCRIPT = Path(__file__).parent.parent / "refseq_longest_isoform.py"

# A real defline, wrapped here only for readability.
REAL_DEFLINE = (
    "lcl|NC_133200.1_cds_NP_001001398.2_93141 [gene=fgf6a] "
    "[db_xref=GeneID:30560] [protein=fibroblast growth factor 6a precursor] "
    "[protein_id=NP_001001398.2] "
    "[location=join(17881715..17882051,17882586..17882689,17887110..17887286)] [gbkey=CDS]"
)


class TestDeflineParsing(unittest.TestCase):
    def test_real_defline(self):
        gene, protein, symbol = rli.parse_cds_defline(
            REAL_DEFLINE, "lcl|NC_133200.1_cds_NP_001001398.2_93141"
        )
        self.assertEqual(gene, "GeneID:30560")
        self.assertEqual(protein, "NP_001001398.2")
        self.assertEqual(symbol, "fgf6a")

    def test_geneid_extracted_from_multi_source_db_xref(self):
        line = "lcl|x_cds_NP_1.1_1 [gene=abc] [db_xref=ZFIN:ZDB-GENE-1,GeneID:12345] [protein_id=NP_1.1]"
        gene, protein, _ = rli.parse_cds_defline(line, "lcl|x_cds_NP_1.1_1")
        self.assertEqual(gene, "GeneID:12345")
        self.assertEqual(protein, "NP_1.1")

    def test_repeated_db_xref_fields_are_accumulated(self):
        line = "lcl|x_cds_NP_1.1_1 [gene=abc] [db_xref=ZFIN:ZDB-1] [db_xref=GeneID:999] [protein_id=NP_1.1]"
        gene, _, _ = rli.parse_cds_defline(line, "lcl|x_cds_NP_1.1_1")
        self.assertEqual(gene, "GeneID:999")

    def test_location_value_with_commas_does_not_break_fields(self):
        # [location=join(1..2,3..4)] must not swallow the fields after it.
        _, protein, _ = rli.parse_cds_defline(REAL_DEFLINE, "irrelevant")
        self.assertEqual(protein, "NP_001001398.2")

    def test_protein_id_falls_back_to_lcl_record_id(self):
        rid = "lcl|NC_133200.1_cds_NP_001001398.2_93141"
        _, protein, _ = rli.parse_cds_defline(f"{rid} [gene=fgf6a] [db_xref=GeneID:30560]", rid)
        self.assertEqual(protein, "NP_001001398.2")

    def test_gene_key_falls_back_to_locus_tag_then_symbol(self):
        line = "lcl|x_cds_NP_1.1_1 [locus_tag=ABC_001] [protein_id=NP_1.1]"
        self.assertEqual(rli.parse_cds_defline(line, "x")[0], "locus_tag:ABC_001")
        line = "lcl|x_cds_NP_1.1_1 [gene=onlysymbol] [protein_id=NP_1.1]"
        self.assertEqual(rli.parse_cds_defline(line, "x")[0], "gene:onlysymbol")

    def test_pseudogene_without_protein_id_yields_none(self):
        line = "lcl|x_cds_1 [gene=psg] [db_xref=GeneID:77] [pseudo=true]"
        gene, protein, _ = rli.parse_cds_defline(line, "lcl|x_cds_1")
        self.assertEqual(gene, "GeneID:77")
        self.assertIsNone(protein)


def write(path: Path, records: list[tuple[str, str]]) -> None:
    path.write_text("".join(f">{h}\n{s}\n" for h, s in records))


def cds_defline(protein_id: str, geneid: str, symbol: str) -> str:
    return (
        f"lcl|NC_1.1_cds_{protein_id}_1 [gene={symbol}] [db_xref=GeneID:{geneid}] "
        f"[protein_id={protein_id}] [location=1..9] [gbkey=CDS]"
    )


class TestSelection(unittest.TestCase):
    """End-to-end runs of the script over small synthetic RefSeq-shaped files."""

    def run_script(self, cds, faa, *extra):
        tmp = Path(tempfile.mkdtemp())
        write(tmp / "cds.fna", cds)
        write(tmp / "prot.faa", faa)
        result = subprocess.run(
            [sys.executable, str(SCRIPT),
             "--cds", str(tmp / "cds.fna"), "--proteins", str(tmp / "prot.faa"),
             "--out-proteins", str(tmp / "out.faa"), "--out-cds", str(tmp / "out.fna"),
             "--report", str(tmp / "report.tsv"), *extra],
            capture_output=True, text=True,
        )
        return result, tmp

    @staticmethod
    def ids(path: Path) -> list[str]:
        return [l[1:].strip() for l in path.read_text().splitlines() if l.startswith(">")]

    def test_longest_isoform_wins_and_ids_match(self):
        cds = [
            (cds_defline("NP_1.1", "100", "geneA"), "ATGAAATAA"),
            (cds_defline("NP_2.1", "100", "geneA"), "ATGAAAAAAAAATAA"),  # longest
            (cds_defline("NP_3.1", "200", "geneB"), "ATGTGGTAA"),
        ]
        faa = [("NP_1.1 short [Danio rerio]", "MK"),
               ("NP_2.1 long [Danio rerio]", "MKKK"),
               ("NP_3.1 other [Danio rerio]", "MW")]
        result, tmp = self.run_script(cds, faa)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.ids(tmp / "out.faa"), ["NP_2.1", "NP_3.1"])
        # The pipeline matches CDS to protein on exact ID equality.
        self.assertEqual(self.ids(tmp / "out.fna"), self.ids(tmp / "out.faa"))
        self.assertNotIn("Danio rerio", (tmp / "out.faa").read_text())

    def test_cds_keeps_its_stop_codon(self):
        cds = [(cds_defline("NP_1.1", "100", "geneA"), "ATGAAATAA")]
        faa = [("NP_1.1 x", "MK")]
        _, tmp = self.run_script(cds, faa)
        self.assertTrue((tmp / "out.fna").read_text().strip().endswith("TAA"))

    def test_gene_survives_when_longest_isoform_has_no_protein(self):
        # NP_2.1 is longest but absent from the .faa; the gene must fall back
        # to NP_1.1 rather than disappearing.
        cds = [
            (cds_defline("NP_1.1", "100", "geneA"), "ATGAAATAA"),
            (cds_defline("NP_2.1", "100", "geneA"), "ATGAAAAAAAAATAA"),
        ]
        faa = [("NP_1.1 short", "MK")]
        result, tmp = self.run_script(cds, faa)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.ids(tmp / "out.faa"), ["NP_1.1"])

    def test_pseudogene_records_are_skipped(self):
        cds = [
            ("lcl|NC_1.1_cds_1 [gene=psg] [db_xref=GeneID:900] [pseudo=true]", "ATGTAA"),
            (cds_defline("NP_1.1", "100", "geneA"), "ATGAAATAA"),
        ]
        faa = [("NP_1.1 x", "MK")]
        result, tmp = self.run_script(cds, faa)
        self.assertEqual(self.ids(tmp / "out.faa"), ["NP_1.1"])
        self.assertIn("no protein_id", result.stderr)

    def test_duplicate_protein_accession_emitted_once(self):
        # Same protein under two gene keys (alt locus): IDs must stay unique,
        # since ks.py treats a duplicate CDS ID as a hard error.
        cds = [
            (cds_defline("NP_1.1", "100", "geneA"), "ATGAAATAA"),
            (cds_defline("NP_1.1", "101", "geneA_alt"), "ATGAAATAA"),
        ]
        faa = [("NP_1.1 x", "MK")]
        result, tmp = self.run_script(cds, faa)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.ids(tmp / "out.faa"), ["NP_1.1"])

    def test_translation_mismatch_reported_but_kept_by_default(self):
        cds = [(cds_defline("NP_1.1", "100", "geneA"), "ATGAAATAA")]  # MK
        faa = [("NP_1.1 wrong", "MWWW")]
        result, tmp = self.run_script(cds, faa)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("did not translate", result.stderr)
        self.assertEqual(self.ids(tmp / "out.faa"), ["NP_1.1"])

    def test_require_translation_fails(self):
        cds = [(cds_defline("NP_1.1", "100", "geneA"), "ATGAAATAA")]
        faa = [("NP_1.1 wrong", "MWWW")]
        result, _ = self.run_script(cds, faa, "--require-translation")
        self.assertEqual(result.returncode, 1)

    def test_drop_untranslatable_removes_the_pair(self):
        cds = [
            (cds_defline("NP_1.1", "100", "geneA"), "ATGAAATAA"),
            (cds_defline("NP_2.1", "200", "geneB"), "ATGTGGTAA"),
        ]
        faa = [("NP_1.1 wrong", "MWWW"), ("NP_2.1 right", "MW")]
        result, tmp = self.run_script(cds, faa, "--drop-untranslatable")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.ids(tmp / "out.faa"), ["NP_2.1"])

    def test_report_columns(self):
        cds = [(cds_defline("NP_1.1", "100", "geneA"), "ATGAAATAA"),
               (cds_defline("NP_2.1", "100", "geneA"), "ATGAAAAAAAAATAA")]
        faa = [("NP_1.1 s", "MK"), ("NP_2.1 l", "MKKK")]
        _, tmp = self.run_script(cds, faa)
        lines = (tmp / "report.tsv").read_text().splitlines()
        self.assertEqual(
            lines[0].split("\t"),
            ["gene_key", "gene_symbol", "protein_id", "cds_length", "n_isoforms", "translation_ok"],
        )
        row = lines[1].split("\t")
        self.assertEqual(row[0], "GeneID:100")
        self.assertEqual(row[2], "NP_2.1")
        self.assertEqual(row[4], "2")   # both isoforms counted
        self.assertEqual(row[5], "yes")
        # Written LF, not CRLF -- the same trap ks.py hit.
        self.assertNotIn("\r", (tmp / "report.tsv").read_text())


if __name__ == "__main__":
    unittest.main()
