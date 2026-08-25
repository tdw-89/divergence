"""Unit tests for bin/extract_paralogs.py."""

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
import extract_paralogs as ep


class TestSpeciesColumn(unittest.TestCase):
    HEADER = ["Orthogroup", "Danio_rerio", "Oryzias_latipes", "Danio_rerio_alt"]

    def test_prefix_match(self):
        self.assertEqual(ep.find_species_column(self.HEADER, "Oryzias_latipes"), 2)

    def test_missing_species_raises(self):
        with self.assertRaises(ValueError):
            ep.find_species_column(self.HEADER, "Gallus_gallus")

    def test_ambiguous_prefix_raises(self):
        # 'Danio_rerio' prefixes two columns, which must not be resolved silently.
        with self.assertRaises(ValueError):
            ep.find_species_column(self.HEADER, "Danio_rerio")


class TestMultipleSpeciesColumns(unittest.TestCase):
    HEADER = ["Orthogroup", "Danio_rerio", "Oryzias_latipes", "Lepisosteus_oculatus"]

    def test_single_name(self):
        self.assertEqual(ep.find_species_columns(self.HEADER, "Danio_rerio"),
                         {"Danio_rerio": 1})

    def test_comma_separated_list(self):
        self.assertEqual(
            ep.find_species_columns(self.HEADER, "Danio_rerio,Oryzias_latipes"),
            {"Danio_rerio": 1, "Oryzias_latipes": 2},
        )

    def test_whitespace_around_names_tolerated(self):
        self.assertEqual(
            ep.find_species_columns(self.HEADER, " Danio_rerio , Oryzias_latipes "),
            {"Danio_rerio": 1, "Oryzias_latipes": 2},
        )

    def test_all_selects_every_species(self):
        self.assertEqual(
            ep.find_species_columns(self.HEADER, "all"),
            {"Danio_rerio": 1, "Oryzias_latipes": 2, "Lepisosteus_oculatus": 3},
        )
        self.assertEqual(ep.find_species_columns(self.HEADER, "ALL"),
                         ep.find_species_columns(self.HEADER, "all"))

    def test_unknown_name_raises(self):
        with self.assertRaises(ValueError):
            ep.find_species_columns(self.HEADER, "Danio_rerio,Gallus_gallus")

    def test_empty_selection_raises(self):
        with self.assertRaises(ValueError):
            ep.find_species_columns(self.HEADER, " , ")

    def test_duplicate_selection_raises(self):
        with self.assertRaises(ValueError):
            ep.find_species_columns(self.HEADER, "Danio_rerio,Danio_rerio")


class TestGeneParsing(unittest.TestCase):
    def test_splits_on_comma_space(self):
        self.assertEqual(ep.parse_genes("g1, g2, g3"), ["g1", "g2", "g3"])

    def test_empty_cell(self):
        self.assertEqual(ep.parse_genes("   "), [])

    def test_single_gene(self):
        self.assertEqual(ep.parse_genes("g1"), ["g1"])

    def test_tolerates_missing_space_after_comma(self):
        self.assertEqual(ep.parse_genes("g1,g2"), ["g1", "g2"])


class TestCombinedTreeFile(unittest.TestCase):
    def test_reads_orthofinder3_layout(self):
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "Resolved_Gene_Trees.txt").write_text(
                "OG0000000: (a,(b,c));\nOG0000001: ((d,e),f);\n\n"
            )
            trees = ep.load_combined_trees(tmp)
        self.assertEqual(sorted(trees), ["OG0000000", "OG0000001"])
        self.assertEqual(trees["OG0000001"], "((d,e),f);")

    def test_missing_directory_is_not_fatal(self):
        self.assertEqual(ep.load_combined_trees(None), {})
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(ep.load_combined_trees(tmp), {})

    def test_ids_containing_colons_survive(self):
        # Only the first colon separates the orthogroup id from the newick.
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "Resolved_Gene_Trees.txt").write_text("OG0: (a:0.1,b:0.2);\n")
            trees = ep.load_combined_trees(tmp)
        self.assertEqual(trees["OG0"], "(a:0.1,b:0.2);")


class TestTreeLabelNormalisation(unittest.TestCase):
    SPECIES = ["Danio_rerio", "Oryzias_latipes"]

    def test_strips_species_prefixes(self):
        genes = {"gene1", "gene2", "gene3"}
        newick = "((Danio_rerio_gene1,Danio_rerio_gene2),Oryzias_latipes_gene3);"
        tree, unmatched = ep.normalise_tree(newick, genes, self.SPECIES)
        self.assertEqual(unmatched, [])
        self.assertEqual(sorted(t.name for t in tree.get_terminals()),
                         ["gene1", "gene2", "gene3"])

    def test_leaves_already_bare_labels_alone(self):
        genes = {"gene1", "gene2"}
        tree, unmatched = ep.normalise_tree("(gene1,gene2);", genes, self.SPECIES)
        self.assertEqual(unmatched, [])
        self.assertEqual(sorted(t.name for t in tree.get_terminals()), ["gene1", "gene2"])

    def test_reports_tips_that_do_not_belong(self):
        genes = {"gene1"}
        tree, unmatched = ep.normalise_tree("(Danio_rerio_gene1,mystery);", genes, self.SPECIES)
        self.assertEqual(unmatched, ["mystery"])

    def test_prefix_only_stripped_when_remainder_is_a_real_gene(self):
        # 'Danio_rerio_x' is itself a gene id here, so it must be left intact
        # rather than truncated to 'x'.
        genes = {"Danio_rerio_x"}
        tree, unmatched = ep.normalise_tree("(Danio_rerio_x,Danio_rerio_x);",
                                            genes, self.SPECIES)
        self.assertEqual(unmatched, [])
        self.assertEqual([t.name for t in tree.get_terminals()],
                         ["Danio_rerio_x", "Danio_rerio_x"])

    def test_handles_pipe_characters_in_ids(self):
        genes = {"gi|123|emb|ABC.1|"}
        newick = "(Danio_rerio_gi|123|emb|ABC.1|,Danio_rerio_gi|123|emb|ABC.1|);"
        tree, unmatched = ep.normalise_tree(newick, genes, self.SPECIES)
        self.assertEqual(unmatched, [])
        self.assertEqual(tree.get_terminals()[0].name, "gi|123|emb|ABC.1|")


class TestPerOrthogroupTreeFallback(unittest.TestCase):
    def test_finds_legacy_per_og_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "OG0000001_tree.txt").write_text("(a,b);")
            self.assertIsNotNone(ep.locate_tree(tmp, "OG0000001"))

    def test_returns_none_when_absent(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertIsNone(ep.locate_tree(tmp, "OG0000001"))


class TestSpeciesPrefixPattern(unittest.TestCase):
    """OrthoFinder spells the species prefix differently in the two files it writes.

    Orthogroups.tsv keeps the name verbatim; gene-tree tips have *some*
    punctuation replaced with underscores -- dots are rewritten, hyphens are
    not. A run whose species names contain dots (RefSeq accessions such as
    GCF_049306965.2_GRCz12tu) lost every gene tree, silently degrading dS to
    pairwise-only. So did the first fix, which over-corrected by rewriting
    hyphens too and so missed GCF_053564925.1_Olatipes_Hd-rR_3.1.
    """

    def test_plain_name_matches_itself(self):
        self.assertTrue(ep.prefix_pattern("Danio_rerio").match("Danio_rerio_geneA"))

    def test_dotted_name_matches_both_spellings(self):
        pattern = ep.prefix_pattern("GCF_049306965.2_GRCz12tu")
        self.assertTrue(pattern.match("GCF_049306965.2_GRCz12tu_XP_1.1"))
        self.assertTrue(pattern.match("GCF_049306965_2_GRCz12tu_XP_1.1"))

    def test_hyphen_survives_but_dots_are_rewritten(self):
        # The real medaka case: OrthoFinder rewrote the dots and kept the hyphen.
        pattern = ep.prefix_pattern("GCF_053564925.1_Olatipes_Hd-rR_3.1_protein_longest")
        self.assertTrue(
            pattern.match("GCF_053564925_1_Olatipes_Hd-rR_3_1_protein_longest_XP_078802286.1")
        )

    def test_prefix_must_be_followed_by_a_separator(self):
        self.assertIsNone(ep.prefix_pattern("Danio_rerio").match("Danio_rerioX"))

    def test_tips_with_sanitised_prefixes_are_matched(self):
        genes = {"gi|290752548|emb|CBH40520.1|", "gi|290752728|emb|CBH40702.1|"}
        newick = ("(GCF_049306965_2_GRCz12tu_protein_longest_gi|290752548|emb|CBH40520.1|:0.1,"
                  "GCF_049306965_2_GRCz12tu_protein_longest_gi|290752728|emb|CBH40702.1|:0.2);")
        tree, unmatched = ep.normalise_tree(
            newick, genes, ["GCF_049306965.2_GRCz12tu_protein_longest"])
        self.assertEqual(unmatched, [])
        self.assertEqual({t.name for t in tree.get_terminals()}, genes)

    def test_verbatim_prefixes_still_match(self):
        genes = {"geneA", "geneB"}
        newick = "(Danio_rerio_geneA:0.1,Danio_rerio_geneB:0.2);"
        tree, unmatched = ep.normalise_tree(newick, genes, ["Danio_rerio"])
        self.assertEqual(unmatched, [])
        self.assertEqual({t.name for t in tree.get_terminals()}, genes)

    def test_mixed_punctuation_tips_are_matched(self):
        genes = {"XP_078802286.1", "XP_078794068.1"}
        species = "GCF_053564925.1_Olatipes_Hd-rR_3.1_protein_longest"
        newick = ("(GCF_053564925_1_Olatipes_Hd-rR_3_1_protein_longest_XP_078802286.1:0.1,"
                  "GCF_053564925_1_Olatipes_Hd-rR_3_1_protein_longest_XP_078794068.1:0.2);")
        tree, unmatched = ep.normalise_tree(newick, genes, [species])
        self.assertEqual(unmatched, [])
        self.assertEqual({t.name for t in tree.get_terminals()}, genes)

    def test_ids_that_merely_look_prefixed_are_not_mangled(self):
        # The remainder must be a real gene, or the label is left alone.
        genes = {"Danio_rerio_geneA"}
        tree, unmatched = ep.normalise_tree(
            "(Danio_rerio_geneA:0.1,Danio_rerio_geneZ:0.2);", genes, ["Danio_rerio"])
        self.assertEqual(unmatched, ["Danio_rerio_geneZ"])
        self.assertIn("Danio_rerio_geneA", {t.name for t in tree.get_terminals()})


if __name__ == "__main__":
    unittest.main()
