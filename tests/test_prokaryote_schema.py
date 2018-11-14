""" Tests of the knowledge base schema for prokaryotes

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Author: Bilal Shaikh <bilal.shaikh@columbia.edu>
:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2018-02-07
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_kb import core, prokaryote_schema
from wc_utils.util import chem
import Bio.Alphabet
import Bio.Seq
import Bio.SeqIO
import Bio.SeqUtils
import mendeleev
import unittest


class CellTestCase(unittest.TestCase):
    def test_constructor(self):
        cell = core.Cell(taxon=2104)

        self.assertEqual(cell.knowledge_base, None)
        self.assertEqual(cell.taxon, 2104)
        self.assertEqual(cell.observables, [])
        self.assertEqual(cell.species_types, [])
        self.assertEqual(cell.compartments, [])
        self.assertEqual(cell.reactions, [])
        self.assertEqual(cell.loci, [])

        self.assertEqual(cell.species_types.get(
            __type=core.DnaSpeciesType), [])
        self.assertEqual(cell.loci.get(__type=prokaryote_schema.PromoterLocus), [])


class RnaSpeciesTypeTestCase(unittest.TestCase):

    def test_constructor(self):
        dna1 = core.DnaSpeciesType(id='dna1', seq=Bio.Seq.Seq(
            'ACGTACGTACGTACG', alphabet=Bio.Alphabet.DNAAlphabet()))
        tu1 = prokaryote_schema.TranscriptionUnitLocus(
            id='tu1', polymer=dna1, start=1, end=15)
        rna1 = prokaryote_schema.RnaSpeciesType(id='rna1', name='rna1', transcription_units=[
                                   tu1], type=1, half_life=2)

        self.assertEqual(rna1.id, 'rna1')
        self.assertEqual(rna1.name, 'rna1')
        self.assertEqual(rna1.transcription_units, [tu1])
        self.assertEqual(rna1.type, 1)
        self.assertEqual(rna1.half_life, 2)

    def test_get_empirical_formula(self):
        dna1 = core.DnaSpeciesType(id='dna1', seq=Bio.Seq.Seq(
            'A', alphabet=Bio.Alphabet.DNAAlphabet()))
        tu1 = prokaryote_schema.TranscriptionUnitLocus(
            id='tu1', polymer=dna1, start=1, end=1)
        rna1 = prokaryote_schema.RnaSpeciesType(
            id='rna1', name='rna1', transcription_units=[tu1])
        self.assertEqual(rna1.get_empirical_formula(),
                         chem.EmpiricalFormula('C10H12N5O7P'))

        dna1 = core.DnaSpeciesType(id='dna1', seq=Bio.Seq.Seq(
            'C', alphabet=Bio.Alphabet.DNAAlphabet()))
        tu1 = prokaryote_schema.TranscriptionUnitLocus(
            id='tu1', polymer=dna1, start=1, end=1)
        rna1 = prokaryote_schema.RnaSpeciesType(
            id='rna1', name='rna1', transcription_units=[tu1])
        self.assertEqual(rna1.get_empirical_formula(),
                         chem.EmpiricalFormula('C9H12N3O8P'))

        dna1 = core.DnaSpeciesType(id='dna1', seq=Bio.Seq.Seq(
            'G', alphabet=Bio.Alphabet.DNAAlphabet()))
        tu1 = prokaryote_schema.TranscriptionUnitLocus(
            id='tu1', polymer=dna1, start=1, end=1)
        rna1 = prokaryote_schema.RnaSpeciesType(
            id='rna1', name='rna1', transcription_units=[tu1])
        self.assertEqual(rna1.get_empirical_formula(),
                         chem.EmpiricalFormula('C10H12N5O8P'))

        dna1 = core.DnaSpeciesType(id='dna1', seq=Bio.Seq.Seq(
            'T', alphabet=Bio.Alphabet.DNAAlphabet()))
        tu1 = prokaryote_schema.TranscriptionUnitLocus(
            id='tu1', polymer=dna1, start=1, end=1)
        rna1 = prokaryote_schema.RnaSpeciesType(
            id='rna1', name='rna1', transcription_units=[tu1])
        self.assertEqual(rna1.get_empirical_formula(),
                         chem.EmpiricalFormula('C9H11N2O9P'))

        dna1 = core.DnaSpeciesType(id='dna1', seq=Bio.Seq.Seq(
            'AAAA', alphabet=Bio.Alphabet.DNAAlphabet()))
        tu1 = prokaryote_schema.TranscriptionUnitLocus(
            id='tu1', polymer=dna1, start=1, end=2)
        rna1 = prokaryote_schema.RnaSpeciesType(
            id='rna1', name='rna1', transcription_units=[tu1])
        self.assertEqual(rna1.get_empirical_formula(),
                         chem.EmpiricalFormula('C20H23N10O13P2'))

    def test_get_charge(self):
        dna1 = core.DnaSpeciesType(id='dna1', seq=Bio.Seq.Seq(
            'AAAA', alphabet=Bio.Alphabet.DNAAlphabet()))
        tu1 = prokaryote_schema.TranscriptionUnitLocus(
            id='tu1', polymer=dna1, start=1, end=1)
        rna1 = prokaryote_schema.RnaSpeciesType(
            id='rna1', name='rna1', transcription_units=[tu1])
        self.assertEqual(rna1.get_charge(), -2)

        dna1 = core.DnaSpeciesType(id='dna1', seq=Bio.Seq.Seq(
            'AAAA', alphabet=Bio.Alphabet.DNAAlphabet()))
        tu1 = prokaryote_schema.TranscriptionUnitLocus(
            id='tu1', polymer=dna1, start=1, end=2)
        rna1 = prokaryote_schema.RnaSpeciesType(
            id='rna1', name='rna1', transcription_units=[tu1])
        self.assertEqual(rna1.get_charge(), -3)

    def test_get_mol_wt(self):
        dna1 = core.DnaSpeciesType(id='dna1', seq=Bio.Seq.Seq(
            'AACCGGTT', alphabet=Bio.Alphabet.DNAAlphabet()))
        tu1 = prokaryote_schema.TranscriptionUnitLocus(
            id='tu1', polymer=dna1, start=1, end=1)
        rna1 = prokaryote_schema.RnaSpeciesType(
            id='rna1', name='rna1', transcription_units=[tu1])
        exp_mol_wt = \
            + Bio.SeqUtils.molecular_weight(rna1.get_seq()) \
            - (rna1.get_len() + 1) * mendeleev.element('H').atomic_weight
        self.assertAlmostEqual(rna1.get_mol_wt(), exp_mol_wt, places=1)

        dna1 = core.DnaSpeciesType(id='dna1', seq=Bio.Seq.Seq(
            'AACCGGTT', alphabet=Bio.Alphabet.DNAAlphabet()))
        tu1 = prokaryote_schema.TranscriptionUnitLocus(
            id='tu1', polymer=dna1, start=3, end=3)
        rna1 = prokaryote_schema.RnaSpeciesType(
            id='rna1', name='rna1', transcription_units=[tu1])
        exp_mol_wt = \
            + Bio.SeqUtils.molecular_weight(rna1.get_seq()) \
            - (rna1.get_len() + 1) * mendeleev.element('H').atomic_weight
        self.assertAlmostEqual(rna1.get_mol_wt(), exp_mol_wt, places=1)

        dna1 = core.DnaSpeciesType(id='dna1', seq=Bio.Seq.Seq(
            'AACCGGTT', alphabet=Bio.Alphabet.DNAAlphabet()))
        tu1 = prokaryote_schema.TranscriptionUnitLocus(
            id='tu1', polymer=dna1, start=5, end=5)
        rna1 = prokaryote_schema.RnaSpeciesType(
            id='rna1', name='rna1', transcription_units=[tu1])
        exp_mol_wt = \
            + Bio.SeqUtils.molecular_weight(rna1.get_seq()) \
            - (rna1.get_len() + 1) * mendeleev.element('H').atomic_weight
        self.assertAlmostEqual(rna1.get_mol_wt(), exp_mol_wt, places=1)


class ProteinSpeciesTypeTestCase(unittest.TestCase):
    def setUp(self):
        # Mycoplasma Genintalium Genome
        records = Bio.SeqIO.parse('tests/fixtures/seq.fna', 'fasta')
        dna_seq = next(records).seq
        dna1 = core.DnaSpeciesType(seq=dna_seq)

        cell1 = dna1.cell = core.Cell()
        cell1.knowledge_base = core.KnowledgeBase(
            translation_table=4)  # Table 4 is for mycoplasma

        # MPN001
        gene1 = prokaryote_schema.GeneLocus(id='gene1', cell=cell1,
                               polymer=dna1, start=692, end=1834)
        tu1 = prokaryote_schema.TranscriptionUnitLocus(
            id='tu1', genes=[gene1], polymer=dna1)
        self.prot1 = prokaryote_schema.ProteinSpeciesType(
            id='prot1', gene=gene1, cell=cell1)

        # MPN011
        gene2 = prokaryote_schema.GeneLocus(id='gene2', cell=cell1, polymer=dna1,
                               start=12838, end=13533, strand=core.PolymerStrand.negative)
        tu2 = prokaryote_schema.TranscriptionUnitLocus(
            id='tu2', genes=[gene2], polymer=dna1)
        self.prot2 = prokaryote_schema.ProteinSpeciesType(
            id='prot2', gene=gene2, cell=cell1)

    def test_constructor(self):
        protein = prokaryote_schema.ProteinSpeciesType(
            id='prot1', name='prot1', half_life=2)
        # attribute_order = ('id', 'cell', 'name', 'gene', 'rna', 'half_life', 'comments')

        self.assertEqual(protein.id, 'prot1')
        self.assertEqual(protein.name, 'prot1')
        self.assertEqual(protein.half_life, 2)
        self.assertEqual(protein.cell, None)

    def test_get_seq(self):
        # Use translation table 4 since example genes are from 
        # Mycoplasma genitallium

        # MPN001
        self.assertEqual(self.prot1.get_seq()[0:10], 'MKVLINKNEL')
        self.assertEqual(self.prot1.get_seq()[-10:], 'ELKEILVPSK')

        # MPN011
        self.assertEqual(self.prot2.get_seq()[0:10], 'MKFKFLLTPL')
        self.assertEqual(self.prot2.get_seq()[-10:], 'LFRYLVYLIE')

    def test_get_empirical_formula(self):
        # MPN001
        self.assertEqual(self.prot1.get_empirical_formula(),
                         chem.EmpiricalFormula('C1980H3146N510O596S7'))
        # MPN011
        self.assertEqual(self.prot2.get_empirical_formula(),
                         chem.EmpiricalFormula('C1246H1928N306O352S3'))

    def test_get_mol_wt(self):
        # MPN001
        self.assertAlmostEqual(self.prot1.get_mol_wt(), 43856.342, delta=0.3)
        # MNP011
        self.assertAlmostEqual(self.prot2.get_mol_wt(), 26923.100, delta=0.3)

    def test_get_charge(self):
        self.assertEqual(self.prot1.get_charge(), 1)

        self.assertEqual(self.prot2.get_charge(), 12)


class PromoterLocusTestCase(unittest.TestCase):
    def test_constructor(self):
        promoter = prokaryote_schema.PromoterLocus(
            id='promoter1', name='promoter1', pribnow_start=1, pribnow_end=2)
        self.assertEqual(promoter.id, 'promoter1')
        self.assertEqual(promoter.name, 'promoter1')
        self.assertEqual(promoter.pribnow_start, 1)
        self.assertEqual(promoter.pribnow_end, 2)


class GeneLocusTestCase(unittest.TestCase):
    def test(self):
        gene = prokaryote_schema.GeneLocus(id='gene1', name='gene1',
                              symbol='gene_1', start=1, end=2)
        self.assertEqual(gene.id, 'gene1')
        self.assertEqual(gene.name, 'gene1')
        self.assertEqual(gene.symbol, 'gene_1')
        self.assertEqual(gene.start, 1)
        self.assertEqual(gene.end, 2)


class TranscriptionUnitLocusTestCase(unittest.TestCase):
    def test(self):
        dna1 = core.DnaSpeciesType(id='dna1', seq=Bio.Seq.Seq('ACGTACGTACGTACG', alphabet=Bio.Alphabet.DNAAlphabet()),
                                   circular=False, double_stranded=False)

        tu1 = prokaryote_schema.TranscriptionUnitLocus(
            id='tu1', name='tu1', polymer=dna1, strand=core.PolymerStrand.positive, start=1, end=15)

        # test constructor
        self.assertEqual(tu1.id, 'tu1')
        self.assertEqual(tu1.name, 'tu1')
        self.assertEqual(tu1.polymer, dna1)
        self.assertEqual(tu1.strand, core.PolymerStrand.positive)
        self.assertEqual(tu1.start, 1)
        self.assertEqual(tu1.end, 15)

        # test methods
        self.assertEqual(tu1.get_3_prime(), 15)
        self.assertEqual(tu1.get_5_prime(), 1)

        # flip strand; test methods
        rev_comp_seq = tu1.get_seq().reverse_complement()
        tu1.strand = core.PolymerStrand.negative
        self.assertEqual(tu1.get_3_prime(), 1)
        self.assertEqual(tu1.get_5_prime(), 15)


class ComplexSpeciesTypeTestCase(unittest.TestCase):
    def test_ComplexSpeciesType(self):

        # Test constructor
        complex1 = core.ComplexSpeciesType()

        # Generate test proteins from  Mycoplasma Genintalium Genome
        records = Bio.SeqIO.parse('tests/fixtures/seq.fna', 'fasta')
        dna_seq = next(records).seq
        dna1 = core.DnaSpeciesType(seq=dna_seq)

        cell1 = dna1.cell = core.Cell()
        cell1.knowledge_base = core.KnowledgeBase(
            translation_table=4)  # Table 4 is for mycoplasma

        # Protein 1,  MPN001
        gene1 = prokaryote_schema.GeneLocus(id='gene1', cell=cell1,
                               polymer=dna1, start=692, end=1834)
        tu1 = prokaryote_schema.TranscriptionUnitLocus(
            id='tu1', genes=[gene1], polymer=dna1)
        prot1 = prokaryote_schema.ProteinSpeciesType(
            id='prot1', gene=gene1, cell=cell1)

        # Protein 2, MPN011
        gene2 = prokaryote_schema.GeneLocus(id='gene2', cell=cell1, polymer=dna1,
                               start=12838, end=13533, strand=core.PolymerStrand.negative)
        tu2 = prokaryote_schema.TranscriptionUnitLocus(
            id='tu2', genes=[gene2], polymer=dna1)
        prot2 = prokaryote_schema.ProteinSpeciesType(
            id='prot2', gene=gene2, cell=cell1)

        # Test adding complexation
        # Add formation reaction: (2) prot1 + (3) prot2 ==> complex1
        species_coeff1 = core.SpeciesTypeCoefficient(
            species_type=prot1, coefficient=2)
        species_coeff2 = core.SpeciesTypeCoefficient(
            species_type=prot2, coefficient=3)
        complex1.subunits = [species_coeff1, species_coeff2]

        self.assertEqual(complex1.get_charge(), 38)
        self.assertAlmostEqual(complex1.get_mol_wt(),
                               (2*prot1.get_mol_wt() + 3 * prot2.get_mol_wt()))
        self.assertEqual(complex1.get_empirical_formula(),
                         chem.EmpiricalFormula('C7698H12076N1938O2248S23'))
