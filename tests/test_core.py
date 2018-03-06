""" Tests of the knowledge base

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-02-07
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_kb import core
from wc_utils.util import chem
import Bio.Alphabet
import Bio.Seq
import Bio.SeqIO
import Bio.SeqUtils
import mendeleev
import unittest


class TestCore(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        records = Bio.SeqIO.parse('tests/fixtures/seq.fna', 'fasta')
        cls.chr_seq = next(records).seq

    def test_PolymerStrand(self):
        self.assertEqual(core.PolymerStrand.negative.value, -1)
        self.assertEqual(core.PolymerStrand.positive.value, 1)

    def test_RnaType(self):
        self.assertEqual(core.RnaType.mRna.value, 0)
        self.assertEqual(core.RnaType.rRna.value, 1)
        self.assertEqual(core.RnaType.sRna.value, 2)
        self.assertEqual(core.RnaType.tRna.value, 3)
        self.assertEqual(core.RnaType.mixed.value, 4)

    def test_GeneType(self):
        self.assertEqual(core.GeneType.mRna.value, 0)
        self.assertEqual(core.GeneType.rRna.value, 1)
        self.assertEqual(core.GeneType.sRna.value, 2)
        self.assertEqual(core.GeneType.tRna.value, 3)


class KnowledgeBaseTestCase(unittest.TestCase):
    def test_constructor(self):
        kb = core.KnowledgeBase()
        self.assertEqual(kb.cell, None)


class CellTestCase(unittest.TestCase):
    def test_constructor(self):
        cell = core.Cell()
        self.assertEqual(cell.knowledge_base, None)
        self.assertEqual(cell.species_types, [])
        self.assertEqual(list(cell.get_dna_species_types()), [])

        dna = [
            core.DnaSpeciesType(id='chr1', seq='AAA'),
            core.DnaSpeciesType(id='chr2', seq='CCC'),
        ]

        cell = core.Cell(species_types=dna)
        self.assertEqual(cell.species_types, dna)


class CompartmentTestCase(unittest.TestCase):
    def test_constructor(self):
        comp = core.Compartment(volume=2.)
        self.assertEqual(comp.volume, 2.)


class SpeciesTypeTestCase(unittest.TestCase):
    def test_constructor(self):
        with self.assertRaisesRegexp(TypeError, 'Can\'t instantiate abstract class'):
            core.SpeciesType()

        class ConcreteSpeciesType(core.SpeciesType):
            def get_structure(self): pass

            def get_empirical_formula(self): pass

            def get_charge(self): pass

            def get_mol_wt(self): pass

        species_type = ConcreteSpeciesType(concentration=2., half_life=3.)
        self.assertEqual(species_type.concentration, 2.)
        self.assertEqual(species_type.half_life, 3.)


class MetaboliteSpeciesTypeTestCase(unittest.TestCase):
    def test_constructor(self):
        met = core.MetaboliteSpeciesType(structure=(
            'InChI=1S'
            '/C10H14N5O7P'
            '/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(22-10)1-21-23(18,19)20'
            '/h2-4,6-7,10,16-17H,1H2,(H2,11,12,13)(H2,18,19,20)'
            '/p-2/t4-,6-,7-,10-'
            '/m1'
            '/s1'
        ))
        self.assertEqual(met.get_structure(), met.structure)
        self.assertEqual(met.get_empirical_formula(), chem.EmpiricalFormula('C10H12N5O7P'))
        self.assertEqual(met.get_charge(), -2)
        self.assertAlmostEqual(met.get_mol_wt(), 345.20530, places=4)


class PolymerSpeciesTypeTestCase(unittest.TestCase):
    def setUp(self):
        class ConcretePolymerSpeciesType(core.PolymerSpeciesType):
            def get_seq(self): return Bio.Seq.Seq('AAATGCCC', alphabet=Bio.Alphabet.DNAAlphabet())

            def get_empirical_formula(self): pass

            def get_charge(self): pass

            def get_mol_wt(self): pass
        self.ConcretePolymerSpeciesType = ConcretePolymerSpeciesType

    def test_constructor(self):
        species_type = self.ConcretePolymerSpeciesType()
        species_type.circular = True
        species_type.double_stranded = True

    def test_get_len(self):
        species_type = self.ConcretePolymerSpeciesType()
        self.assertEqual(species_type.get_len(), 8)

    def test_get_subseq_linear(self):
        species_type = self.ConcretePolymerSpeciesType()
        species_type.circular = False

        self.assertEqual(species_type.get_subseq(1, 2), 'AA')
        self.assertEqual(species_type.get_subseq(2, 4), 'AAT')

        with self.assertRaisesRegexp(ValueError, 'Start and end coordinates for linear polymers must be'):
            species_type.get_subseq(0, 2)

        with self.assertRaisesRegexp(ValueError, 'Start and end coordinates for linear polymers must be'):
            species_type.get_subseq(1, 9)

    def test_get_subseq_circular(self):
        species_type = self.ConcretePolymerSpeciesType()
        species_type.circular = True

        self.assertEqual(species_type.get_subseq(1, 2), 'AA')
        self.assertEqual(species_type.get_subseq(2, 4), 'AAT')

        self.assertEqual(species_type.get_subseq(3, 6), 'ATGC')
        self.assertEqual(species_type.get_subseq(3, 6, strand=core.PolymerStrand.positive), 'ATGC')
        self.assertEqual(species_type.get_subseq(3, 6, strand=core.PolymerStrand.negative), 'GCAT')

        self.assertEqual(species_type.get_subseq(0, 1), 'CA')
        self.assertEqual(species_type.get_subseq(-1, 1), 'CCA')
        self.assertEqual(species_type.get_subseq(-3, 1), 'GCCCA')

        self.assertEqual(species_type.get_subseq(6, 10), 'CCCAA')
        self.assertEqual(species_type.get_subseq(6, 10, strand=core.PolymerStrand.negative), 'TTGGG')

        self.assertEqual(species_type.get_subseq(6, 26), 'CCCAAATGCCCAAATGCCCAA')
        self.assertEqual(species_type.get_subseq(-2, 18), 'CCCAAATGCCCAAATGCCCAA')
        self.assertEqual(species_type.get_subseq(-10, 10), 'CCCAAATGCCCAAATGCCCAA')


class DnaSpeciesTypeTestCase(unittest.TestCase):

    def test(self):
        L = 8
        dna = core.DnaSpeciesType(seq=Bio.Seq.Seq('ACGTACGT', alphabet=Bio.Alphabet.DNAAlphabet()))

        dna.circular = False
        dna.double_stranded = False
        self.assertEqual(dna.get_empirical_formula(),
                         chem.EmpiricalFormula('C10H12N5O6P') * 2
                         + chem.EmpiricalFormula('C9H12N3O7P') * 2
                         + chem.EmpiricalFormula('C10H12N5O7P') * 2
                         + chem.EmpiricalFormula('C10H13N2O8P') * 2
                         - chem.EmpiricalFormula('OH') * (L - 1)
                         )
        self.assertEqual(dna.get_charge(), -L - 1)
        exp_mol_wt = \
            + Bio.SeqUtils.molecular_weight(dna.get_seq(),
                                            seq_type='DNA',
                                            circular=dna.circular,
                                            double_stranded=dna.double_stranded) \
            - 9 * mendeleev.element('H').atomic_weight
        self.assertAlmostEqual(dna.get_mol_wt(), exp_mol_wt, places=0)

        dna.circular = True
        dna.double_stranded = False
        self.assertEqual(dna.get_empirical_formula(),
                         chem.EmpiricalFormula('C10H12N5O6P') * 2
                         + chem.EmpiricalFormula('C9H12N3O7P') * 2
                         + chem.EmpiricalFormula('C10H12N5O7P') * 2
                         + chem.EmpiricalFormula('C10H13N2O8P') * 2
                         - chem.EmpiricalFormula('OH') * L
                         )
        self.assertEqual(dna.get_charge(), -L)
        exp_mol_wt = \
            + Bio.SeqUtils.molecular_weight(dna.get_seq(),
                                            seq_type='DNA',
                                            circular=dna.circular,
                                            double_stranded=dna.double_stranded) \
            - 8 * mendeleev.element('H').atomic_weight
        self.assertAlmostEqual(dna.get_mol_wt(), exp_mol_wt, places=0)

        dna.circular = False
        dna.double_stranded = True
        self.assertEqual(dna.get_empirical_formula(),
                         chem.EmpiricalFormula('C10H12N5O6P') * 2 * 2
                         + chem.EmpiricalFormula('C9H12N3O7P') * 2 * 2
                         + chem.EmpiricalFormula('C10H12N5O7P') * 2 * 2
                         + chem.EmpiricalFormula('C10H13N2O8P') * 2 * 2
                         - chem.EmpiricalFormula('OH') * (L - 1) * 2
                         )
        self.assertEqual(dna.get_charge(), 2 * (-L - 1))
        exp_mol_wt = \
            + Bio.SeqUtils.molecular_weight(dna.get_seq(),
                                            seq_type='DNA',
                                            circular=dna.circular,
                                            double_stranded=dna.double_stranded) \
            - 18 * mendeleev.element('H').atomic_weight
        self.assertAlmostEqual(dna.get_mol_wt(), exp_mol_wt, places=0)

        dna.circular = True
        dna.double_stranded = True
        self.assertEqual(dna.get_empirical_formula(),
                         chem.EmpiricalFormula('C10H12N5O6P') * 2 * 2
                         + chem.EmpiricalFormula('C9H12N3O7P') * 2 * 2
                         + chem.EmpiricalFormula('C10H12N5O7P') * 2 * 2
                         + chem.EmpiricalFormula('C10H13N2O8P') * 2 * 2
                         - chem.EmpiricalFormula('OH') * L * 2
                         )
        self.assertEqual(dna.get_charge(), 2 * -L)
        exp_mol_wt = \
            + Bio.SeqUtils.molecular_weight(dna.get_seq(),
                                            seq_type='DNA',
                                            circular=dna.circular,
                                            double_stranded=dna.double_stranded) \
            - 16 * mendeleev.element('H').atomic_weight
        self.assertAlmostEqual(dna.get_mol_wt(), exp_mol_wt, places=0)


class RnaSpeciesTypeTestCase(unittest.TestCase):
    def test_constructor(self):
        rna = core.RnaSpeciesType(id='rna_1', start=1, end=2, strand=core.PolymerStrand.positive)

        self.assertEqual(rna.id, 'rna_1')
        self.assertEqual(rna.start, 1)
        self.assertEqual(rna.end, 2)
        self.assertEqual(rna.strand, core.PolymerStrand.positive)

    def test_get_3_prime(self):
        rna = core.RnaSpeciesType(start=100, end=200, strand=core.PolymerStrand.positive)
        self.assertEqual(rna.get_3_prime(), 200)

        rna = core.RnaSpeciesType(start=100, end=200, strand=core.PolymerStrand.negative)
        self.assertEqual(rna.get_3_prime(), 100)

    def test_get_5_prime(self):
        rna = core.RnaSpeciesType(start=100, end=200, strand=core.PolymerStrand.positive)
        self.assertEqual(rna.get_5_prime(), 100)

        rna = core.RnaSpeciesType(start=100, end=200, strand=core.PolymerStrand.negative)
        self.assertEqual(rna.get_5_prime(), 200)

    def test_get_len(self):
        rna = core.RnaSpeciesType(start=100, end=200)
        self.assertEqual(rna.get_len(), 101)

    def test_get_seq(self):
        dna = core.DnaSpeciesType(seq=Bio.Seq.Seq('AAATGCCC', alphabet=Bio.Alphabet.DNAAlphabet()))

        rna = core.RnaSpeciesType(dna=dna, start=3, end=6)
        self.assertEqual(rna.get_seq(), 'AUGC')

        rna = core.RnaSpeciesType(dna=dna, start=3, end=6, strand=core.PolymerStrand.positive)
        self.assertEqual(rna.get_seq(), 'AUGC')

        rna = core.RnaSpeciesType(dna=dna, start=3, end=6, strand=core.PolymerStrand.negative)
        self.assertEqual(rna.get_seq(), 'GCAU')

    def test_get_empirical_formula(self):
        dna = core.DnaSpeciesType()

        rna = core.RnaSpeciesType(dna=dna, start=1, end=1)
        dna.seq = Bio.Seq.Seq('A', alphabet=Bio.Alphabet.DNAAlphabet())
        self.assertEqual(rna.get_empirical_formula(), chem.EmpiricalFormula('C10H12N5O7P'))

        rna = core.RnaSpeciesType(dna=dna, start=1, end=1)
        dna.seq = Bio.Seq.Seq('C', alphabet=Bio.Alphabet.DNAAlphabet())
        self.assertEqual(rna.get_empirical_formula(), chem.EmpiricalFormula('C9H12N3O8P'))

        rna = core.RnaSpeciesType(dna=dna, start=1, end=1)
        dna.seq = Bio.Seq.Seq('G', alphabet=Bio.Alphabet.DNAAlphabet())
        self.assertEqual(rna.get_empirical_formula(), chem.EmpiricalFormula('C10H12N5O8P'))

        rna = core.RnaSpeciesType(dna=dna, start=1, end=1)
        dna.seq = Bio.Seq.Seq('T', alphabet=Bio.Alphabet.DNAAlphabet())
        self.assertEqual(rna.get_empirical_formula(), chem.EmpiricalFormula('C9H11N2O9P'))

        rna = core.RnaSpeciesType(dna=dna, start=1, end=2)
        dna.seq = Bio.Seq.Seq('AA', alphabet=Bio.Alphabet.DNAAlphabet())
        self.assertEqual(rna.get_empirical_formula(), chem.EmpiricalFormula('C20H23N10O13P2'))

    def test_get_charge(self):
        rna = core.RnaSpeciesType(start=1, end=1)
        self.assertEqual(rna.get_charge(), -2)

        rna = core.RnaSpeciesType(start=1, end=2)
        self.assertEqual(rna.get_charge(), -3)

    def test_get_mol_wt(self):
        dna = core.DnaSpeciesType(seq=Bio.Seq.Seq('AACCGGTT'))

        rna = core.RnaSpeciesType(dna=dna, start=1, end=1)
        exp_mol_wt = \
            + Bio.SeqUtils.molecular_weight(rna.get_seq()) \
            - (rna.get_len() + 1) * mendeleev.element('H').atomic_weight
        self.assertAlmostEqual(rna.get_mol_wt(), exp_mol_wt, places=1)

        rna = core.RnaSpeciesType(dna=dna, start=3, end=3)
        exp_mol_wt = \
            + Bio.SeqUtils.molecular_weight(rna.get_seq()) \
            - (rna.get_len() + 1) * mendeleev.element('H').atomic_weight
        self.assertAlmostEqual(rna.get_mol_wt(), exp_mol_wt, places=1)

        rna = core.RnaSpeciesType(dna=dna, start=5, end=5)
        exp_mol_wt = \
            + Bio.SeqUtils.molecular_weight(rna.get_seq()) \
            - (rna.get_len() + 1) * mendeleev.element('H').atomic_weight
        self.assertAlmostEqual(rna.get_mol_wt(), exp_mol_wt, places=1)

        rna = core.RnaSpeciesType(dna=dna, start=7, end=7)
        exp_mol_wt = \
            + Bio.SeqUtils.molecular_weight(rna.get_seq()) \
            - (rna.get_len() + 1) * mendeleev.element('H').atomic_weight
        self.assertAlmostEqual(rna.get_mol_wt(), exp_mol_wt, places=1)

        rna = core.RnaSpeciesType(dna=dna, start=1, end=2)
        exp_mol_wt = \
            + Bio.SeqUtils.molecular_weight(rna.get_seq()) \
            - (rna.get_len() + 1) * mendeleev.element('H').atomic_weight
        self.assertAlmostEqual(rna.get_mol_wt(), exp_mol_wt, places=1)


class ProteinSpeciesTypeTestCase(unittest.TestCase):
    def test_constructor(self):
        protein = core.ProteinSpeciesType(id='prot1', name='protein_1', concentration = 1., half_life =5.)
        self.assertEqual(protein.id, 'prot1')
        self.assertEqual(protein.name, 'protein_1')
        self.assertEqual(protein.reaction_participants, [])
        self.assertEqual(protein.concentration, 1.)
        self.assertEqual(protein.half_life, 5.)
        self.assertEqual(protein.circular, False)
        self.assertEqual(protein.cell, None)
        self.assertEqual(protein.orfs, [])
        self.assertEqual(protein.loci, [])

    def test_get_seq(self):
        records = Bio.SeqIO.parse('tests/fixtures/seq.fna', 'fasta')
        dna_seq = next(records).seq
        dna = core.DnaSpeciesType(seq=dna_seq)
        cell = dna.cell = core.Cell()
        cell.knowledge_base = core.KnowledgeBase(translation_table=4)

        # MPN001
        rna = core.RnaSpeciesType(dna=dna, start=692, end=1834, strand=core.PolymerStrand.positive)
        self.assertEqual(rna.get_seq()[0:10], 'AUGAAAGUUU')
        self.assertEqual(rna.get_seq()[-10:], 'UUCCAAGUAA')

        orf = core.OpenReadingFrameLocus(polymer=rna, start=1, end=rna.get_len())
        self.assertEqual(orf.get_seq()[0:10], 'AUGAAAGUUU')
        self.assertEqual(orf.get_seq()[-10:], 'UUCCAAGUAA')

        prot = core.ProteinSpeciesType(orfs=[orf])
        self.assertEqual(prot.get_seq()[0:10], 'MKVLINKNEL')

        # MPN011
        rna = core.RnaSpeciesType(dna=dna, start=12838, end=13533, strand=core.PolymerStrand.negative)
        self.assertEqual(rna.get_seq()[0:10], 'AUGAAAUUUA')
        self.assertEqual(rna.get_seq()[-10:], 'AAUUGAGUAA')

        orf = core.OpenReadingFrameLocus(polymer=rna, start=1, end=rna.get_len())
        self.assertEqual(orf.get_seq()[0:10], 'AUGAAAUUUA')
        self.assertEqual(orf.get_seq()[-10:], 'AAUUGAGUAA')

        prot = core.ProteinSpeciesType(orfs=[orf])
        self.assertEqual(prot.get_seq()[0:10], 'MKFKFLLTPL')

    @unittest.skip('todo')
    def test_get_empirical_formula(self):
        # Test is based on Collagen Type IV a3 (https://pubchem.ncbi.nlm.nih.gov/compound/44511378)
        dna1 = core.DnaSpeciesType(seq=Bio.Seq.Seq('TGTAATTATTATTCTAATTCTTATTCTTTTTGGTTAGCTTCTTTAAATCCTGAACGT', alphabet=Bio.Alphabet.DNAAlphabet()))
        rna1 = core.RnaSpeciesType(dna=dna1)
        prot1 = core.ProteinSpeciesType(rna=rna1)
        self.assertEqual(prot1.get_empirical_formula(),chem.EmpiricalFormula('C105H144N26O32S'))

        # Test is based on Histone 7 (https://pubchem.ncbi.nlm.nih.gov/compound/22461943)
        dna2 = core.DnaSpeciesType(seq=Bio.Seq.Seq('TCGCCGCAAAGCTTCTGGTCCT', alphabet=Bio.Alphabet.DNAAlphabet()))
        rna2 = core.RnaSpeciesType(dna=dna2)
        prot2 = core.ProteinSpeciesType(rna=rna2)
        self.assertEqual(prot2.get_empirical_formula(),chem.EmpiricalFormula('C31H58N14O9'))

    @unittest.skip('todo')
    def test_get_mol_wt(self):
        # Test is based on Collagen Type IV a3 (https://pubchem.ncbi.nlm.nih.gov/compound/44511378)
        dna1 = core.DnaSpeciesType(seq=Bio.Seq.Seq('TGTAATTATTATTCTAATTCTTATTCTTTTTGGTTAGCTTCTTTAAATCCTGAACGT', alphabet=Bio.Alphabet.DNAAlphabet()))
        rna1 = core.RnaSpeciesType(dna=dna1)
        prot1 = core.ProteinSpeciesType(rna=rna1)
        self.assertAlmostEqual(prot1.get_mol_wt(),2314.517) # Double check units match, ref from DB is in g/mol

        # Test is based on Histone 7 (https://pubchem.ncbi.nlm.nih.gov/compound/22461943)
        dna2 = core.DnaSpeciesType(seq=Bio.Seq.Seq('TCGCCGCAAAGCTTCTGGTCCT', alphabet=Bio.Alphabet.DNAAlphabet()))
        rna2 = core.RnaSpeciesType(dna=dna2)
        prot2 = core.ProteinSpeciesType(rna=rna2)
        self.assertAlmostEqual(prot2.get_mol_wt(),770.894)

    @unittest.skip('todo')
    def test_get_charge(self):
        pass


class PolymerLocusTestCase(unittest.TestCase):
    def test_constructor(self):
        dna = core.DnaSpeciesType()
        locus = core.PolymerLocus(polymer=dna, start=10, end=20, strand=core.PolymerStrand.positive)
        self.assertEqual(locus.start, 10)
        self.assertEqual(locus.end, 20)


class GeneLocusTestCase(unittest.TestCase):
    def test_constructor(self):
        gene = core.GeneLocus(id='gene_1', name='gene_1', symbol='gene_1',
                              start=1, end=2, strand=core.PolymerStrand.positive,
                              type=core.GeneType.mRna)
        rna = core.RnaSpeciesType(id='rna_1', name='rna_1', genes=[gene],
                                  type=core.RnaType.mRna)

        self.assertEqual(gene.id, 'gene_1')
        self.assertEqual(gene.name, 'gene_1')
        self.assertEqual(gene.rnas, [rna])
        self.assertEqual(gene.type.name, 'mRna')
        self.assertEqual(gene.symbol, 'gene_1')
    def test_transcription_unit_constructor(self):
        locus = core.PromoterLocus(id='tu_1', start=3, end=4, strand=core.PolymerStrand.positive)

        self.assertEqual(locus.start, 3)
        self.assertEqual(locus.end, 4)
        self.assertEqual(locus.strand, core.PolymerStrand.positive)

    def test_get_len(self):
        locus = core.PromoterLocus(start=-40, end=-50)
        self.assertEqual(locus.get_len(), 11)

    def test_get_seq(self):
        dna = core.DnaSpeciesType(seq=Bio.Seq.Seq('ACGTACGTACGTACGT', alphabet=Bio.Alphabet.DNAAlphabet()),
                                  circular=True, double_stranded=True)
        rna = core.RnaSpeciesType(dna=dna, start=7, end=8, strand=core.PolymerStrand.positive)
        locus = core.PromoterLocus(polymer=dna, rnas=[rna], pribnow_start=-2, pribnow_end=-4)
        self.assertEqual(locus.get_pribnow_seq(), 'ATG')

        dna = core.DnaSpeciesType(seq=Bio.Seq.Seq('ACGTACGTACGTACGT', alphabet=Bio.Alphabet.DNAAlphabet()),
                                  circular=True, double_stranded=True)
        rna = core.RnaSpeciesType(dna=dna, start=7, end=8, strand=core.PolymerStrand.negative)
        locus = core.PromoterLocus(polymer=dna, rnas=[rna], pribnow_start=-2, pribnow_end=-4)
        self.assertEqual(locus.get_pribnow_seq(), 'GCA')

        records = Bio.SeqIO.parse('tests/fixtures/seq.fna', 'fasta')
        dna_seq = next(records).seq
        dna = core.DnaSpeciesType(seq=dna_seq, circular=True, double_stranded=True)

        rna = core.RnaSpeciesType(dna=dna, start=652, end=1900, strand=core.PolymerStrand.positive)
        locus = core.PromoterLocus(polymer=dna, rnas=[rna], pribnow_start=-35, pribnow_end=-40)
        self.assertEqual(locus.get_pribnow_seq(), 'TAAAAC')

        rna = core.RnaSpeciesType(dna=dna, start=1252, end=1377, strand=core.PolymerStrand.negative)
        locus = core.PromoterLocus(polymer=dna, rnas=[rna], pribnow_start=-125, pribnow_end=-130)
        self.assertEqual(locus.get_pribnow_seq(), 'TAAGTT')


class OpenReadingFrameLocusTestCase(unittest.TestCase):
    def test(self):
        orf = core.OpenReadingFrameLocus()
        prot = core.ProteinSpeciesType(orfs=[orf])


@unittest.skip('todo')
class ReactionTestCase(unittest.TestCase):
    pass
