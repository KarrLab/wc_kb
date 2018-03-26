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
        self.assertEqual(list(cell.get_species_types(core.DnaSpeciesType)), [])

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
            def get_empirical_formula(self): pass
            def get_charge(self): pass
            def get_mol_wt(self): pass

        species_type = ConcreteSpeciesType(id='species1', name = 'species1', concentration=2., half_life=3.)
        self.assertEqual(species_type.id, 'species1')
        self.assertEqual(species_type.name, 'species1')
        self.assertEqual(species_type.concentration, 2.)
        self.assertEqual(species_type.half_life, 3.)


class PolymerSpeciesTypeTestCase(unittest.TestCase):
    def test_constructor(self):
        with self.assertRaisesRegexp(TypeError, 'Can\'t instantiate abstract class'):
            core.PolymerSpeciesType()

        class ConcretePolymerSpeciesType(core.PolymerSpeciesType):
            def get_empirical_formula(self): pass
            def get_charge(self): pass
            def get_mol_wt(self): pass
            def get_seq(self): return Bio.Seq.Seq('AAATGCCC', alphabet=Bio.Alphabet.DNAAlphabet())

            # Methods from PolymerSpeciesType
            #def get_subseq(self): pass
            #def get_len(self): pass

        pst1 = ConcretePolymerSpeciesType(id='pst1', name = 'pst1', concentration=1, half_life=2)

        #Test constructor
        self.assertEqual(pst1.id, 'pst1')
        self.assertEqual(pst1.name, 'pst1')
        self.assertEqual(pst1.concentration, 1)
        self.assertEqual(pst1.half_life, 2)

        # Test methods: linear, single stranded case
        pst1.circular = False
        pst1.double_stranded = False

        self.assertEqual(pst1.get_len(), 8)
        self.assertEqual(pst1.get_subseq(1, 3), 'AAA')
        self.assertEqual(pst1.get_subseq(2, 4), 'AAT')

        # Test methods: circular, single stranded case
        pst1.circular = True
        pst1.double_stranded = False

        self.assertEqual(pst1.get_subseq(2, 4), 'AAT')
        self.assertEqual(pst1.get_subseq(0, 1), 'CA')
        self.assertEqual(pst1.get_subseq(-3, 1), 'GCCCA')
        self.assertEqual(pst1.get_subseq(6, 10), 'CCCAA')
        self.assertEqual(pst1.get_subseq(-10, 10), 'CCCAAATGCCCAAATGCCCAA')

        # Test methods: circular, double stranded case
        pst1.circular = True
        pst1.double_stranded = True

        self.assertEqual(pst1.get_subseq(3,  6, strand=core.PolymerStrand.positive), 'ATGC')
        self.assertEqual(pst1.get_subseq(3,  6, strand=core.PolymerStrand.negative), 'GCAT')
        self.assertEqual(pst1.get_subseq(6, 26, strand=core.PolymerStrand.positive), 'CCCAAATGCCCAAATGCCCAA')
        self.assertEqual(pst1.get_subseq(6, 26, strand=core.PolymerStrand.negative), 'TTGGGCATTTGGGCATTTGGG')


class DnaSpeciesTypeTestCase(unittest.TestCase):
    def test(self):
        dna = core.DnaSpeciesType(id = 'dna1', name = 'dna1', seq=Bio.Seq.Seq('ACGTACGT', alphabet=Bio.Alphabet.DNAAlphabet()),
              circular = False, double_stranded = False)

        self.assertEqual(dna.id, 'dna1')
        self.assertEqual(dna.name, 'dna1')
        self.assertEqual(dna.circular, False)
        self.assertEqual(dna.double_stranded, False)

        L = dna.get_len()
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

        # Make DNA circular, single stranded
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

        # Make DNA linear, double stranded
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

        # Make DNA circular, double stranded
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
        dna1 = core.DnaSpeciesType(id='dna1', seq=Bio.Seq.Seq('ACGTACGTACGTACG', alphabet=Bio.Alphabet.DNAAlphabet()))
        tu1 = core.TranscriptionUnitLocus(id='tu1', polymer=dna1, start=1, end=15)
        rna1 = core.RnaSpeciesType(id='rna1', name='rna1', transcription_unit=[tu1], type=1, concentration=1, half_life=2)

        self.assertEqual(rna1.id, 'rna1')
        self.assertEqual(rna1.name, 'rna1')
        self.assertEqual(rna1.transcription_unit, [tu1])
        self.assertEqual(rna1.type, 1)
        self.assertEqual(rna1.concentration, 1)
        self.assertEqual(rna1.half_life, 2)

    def test_get_empirical_formula(self):
        dna1 = core.DnaSpeciesType(id='dna1', seq=Bio.Seq.Seq('A', alphabet=Bio.Alphabet.DNAAlphabet()))
        tu1 = core.TranscriptionUnitLocus(id='tu1', polymer=dna1, start=1, end=1)
        rna1 = core.RnaSpeciesType(id='rna1', name='rna1', transcription_unit=[tu1])
        self.assertEqual(rna1.get_empirical_formula(), chem.EmpiricalFormula('C10H12N5O7P'))

        dna1 = core.DnaSpeciesType(id='dna1', seq=Bio.Seq.Seq('C', alphabet=Bio.Alphabet.DNAAlphabet()))
        tu1 = core.TranscriptionUnitLocus(id='tu1', polymer=dna1, start=1, end=1)
        rna1 = core.RnaSpeciesType(id='rna1', name='rna1', transcription_unit=[tu1])
        self.assertEqual(rna1.get_empirical_formula(), chem.EmpiricalFormula('C9H12N3O8P'))

        dna1 = core.DnaSpeciesType(id='dna1', seq=Bio.Seq.Seq('G', alphabet=Bio.Alphabet.DNAAlphabet()))
        tu1 = core.TranscriptionUnitLocus(id='tu1', polymer=dna1, start=1, end=1)
        rna1 = core.RnaSpeciesType(id='rna1', name='rna1', transcription_unit=[tu1])
        self.assertEqual(rna1.get_empirical_formula(), chem.EmpiricalFormula('C10H12N5O8P'))

        dna1 = core.DnaSpeciesType(id='dna1', seq=Bio.Seq.Seq('T', alphabet=Bio.Alphabet.DNAAlphabet()))
        tu1 = core.TranscriptionUnitLocus(id='tu1', polymer=dna1, start=1, end=1)
        rna1 = core.RnaSpeciesType(id='rna1', name='rna1', transcription_unit=[tu1])
        self.assertEqual(rna1.get_empirical_formula(), chem.EmpiricalFormula('C9H11N2O9P'))

        dna1 = core.DnaSpeciesType(id='dna1', seq=Bio.Seq.Seq('AAAA', alphabet=Bio.Alphabet.DNAAlphabet()))
        tu1 = core.TranscriptionUnitLocus(id='tu1', polymer=dna1, start=1, end=2)
        rna1 = core.RnaSpeciesType(id='rna1', name='rna1', transcription_unit=[tu1])
        self.assertEqual(rna1.get_empirical_formula(), chem.EmpiricalFormula('C20H23N10O13P2'))

    def test_get_charge(self):
        dna1 = core.DnaSpeciesType(id='dna1', seq=Bio.Seq.Seq('AAAA', alphabet=Bio.Alphabet.DNAAlphabet()))
        tu1 = core.TranscriptionUnitLocus(id='tu1', polymer=dna1, start=1, end=1)
        rna1 = core.RnaSpeciesType(id='rna1', name='rna1', transcription_unit=[tu1])
        self.assertEqual(rna1.get_charge(), -2)

        dna1 = core.DnaSpeciesType(id='dna1', seq=Bio.Seq.Seq('AAAA', alphabet=Bio.Alphabet.DNAAlphabet()))
        tu1 = core.TranscriptionUnitLocus(id='tu1', polymer=dna1, start=1, end=2)
        rna1 = core.RnaSpeciesType(id='rna1', name='rna1', transcription_unit=[tu1])
        self.assertEqual(rna1.get_charge(), -3)

    def test_get_mol_wt(self):
        dna1 = core.DnaSpeciesType(id='dna1', seq=Bio.Seq.Seq('AACCGGTT', alphabet=Bio.Alphabet.DNAAlphabet()))
        tu1 = core.TranscriptionUnitLocus(id='tu1', polymer=dna1, start=1, end=1)
        rna1 = core.RnaSpeciesType(id='rna1', name='rna1', transcription_unit=[tu1])
        exp_mol_wt = \
            + Bio.SeqUtils.molecular_weight(rna1.get_seq()) \
            - (rna1.get_len() + 1) * mendeleev.element('H').atomic_weight
        self.assertAlmostEqual(rna1.get_mol_wt(), exp_mol_wt, places=1)

        dna1 = core.DnaSpeciesType(id='dna1', seq=Bio.Seq.Seq('AACCGGTT', alphabet=Bio.Alphabet.DNAAlphabet()))
        tu1 = core.TranscriptionUnitLocus(id='tu1', polymer=dna1, start=3, end=3)
        rna1 = core.RnaSpeciesType(id='rna1', name='rna1', transcription_unit=[tu1])
        exp_mol_wt = \
            + Bio.SeqUtils.molecular_weight(rna1.get_seq()) \
            - (rna1.get_len() + 1) * mendeleev.element('H').atomic_weight
        self.assertAlmostEqual(rna1.get_mol_wt(), exp_mol_wt, places=1)

        dna1 = core.DnaSpeciesType(id='dna1', seq=Bio.Seq.Seq('AACCGGTT', alphabet=Bio.Alphabet.DNAAlphabet()))
        tu1 = core.TranscriptionUnitLocus(id='tu1', polymer=dna1, start=5, end=5)
        rna1 = core.RnaSpeciesType(id='rna1', name='rna1', transcription_unit=[tu1])
        exp_mol_wt = \
            + Bio.SeqUtils.molecular_weight(rna1.get_seq()) \
            - (rna1.get_len() + 1) * mendeleev.element('H').atomic_weight
        self.assertAlmostEqual(rna1.get_mol_wt(), exp_mol_wt, places=1)


class ProteinSpeciesTypeTestCase(unittest.TestCase):
    def test_constructor(self):
        protein = core.ProteinSpeciesType(id='prot1', name='prot1', concentration=1, half_life=2)
        # attribute_order = ('id', 'cell', 'name', 'gene', 'rna', 'concentration', 'half_life', 'comments')

        self.assertEqual(protein.id, 'prot1')
        self.assertEqual(protein.name, 'prot1')
        self.assertEqual(protein.concentration, 1)
        self.assertEqual(protein.half_life, 2)
        self.assertEqual(protein.cell, None)

    @unittest.skip('Work in progress')
    def test_get_seq(self):
        records = Bio.SeqIO.parse('tests/fixtures/seq.fna', 'fasta')
        dna_seq = next(records).seq
        dna = core.DnaSpeciesType(seq=dna_seq)
        cell.knowledge_base = core.KnowledgeBase(translation_table=1)
        cell = dna.cell = core.Cell()

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

    @unittest.skip('Work in progress')
    def test_get_empirical_formula(self):
        # Test is based on Collagen Type IV a3 (https://pubchem.ncbi.nlm.nih.gov/compound/44511378)
        dna1 = core.DnaSpeciesType(seq=Bio.Seq.Seq('TGTAATTATTATTCTAATTCTTATTCTTTTTGGTTAGCTTCTTTAAATCCTGAACGT', alphabet=Bio.Alphabet.DNAAlphabet()))
        cell = dna1.cell = core.Cell()
        cell.knowledge_base = core.KnowledgeBase(translation_table=1)
        rna1 = core.RnaSpeciesType(dna=dna1, start=1, end=dna1.get_len(), strand=core.PolymerStrand.positive)
        orf1 = core.OpenReadingFrameLocus(polymer=rna1, start=1, end=rna1.get_len())
        prot1 = core.ProteinSpeciesType(orfs=[orf1])
        self.assertEqual(prot1.get_empirical_formula(), chem.EmpiricalFormula('C105H144N26O32S'))

        # Test is based on Tuftsin (hhttps://pubchem.ncbi.nlm.nih.gov/compounds/156080)
        dna2 = core.DnaSpeciesType(seq=Bio.Seq.Seq('ACTAAACCTCGT', alphabet=Bio.Alphabet.DNAAlphabet()))
        cell = dna2.cell = core.Cell()
        cell.knowledge_base = core.KnowledgeBase(translation_table=1)
        rna2 = core.RnaSpeciesType(dna=dna2, start=1, end=dna2.get_len(), strand=core.PolymerStrand.positive)
        orf2 = core.OpenReadingFrameLocus(polymer=rna2, start=1, end=rna2.get_len())
        prot2 = core.ProteinSpeciesType(orfs=[orf2])
        self.assertEqual(prot2.get_empirical_formula(), chem.EmpiricalFormula('C21H40N8O6'))

    @unittest.skip('Work in progress')
    def test_get_mol_wt(self):
        # Test is based on Collagen Type IV a3 (https://pubchem.ncbi.nlm.nih.gov/compound/44511378)
        dna1 = core.DnaSpeciesType(seq=Bio.Seq.Seq(
            'TGTAATTATTATTCTAATTCTTATTCTTTTTGGTTAGCTTCTTTAAATCCTGAACGT', alphabet=Bio.Alphabet.DNAAlphabet()))
        cell = dna1.cell = core.Cell()
        cell.knowledge_base = core.KnowledgeBase(translation_table=1)
        rna1 = core.RnaSpeciesType(dna=dna1, start=1, end=dna1.get_len(), strand=core.PolymerStrand.positive)
        orf1 = core.OpenReadingFrameLocus(polymer=rna1, start=1, end=rna1.get_len())
        prot1 = core.ProteinSpeciesType(orfs=[orf1])
        self.assertAlmostEqual(prot1.get_mol_wt(), 2314.517)

        # Test is based on Tuftsin (hhttps://pubchem.ncbi.nlm.nih.gov/compounds/156080)
        dna2 = core.DnaSpeciesType(seq=Bio.Seq.Seq('ACTAAACCTCGT', alphabet=Bio.Alphabet.DNAAlphabet()))
        cell = dna2.cell = core.Cell()
        cell.knowledge_base = core.KnowledgeBase(translation_table=1)
        rna2 = core.RnaSpeciesType(dna=dna2, start=1, end=dna2.get_len(), strand=core.PolymerStrand.positive)
        orf2 = core.OpenReadingFrameLocus(polymer=rna2, start=1, end=rna2.get_len())
        prot2 = core.ProteinSpeciesType(orfs=[orf2])
        self.assertAlmostEqual(prot2.get_mol_wt(), 500.601)

    @unittest.skip('Work in progress')
    def test_get_charge(self):
        # Test is based on Collagen Type IV a3 (https://pubchem.ncbi.nlm.nih.gov/compound/44511378)
        dna1 = core.DnaSpeciesType(seq=Bio.Seq.Seq(
            'TGTAATTATTATTCTAATTCTTATTCTTTTTGGTTAGCTTCTTTAAATCCTGAACGT', alphabet=Bio.Alphabet.DNAAlphabet()))
        cell = dna1.cell = core.Cell()
        cell.knowledge_base = core.KnowledgeBase(translation_table=1)
        rna1 = core.RnaSpeciesType(dna=dna1, start=1, end=dna1.get_len(), strand=core.PolymerStrand.positive)
        orf1 = core.OpenReadingFrameLocus(polymer=rna1, start=1, end=rna1.get_len())
        prot1 = core.ProteinSpeciesType(orfs=[orf1])
        self.assertEqual(prot1.get_charge(), 0)

        # Test is based on Tuftsin (hhttps://pubchem.ncbi.nlm.nih.gov/compounds/156080)
        dna2 = core.DnaSpeciesType(seq=Bio.Seq.Seq('ACTAAACCTCGT', alphabet=Bio.Alphabet.DNAAlphabet()))
        cell = dna2.cell = core.Cell()
        cell.knowledge_base = core.KnowledgeBase(translation_table=1)
        rna2 = core.RnaSpeciesType(dna=dna2, start=1, end=dna2.get_len(), strand=core.PolymerStrand.positive)
        orf2 = core.OpenReadingFrameLocus(polymer=rna2, start=1, end=rna2.get_len())
        prot2 = core.ProteinSpeciesType(orfs=[orf2])
        self.assertEqual(prot2.get_charge(), 2)


class PolymerLocusTestCase(unittest.TestCase):
    def test_constructor(self):
        cell1 = core.Cell()
        dna1 = core.DnaSpeciesType(id = 'dna1', seq=Bio.Seq.Seq('ACGTACGTACGTACG', alphabet=Bio.Alphabet.DNAAlphabet()),
               circular = False, double_stranded = False)

        locus1 = core.PolymerLocus(id='locus1', cell=cell1, name='locus1', polymer=dna1, strand=core.PolymerStrand.positive, start=1, end=15)

        # test constructor
        self.assertEqual(locus1.id, 'locus1')
        self.assertEqual(locus1.cell, cell1)
        self.assertEqual(locus1.name, 'locus1')
        self.assertEqual(locus1.polymer, dna1)
        self.assertEqual(locus1.strand, core.PolymerStrand.positive)
        self.assertEqual(locus1.start, 1)
        self.assertEqual(locus1.end, 15)

        # test methods
        self.assertEqual(locus1.get_seq(),'ACGTACGTACGTACG')
        self.assertEqual(locus1.get_len(),15)

        # flip strand; test methods
        rev_comp_seq = locus1.get_seq().reverse_complement()
        locus1.strand = core.PolymerStrand.negative
        self.assertEqual(locus1.get_seq(), rev_comp_seq)
        self.assertEqual(locus1.get_len(),15)


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


class PromoterLocusTestCase(unittest.TestCase):
    def test_constructor(self):
        promoter = core.PromoterLocus(id='promoter1', name='promoter1', pribnow_start=1, pribnow_end=2)
        self.assertEqual(promoter.id, 'promoter1')
        self.assertEqual(promoter.name, 'promoter1')
        self.assertEqual(promoter.pribnow_start, 1)
        self.assertEqual(promoter.pribnow_end, 2)


class GeneLocusTestCase(unittest.TestCase):
    def test(self):
        gene = core.GeneLocus(id='gene1', name='gene1', symbol='gene_1', start=1, end=2)
        self.assertEqual(gene.id, 'gene1')
        self.assertEqual(gene.name, 'gene1')
        self.assertEqual(gene.symbol, 'gene_1')
        self.assertEqual(gene.start, 1)
        self.assertEqual(gene.end, 2)


class TranscriptionUnitLocusTestCase(unittest.TestCase):
    def test(self):
        dna1 = core.DnaSpeciesType(id = 'dna1', seq=Bio.Seq.Seq('ACGTACGTACGTACG', alphabet=Bio.Alphabet.DNAAlphabet()),
               circular = False, double_stranded = False)

        tu1 = core.TranscriptionUnitLocus(id='tu1', name='tu1', polymer = dna1, strand=core.PolymerStrand.positive, start=1, end=15)

        # test constructor
        self.assertEqual(tu1.id, 'tu1')
        self.assertEqual(tu1.name, 'tu1')
        self.assertEqual(tu1.polymer, dna1)
        self.assertEqual(tu1.strand, core.PolymerStrand.positive)
        self.assertEqual(tu1.start, 1)
        self.assertEqual(tu1.end, 15)

        # test methods
        self.assertEqual(tu1.get_3_prime(),15)
        self.assertEqual(tu1.get_5_prime(),1)

        # flip strand; test methods
        rev_comp_seq = tu1.get_seq().reverse_complement()
        tu1.strand = core.PolymerStrand.negative
        self.assertEqual(tu1.get_3_prime(),1)
        self.assertEqual(tu1.get_5_prime(),15)


class ReactionParticipantTestCase(unittest.TestCase):
    def test_constructor(self):
        cell1 = core.Cell()
        compartment1 = core.Compartment(cell = cell1)
        species1 = core.MetaboliteSpeciesType(id ='1')
        species2 = core.MetaboliteSpeciesType(id ='2')

        participant1 = core.ReactionParticipant(species_type=[species1, species2], compartment=[compartment1], coefficient=5)

        self.assertEqual(participant1.species_type, [species1, species2])
        self.assertEqual(participant1.compartment, [compartment1])
        self.assertEqual(participant1.coefficient, 5)


class ReactionTestCase(unittest.TestCase):
        def test_constructor(self):
            cell1 = core.Cell()
            compartment1 = core.Compartment(cell = cell1)
            species1 = core.MetaboliteSpeciesType(id ='1')
            species2 = core.MetaboliteSpeciesType(id ='2')
            participant1 = core.ReactionParticipant(species_type = [species1], compartment = [compartment1], coefficient = 1)
            participant2 = core.ReactionParticipant(species_type = [species2], compartment = [compartment1], coefficient = 1)

            reaction1 = core.Reaction(
                                id ='reaction1',
                                name = 'test_reaction',
                                cell = cell1,
                                participants =[participant1, participant2],
                                k_m = 0.1,
                                v_max = 0.5,
                                reversible=0)

            self.assertEqual(reaction1.id,'reaction1')
            self.assertEqual(reaction1.name,'test_reaction')
            self.assertEqual(reaction1.cell,cell1)
            self.assertEqual(reaction1.participants, [participant1, participant2])
            self.assertEqual(reaction1.k_m, 0.1)
            self.assertEqual(reaction1.v_max, 0.5)
            self.assertEqual(reaction1.reversible, 0)

        def test_constructor(self):
            cell1 = core.Cell()
            compartment1 = core.Compartment(cell=cell1)
            species1 = core.MetaboliteSpeciesType(id='1')
            species2 = core.MetaboliteSpeciesType(id='2')
            participant1 = core.ReactionParticipant(species_type=[species1], compartment=[compartment1], coefficient=1)
            participant2 = core.ReactionParticipant(species_type=[species2], compartment=[compartment1], coefficient=1)

            reaction1 = core.Reaction(
                id='reaction1',
                name='test_reaction',
                cell=cell1,
                participants=[participant1, participant2],
                k_m=0.1,
                v_max=0.5,
                reversible=0)

            self.assertEqual(reaction1.id, 'reaction1')
            self.assertEqual(reaction1.name, 'test_reaction')
            self.assertEqual(reaction1.cell, cell1)
            self.assertEqual(reaction1.participants, [participant1, participant2])
            self.assertEqual(reaction1.k_m, 0.1)
            self.assertEqual(reaction1.v_max, 0.5)
            self.assertEqual(reaction1.reversible, 0)
