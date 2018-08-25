""" Tests of the knowledge base

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Author: Bilal Shaikh <bilal.shaikh@columbia.edu>
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
        self.assertEqual(cell.compartments, [])
        self.assertEqual(cell.reactions, [])
        self.assertEqual(cell.loci, [])

        self.assertEqual(cell.species_types.get(
            __type=core.DnaSpeciesType), [])
        self.assertEqual(cell.loci.get(__type=core.PromoterLocus), [])


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

        species_type = ConcreteSpeciesType(
            id='species1', name='species1', concentration=2., half_life=3.)
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

            def get_seq(self): return Bio.Seq.Seq(
                'AAATGCCC', alphabet=Bio.Alphabet.DNAAlphabet())

        pst1 = ConcretePolymerSpeciesType(
            id='pst1', name='pst1', concentration=1, half_life=2)

        # Test constructor
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

        with self.assertRaisesRegexp(ValueError, 'Start and end coordinates'):
            self.assertEqual(pst1.get_subseq(-1, 3), 'AAA')

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

        self.assertEqual(pst1.get_subseq(
            3,  6, strand=core.PolymerStrand.positive), 'ATGC')
        self.assertEqual(pst1.get_subseq(
            3,  6, strand=core.PolymerStrand.negative), 'GCAT')
        self.assertEqual(pst1.get_subseq(
            6, 26, strand=core.PolymerStrand.positive), 'CCCAAATGCCCAAATGCCCAA')
        self.assertEqual(pst1.get_subseq(
            6, 26, strand=core.PolymerStrand.negative), 'TTGGGCATTTGGGCATTTGGG')


class DnaSpeciesTypeTestCase(unittest.TestCase):
    def test(self):
        dna = core.DnaSpeciesType(id='dna1', name='dna1', seq=Bio.Seq.Seq('ACGTACGT', alphabet=Bio.Alphabet.DNAAlphabet()),
                                  circular=False, double_stranded=False)

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
        dna1 = core.DnaSpeciesType(id='dna1', seq=Bio.Seq.Seq(
            'ACGTACGTACGTACG', alphabet=Bio.Alphabet.DNAAlphabet()))
        tu1 = core.TranscriptionUnitLocus(
            id='tu1', polymer=dna1, start=1, end=15)
        rna1 = core.RnaSpeciesType(id='rna1', name='rna1', transcription_units=[
                                   tu1], type=1, concentration=1, half_life=2)

        self.assertEqual(rna1.id, 'rna1')
        self.assertEqual(rna1.name, 'rna1')
        self.assertEqual(rna1.transcription_units, [tu1])
        self.assertEqual(rna1.type, 1)
        self.assertEqual(rna1.concentration, 1)
        self.assertEqual(rna1.half_life, 2)

    def test_get_empirical_formula(self):
        dna1 = core.DnaSpeciesType(id='dna1', seq=Bio.Seq.Seq(
            'A', alphabet=Bio.Alphabet.DNAAlphabet()))
        tu1 = core.TranscriptionUnitLocus(
            id='tu1', polymer=dna1, start=1, end=1)
        rna1 = core.RnaSpeciesType(
            id='rna1', name='rna1', transcription_units=[tu1])
        self.assertEqual(rna1.get_empirical_formula(),
                         chem.EmpiricalFormula('C10H12N5O7P'))

        dna1 = core.DnaSpeciesType(id='dna1', seq=Bio.Seq.Seq(
            'C', alphabet=Bio.Alphabet.DNAAlphabet()))
        tu1 = core.TranscriptionUnitLocus(
            id='tu1', polymer=dna1, start=1, end=1)
        rna1 = core.RnaSpeciesType(
            id='rna1', name='rna1', transcription_units=[tu1])
        self.assertEqual(rna1.get_empirical_formula(),
                         chem.EmpiricalFormula('C9H12N3O8P'))

        dna1 = core.DnaSpeciesType(id='dna1', seq=Bio.Seq.Seq(
            'G', alphabet=Bio.Alphabet.DNAAlphabet()))
        tu1 = core.TranscriptionUnitLocus(
            id='tu1', polymer=dna1, start=1, end=1)
        rna1 = core.RnaSpeciesType(
            id='rna1', name='rna1', transcription_units=[tu1])
        self.assertEqual(rna1.get_empirical_formula(),
                         chem.EmpiricalFormula('C10H12N5O8P'))

        dna1 = core.DnaSpeciesType(id='dna1', seq=Bio.Seq.Seq(
            'T', alphabet=Bio.Alphabet.DNAAlphabet()))
        tu1 = core.TranscriptionUnitLocus(
            id='tu1', polymer=dna1, start=1, end=1)
        rna1 = core.RnaSpeciesType(
            id='rna1', name='rna1', transcription_units=[tu1])
        self.assertEqual(rna1.get_empirical_formula(),
                         chem.EmpiricalFormula('C9H11N2O9P'))

        dna1 = core.DnaSpeciesType(id='dna1', seq=Bio.Seq.Seq(
            'AAAA', alphabet=Bio.Alphabet.DNAAlphabet()))
        tu1 = core.TranscriptionUnitLocus(
            id='tu1', polymer=dna1, start=1, end=2)
        rna1 = core.RnaSpeciesType(
            id='rna1', name='rna1', transcription_units=[tu1])
        self.assertEqual(rna1.get_empirical_formula(),
                         chem.EmpiricalFormula('C20H23N10O13P2'))

    def test_get_charge(self):
        dna1 = core.DnaSpeciesType(id='dna1', seq=Bio.Seq.Seq(
            'AAAA', alphabet=Bio.Alphabet.DNAAlphabet()))
        tu1 = core.TranscriptionUnitLocus(
            id='tu1', polymer=dna1, start=1, end=1)
        rna1 = core.RnaSpeciesType(
            id='rna1', name='rna1', transcription_units=[tu1])
        self.assertEqual(rna1.get_charge(), -2)

        dna1 = core.DnaSpeciesType(id='dna1', seq=Bio.Seq.Seq(
            'AAAA', alphabet=Bio.Alphabet.DNAAlphabet()))
        tu1 = core.TranscriptionUnitLocus(
            id='tu1', polymer=dna1, start=1, end=2)
        rna1 = core.RnaSpeciesType(
            id='rna1', name='rna1', transcription_units=[tu1])
        self.assertEqual(rna1.get_charge(), -3)

    def test_get_mol_wt(self):
        dna1 = core.DnaSpeciesType(id='dna1', seq=Bio.Seq.Seq(
            'AACCGGTT', alphabet=Bio.Alphabet.DNAAlphabet()))
        tu1 = core.TranscriptionUnitLocus(
            id='tu1', polymer=dna1, start=1, end=1)
        rna1 = core.RnaSpeciesType(
            id='rna1', name='rna1', transcription_units=[tu1])
        exp_mol_wt = \
            + Bio.SeqUtils.molecular_weight(rna1.get_seq()) \
            - (rna1.get_len() + 1) * mendeleev.element('H').atomic_weight
        self.assertAlmostEqual(rna1.get_mol_wt(), exp_mol_wt, places=1)

        dna1 = core.DnaSpeciesType(id='dna1', seq=Bio.Seq.Seq(
            'AACCGGTT', alphabet=Bio.Alphabet.DNAAlphabet()))
        tu1 = core.TranscriptionUnitLocus(
            id='tu1', polymer=dna1, start=3, end=3)
        rna1 = core.RnaSpeciesType(
            id='rna1', name='rna1', transcription_units=[tu1])
        exp_mol_wt = \
            + Bio.SeqUtils.molecular_weight(rna1.get_seq()) \
            - (rna1.get_len() + 1) * mendeleev.element('H').atomic_weight
        self.assertAlmostEqual(rna1.get_mol_wt(), exp_mol_wt, places=1)

        dna1 = core.DnaSpeciesType(id='dna1', seq=Bio.Seq.Seq(
            'AACCGGTT', alphabet=Bio.Alphabet.DNAAlphabet()))
        tu1 = core.TranscriptionUnitLocus(
            id='tu1', polymer=dna1, start=5, end=5)
        rna1 = core.RnaSpeciesType(
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
        gene1 = core.GeneLocus(id='gene1', cell=cell1,
                               polymer=dna1, start=692, end=1834)
        tu1 = core.TranscriptionUnitLocus(
            id='tu1', genes=[gene1], polymer=dna1)
        self.prot1 = core.ProteinSpeciesType(
            id='prot1', gene=gene1, cell=cell1)

        # MPN011
        gene2 = core.GeneLocus(id='gene2', cell=cell1, polymer=dna1,
                               start=12838, end=13533, strand=core.PolymerStrand.negative)
        tu2 = core.TranscriptionUnitLocus(
            id='tu2', genes=[gene2], polymer=dna1)
        self.prot2 = core.ProteinSpeciesType(
            id='prot2', gene=gene2, cell=cell1)

    def test_constructor(self):
        protein = core.ProteinSpeciesType(
            id='prot1', name='prot1', concentration=1, half_life=2)
        # attribute_order = ('id', 'cell', 'name', 'gene', 'rna', 'concentration', 'half_life', 'comments')

        self.assertEqual(protein.id, 'prot1')
        self.assertEqual(protein.name, 'prot1')
        self.assertEqual(protein.concentration, 1)
        self.assertEqual(protein.half_life, 2)
        self.assertEqual(protein.cell, None)

    def test_get_seq(self):

        # Use table 4 since example genes are from mycoplasma genitallium

        # MPN001
        self.assertEqual(self.prot1.get_seq()[0:10], 'MKVLINKNEL')

        # MPN011
        self.assertEqual(self.prot2.get_seq()[0:10], 'MKFKFLLTPL')

    def test_get_empirical_formula(self):
        # MPN001
        self.assertEqual(self.prot1.get_empirical_formula(),
                         chem.EmpiricalFormula('C1980H3146N510O596S7'))
        # MPN011
        self.assertEqual(self.prot2.get_empirical_formula(),
                         chem.EmpiricalFormula('C1246H1928N306O352S3'))

    def test_get_mol_wt(self):

        # MPN001
        self.assertAlmostEqual(self.prot1.get_mol_wt(), 43856.113, delta=0.3)
        # MNP011
        self.assertAlmostEqual(self.prot2.get_mol_wt(), 26922.957, delta=0.3)

    def test_get_charge(self):
        self.assertEqual(self.prot1.get_charge(), 1)

        self.assertEqual(self.prot2.get_charge(), 12)


class PolymerLocusTestCase(unittest.TestCase):
    def test_constructor(self):
        cell1 = core.Cell()
        dna1 = core.DnaSpeciesType(id='dna1', seq=Bio.Seq.Seq('ACGTACGTACGTACG', alphabet=Bio.Alphabet.DNAAlphabet()),
                                   circular=False, double_stranded=False)

        locus1 = core.PolymerLocus(id='locus1', cell=cell1, name='locus1', polymer=dna1,
                                   strand=core.PolymerStrand.positive, start=1, end=15)

        # test constructor
        self.assertEqual(locus1.id, 'locus1')
        self.assertEqual(locus1.cell, cell1)
        self.assertEqual(locus1.name, 'locus1')
        self.assertEqual(locus1.polymer, dna1)
        self.assertEqual(locus1.strand, core.PolymerStrand.positive)
        self.assertEqual(locus1.start, 1)
        self.assertEqual(locus1.end, 15)

        # test methods
        self.assertEqual(locus1.get_seq(), 'ACGTACGTACGTACG')
        self.assertEqual(locus1.get_len(), 15)

        # flip strand; test methods
        rev_comp_seq = locus1.get_seq().reverse_complement()
        locus1.strand = core.PolymerStrand.negative
        self.assertEqual(locus1.get_seq(), rev_comp_seq)
        self.assertEqual(locus1.get_len(), 15)


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
        self.assertEqual(met.get_empirical_formula(),
                         chem.EmpiricalFormula('C10H12N5O7P'))
        self.assertEqual(met.get_charge(), -2)
        self.assertAlmostEqual(met.get_mol_wt(), 345.20530, places=4)


class PromoterLocusTestCase(unittest.TestCase):
    def test_constructor(self):
        promoter = core.PromoterLocus(
            id='promoter1', name='promoter1', pribnow_start=1, pribnow_end=2)
        self.assertEqual(promoter.id, 'promoter1')
        self.assertEqual(promoter.name, 'promoter1')
        self.assertEqual(promoter.pribnow_start, 1)
        self.assertEqual(promoter.pribnow_end, 2)


class GeneLocusTestCase(unittest.TestCase):
    def test(self):
        gene = core.GeneLocus(id='gene1', name='gene1',
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

        tu1 = core.TranscriptionUnitLocus(
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


class ReactionAndRelatedClassesCase(unittest.TestCase):

    def setUp(self):
        self.cell_1 = cell_1 = core.Cell()
        self.compartment_1 = compartment_1 = core.Compartment(cell=cell_1, id='compartment_1')
        self.species_type_1 = species_type_1 = core.MetaboliteSpeciesType(id='species_type_1')
        self.species_type_2 = species_type_2 = core.MetaboliteSpeciesType(id='species_type_2')
        self.species_1 = species_1 = core.Species(
            species_type=species_type_1, compartment=compartment_1)
        self.species_2 = species_2 = core.Species(
            species_type=species_type_2, compartment=compartment_1)
        self.participant_1 = participant_1 = core.SpeciesCoefficient(
            species=species_1, coefficient=1)
        self.participant_2 = participant_2 = core.SpeciesCoefficient(
            species=species_2, coefficient=1)

        self.reaction_1 = reaction_1 = core.Reaction(
            id='reaction_1',
            name='test_reaction',
            cell=cell_1,
            participants=[participant_1, participant_2],
            reversible=False)

        self.rate_law_equation_1 = rate_law_equation_1 = core.RateLawEquation(
            id='rate_law_equation_1',
            expression=species_1.id()
        )

        self.rate_law_1 = rate_law_1 = core.RateLaw(
            id='rate_law_1',
            reaction=reaction_1,
            direction=core.RateLawDirection.forward,
            k_m=0.1,
            k_cat=0.5,
            equation=rate_law_equation_1
        )

        self.objects = {
            core.Species: {
                species_1.id(): species_1,
                species_2.id(): species_2
            }
        }

    def test_Reaction(self):
        self.assertEqual(self.reaction_1.id, 'reaction_1')
        self.assertEqual(self.reaction_1.name, 'test_reaction')
        self.assertEqual(self.reaction_1.cell, self.cell_1)
        self.assertEqual(self.reaction_1.participants, [self.participant_1, self.participant_2])
        self.assertEqual(self.reaction_1.reversible, 0)

    def test_RateLawEquation(self):
        self.assertEqual(self.rate_law_equation_1.id, 'rate_law_equation_1')
        self.assertEqual(self.rate_law_equation_1.serialize(), self.species_1.id())

    def test_RateLaw(self):
        self.assertEqual(self.rate_law_1.id, 'rate_law_1')
        self.assertEqual(self.rate_law_1.direction, core.RateLawDirection.forward)
        self.assertEqual(self.rate_law_1.equation, self.rate_law_equation_1)
        self.assertEqual(self.rate_law_1.k_m, 0.1)
        self.assertEqual(self.rate_law_1.k_cat, 0.5)
        self.assertIn(self.reaction_1.id, self.rate_law_1.serialize())
        self.assertIn('forward', self.rate_law_1.serialize())

    def test_deserialize_RateLawEquation(self):
        attr = core.RateLawEquation.expression

        # add RateLawEquation to self.objects
        rle_1, error = core.RateLawEquation.deserialize(attr, self.species_1.id(), self.objects)
        self.assertTrue(error is None)
        self.assertEqual(rle_1.modifiers, [self.species_1])

        # RateLawEquation is in self.objects, and value is in self.objects[RateLawEquation]
        rle_2, error = core.RateLawEquation.deserialize(attr, self.species_1.id(), self.objects)
        self.assertTrue(error is None)
        self.assertEqual(rle_2.modifiers, [self.species_1])
        self.assertEqual(rle_1, rle_2)

        # serialization errors
        # value not string; re.findall raises exception
        value = 123
        obj, error = core.RateLawEquation.deserialize(attr, value, self.objects)
        self.assertTrue(obj is None)

        # Species.deserialize fails
        value = 'not_species_id[x]'
        obj, error = core.RateLawEquation.deserialize(attr, value, self.objects)
        self.assertTrue(obj is None)


class ComplexSpeciesTypeTestCase(unittest.TestCase):
    def test_ComplexSpeciesType(self):

        # Test constructor
        complex1 = core.ComplexSpeciesType()

        self.assertEqual(complex1.region, '')
        self.assertEqual(complex1.binding, '')
        self.assertEqual(complex1.complex_type, '')
        self.assertEqual(complex1.formation_process, None)
        self.assertEqual(complex1.subunits, [])

        # Generate test proteins from  Mycoplasma Genintalium Genome
        records = Bio.SeqIO.parse('tests/fixtures/seq.fna', 'fasta')
        dna_seq = next(records).seq
        dna1 = core.DnaSpeciesType(seq=dna_seq)

        cell1 = dna1.cell = core.Cell()
        cell1.knowledge_base = core.KnowledgeBase(
            translation_table=4)  # Table 4 is for mycoplasma

        # Protein 1,  MPN001
        gene1 = core.GeneLocus(id='gene1', cell=cell1,
                               polymer=dna1, start=692, end=1834)
        tu1 = core.TranscriptionUnitLocus(
            id='tu1', genes=[gene1], polymer=dna1)
        prot1 = core.ProteinSpeciesType(
            id='prot1', gene=gene1, cell=cell1)

        # Protein 2, MPN011
        gene2 = core.GeneLocus(id='gene2', cell=cell1, polymer=dna1,
                               start=12838, end=13533, strand=core.PolymerStrand.negative)
        tu2 = core.TranscriptionUnitLocus(
            id='tu2', genes=[gene2], polymer=dna1)
        prot2 = core.ProteinSpeciesType(
            id='prot2', gene=gene2, cell=cell1)

        # Test adding formation reaction
        # Add formation reaction: [c]: (2) prot1 + (3) prot2 ==> complex1
        comp1 = core.Compartment(id='comp1')
        species1 = core.Species(species_type=prot1, compartment=comp1)
        species2 = core.Species(species_type=prot2, compartment=comp1)
        species_coeff1 = core.SpeciesCoefficient(
            species=species1, coefficient=2)
        species_coeff2 = core.SpeciesCoefficient(
            species=species2, coefficient=3)
        complex1.subunits = [species_coeff1, species_coeff2]

        self.assertEqual(complex1.get_charge(), 38)
        self.assertAlmostEqual(complex1.get_mol_wt(),
                               (2*prot1.get_mol_wt() + 3 * prot2.get_mol_wt()))
        self.assertEqual(complex1.get_empirical_formula(),
                         chem.EmpiricalFormula('C7698H12076N1938O2248S23'))


class SpeciesTestCase(unittest.TestCase):
    def test_SpeciesType(self):
        comp1 = core.Compartment(id='c')
        prot1 = core.ProteinSpeciesType(id='prot1')
        species1 = core.Species(species_type=prot1, compartment=comp1)

        self.assertEqual(core.Species.gen_id(prot1, comp1), 'prot1[c]')
        self.assertEqual(core.Species.gen_id('prot1', 'c'), 'prot1[c]')
        with self.assertRaisesRegexp(ValueError, 'incorrect species type'):
            core.Species.gen_id(None, 'c')
        with self.assertRaisesRegexp(ValueError, 'incorrect compartment type'):
            core.Species.gen_id('prot1', None)

        self.assertEqual(species1.id(), 'prot1[c]')

    def test_serialize(self):
        comp1 = core.Compartment(id='c')
        prot1 = core.ProteinSpeciesType(id='prot1')
        species1 = core.Species(species_type=prot1, compartment=comp1)

        self.assertEqual(species1.serialize(), 'prot1[c]')

    def test_deserialize(self):
        comp1 = core.Compartment(id='c')
        prot1 = core.ProteinSpeciesType(id='prot1')
        rna1 = core.RnaSpeciesType(id='rna1')

        objects = {
            core.Compartment: {
                'c': comp1,
            },
            core.ProteinSpeciesType: {
                'prot1': prot1,
            },
            core.RnaSpeciesType: {
                'rna1': rna1,
            },
        }

        attr = core.SpeciesCoefficient.species
        result = core.Species.deserialize(attr, 'prot1[c]', objects)
        self.assertEqual(result[0].species_type, prot1)
        self.assertEqual(result[0].compartment, comp1)
        self.assertEqual(result[1], None)

        result2 = core.Species.deserialize(attr, 'prot1[c]', objects)
        self.assertEqual(result2[0], result[0])
        self.assertIn(core.Species, objects)
        self.assertIn('prot1[c]', objects[core.Species])

        self.assertNotIn('rna1[c]', objects[core.Species])
        self.assertEqual(core.Species.deserialize(
            attr, 'rna1[c]', objects)[1], None)
        self.assertIn('rna1[c]', objects[core.Species])

        self.assertNotEqual(core.Species.deserialize(
            attr, 'prot2[c]', objects)[1], None)
        self.assertNotEqual(core.Species.deserialize(
            attr, 'prot1[e]', objects)[1], None)
        self.assertNotEqual(core.Species.deserialize(
            attr, 'prot1', objects)[1], None)


class SpeciesCoefficientTestCase(unittest.TestCase):
    def test_constructor(self):
        comp = core.Compartment(id='c')
        prot = core.ProteinSpeciesType(id='prot')

        spec = core.Species(species_type=prot, compartment=comp)
        spec_coeff = core.SpeciesCoefficient(species=spec, coefficient=3)

    def test_serialize(self):
        comp = core.Compartment(id='c')
        prot = core.ProteinSpeciesType(id='prot')

        spec = core.Species(species_type=prot, compartment=comp)
        spec_coeff = core.SpeciesCoefficient(species=spec, coefficient=3)

        self.assertEqual(spec_coeff.serialize(), '(3) prot[c]')
        self.assertEqual(core.SpeciesCoefficient._serialize(
            species=spec, coefficient=3), '(3) prot[c]')
        self.assertEqual(core.SpeciesCoefficient._serialize(
            species=spec, coefficient=1), 'prot[c]')
        self.assertEqual(core.SpeciesCoefficient._serialize(
            species=spec, coefficient=2000), '(2.000000e+03) prot[c]')
        self.assertEqual(core.SpeciesCoefficient._serialize(
            species=spec, coefficient=3, show_compartment=False), '(3) prot')
        self.assertEqual(core.SpeciesCoefficient._serialize(
            species=spec, coefficient=-1), '(-1) prot[c]')
        self.assertEqual(core.SpeciesCoefficient._serialize(
            species=spec, coefficient=-1, show_coefficient_sign=False), 'prot[c]')

    def test_deserialize(self):
        comp = core.Compartment(id='c')
        rna = core.RnaSpeciesType(id='rna')
        prot = core.ProteinSpeciesType(id='prot')

        spec = core.Species(species_type=rna, compartment=comp)

        objects = {
            core.Compartment: {
                'c': comp,
            },
            core.RnaSpeciesType: {
                'rna': rna,
            },
            core.ProteinSpeciesType: {
                'prot': prot,
            },
            core.Species: {
                'rna[c]': spec,
            },
        }

        attr = core.Reaction.participants

        result = core.SpeciesCoefficient.deserialize(
            attr, '(3) rna[c]', objects)
        self.assertEqual(result[0].species, spec)
        self.assertEqual(result[0].coefficient, 3)
        self.assertEqual(result[1], None)
        self.assertEqual(len(objects[core.Species]), 1)
        self.assertEqual(len(objects[core.SpeciesCoefficient]), 1)

        result2 = core.SpeciesCoefficient.deserialize(
            attr, '(3) rna[c]', objects)
        self.assertEqual(result2[0], result[0])
        self.assertEqual(result2[1], None)
        self.assertEqual(len(objects[core.Species]), 1)
        self.assertEqual(len(objects[core.SpeciesCoefficient]), 1)

        result = core.SpeciesCoefficient.deserialize(
            attr, '(-2) rna[c]', objects)
        self.assertEqual(result[0].species, spec)
        self.assertEqual(result[0].coefficient, -2)
        self.assertEqual(result[1], None)
        self.assertEqual(len(objects[core.Species]), 1)
        self.assertEqual(len(objects[core.SpeciesCoefficient]), 2)

        result = core.SpeciesCoefficient.deserialize(attr, 'rna[c]', objects)
        self.assertEqual(result[0].species, spec)
        self.assertEqual(result[0].coefficient, 1)
        self.assertEqual(result[1], None)
        self.assertEqual(len(objects[core.Species]), 1)
        self.assertEqual(len(objects[core.SpeciesCoefficient]), 3)

        result = core.SpeciesCoefficient.deserialize(attr, 'prot[c]', objects)
        self.assertEqual(result[0].species.species_type, prot)
        self.assertEqual(result[0].species.compartment, comp)
        self.assertEqual(result[0].coefficient, 1)
        self.assertEqual(result[1], None)
        self.assertEqual(len(objects[core.Species]), 2)
        self.assertEqual(len(objects[core.SpeciesCoefficient]), 4)

        result = core.SpeciesCoefficient.deserialize(attr, 'rna2[c]', objects)
        self.assertEqual(result[0], None)
        self.assertNotEqual(result[1], None)

        result = core.SpeciesCoefficient.deserialize(attr, 'rna', objects)
        self.assertEqual(result[0], None)
        self.assertNotEqual(result[1], None)

        result = core.SpeciesCoefficient.deserialize(
            attr, 'rna', objects, compartment=comp)
        self.assertEqual(result[0].species.species_type, rna)
        self.assertEqual(result[0].species.compartment, comp)
        self.assertEqual(result[0].coefficient, 1)
        self.assertEqual(result[1], None)

        result = core.SpeciesCoefficient.deserialize(
            attr, '(2) rna', objects, compartment=comp)
        self.assertEqual(result[0].species.species_type, rna)
        self.assertEqual(result[0].species.compartment, comp)
        self.assertEqual(result[0].coefficient, 2)
        self.assertEqual(result[1], None)

        result = core.SpeciesCoefficient.deserialize(
            attr, '(-3) rna', objects, compartment=comp)
        self.assertEqual(result[0].species.species_type, rna)
        self.assertEqual(result[0].species.compartment, comp)
        self.assertEqual(result[0].coefficient, -3)
        self.assertEqual(result[1], None)


class SubunitAttributeTestCase(unittest.TestCase):
    def test_SubunitAttribute(self):
        compart1 = core.Compartment(id='c')
        compart2 = core.Compartment(id='m')

        prot1 = core.ProteinSpeciesType(id='prot1')
        prot2 = core.ProteinSpeciesType(id='prot2')
        complex1 = core.ComplexSpeciesType(id='complex1')

        species1 = core.Species(species_type=prot1, compartment=compart1)
        species2 = core.Species(species_type=prot2, compartment=compart1)
        species3 = core.Species(species_type=complex1, compartment=compart1)
        species4 = core.Species(species_type=complex1, compartment=compart2)

        species_coeff1 = core.SpeciesCoefficient(
            species=species1, coefficient=2)
        species_coeff2 = core.SpeciesCoefficient(
            species=species2, coefficient=3)
        species_coeff3 = core.SpeciesCoefficient(
            species=species3, coefficient=5)
        species_coeff4 = core.SpeciesCoefficient(
            species=species4, coefficient=7)

        self.assertEqual(
            core.SubunitAttribute().serialize(participants=[]), '')
        self.assertEqual(core.SubunitAttribute().serialize(participants=[species_coeff1]),
                         '[c]: (2) prot1')
        self.assertEqual(core.SubunitAttribute().serialize(participants=[species_coeff1, species_coeff2]),
                         '[c]: (2) prot1 + (3) prot2')
        self.assertEqual(core.SubunitAttribute().serialize(participants=[species_coeff1, species_coeff2, species_coeff3]),
                         '[c]: (5) complex1 + (2) prot1 + (3) prot2')
        self.assertEqual(core.SubunitAttribute().serialize(participants=[species_coeff1, species_coeff2, species_coeff4]),
                         '(7) complex1[m] + (2) prot1[c] + (3) prot2[c]')

        objects = {
            core.DnaSpeciesType: {},
            core.RnaSpeciesType: {},
            core.MetaboliteSpeciesType: {},
            core.ComplexSpeciesType: {
                'complex1': complex1
            },
            core.Compartment: {
                'c': compart1, 'm': compart2
            },
            core.ProteinSpeciesType: {
                'prot1': prot1, 'prot2': prot2
            },
            core.Species: {
                'prot1[c]': species1, 'prot2[c]': species2, 'complex1[c]': species3, 'complex[m]': species4
            },
        }

        result = core.SubunitAttribute().deserialize(
            value='[c]: prot1 + (2) prot2', objects=objects)
        self.assertEqual(result[0][0].species.species_type, prot1)
        self.assertEqual(result[0][1].species.species_type, prot2)
        self.assertEqual(result[0][0].coefficient, 1)
        self.assertEqual(result[0][1].coefficient, 2)
        self.assertEqual(result[0][0].species.compartment, compart1)
        self.assertEqual(result[0][1].species.compartment, compart1)
        self.assertEqual(result[0][0].species.id(), 'prot1[c]')
        self.assertEqual(result[0][1].species.id(), 'prot2[c]')
        self.assertEqual(result[1], None)

        result = core.SubunitAttribute().deserialize(
            value='(2) prot1[c] + (7) complex1[m]', objects=objects)
        self.assertEqual(result[0][0].species.species_type, prot1)
        self.assertEqual(result[0][1].species.species_type, complex1)
        self.assertEqual(result[0][0].coefficient, 2)
        self.assertEqual(result[0][1].coefficient, 7)
        self.assertEqual(result[0][0].species.compartment, compart1)
        self.assertEqual(result[0][1].species.compartment, compart2)
        self.assertEqual(result[0][0].species.id(), 'prot1[c]')
        self.assertEqual(result[0][1].species.id(), 'complex1[m]')
        self.assertEqual(result[1], None)

        result = core.SubunitAttribute().deserialize(
            value='[e]: prot1 + (2) prot2', objects=objects)
        self.assertEqual(result[0], None)
        self.assertEqual(result[1].messages[0], 'Undefined compartment "e"')

        result = core.SubunitAttribute().deserialize(
            value='[c]prot1 + (2) prot2', objects=objects)
        self.assertEqual(result[0], None)
        self.assertEqual(
            result[1].messages[0], 'Incorrectly formatted participants: [c]prot1 + (2) prot2')


class ReactionParticipantAttributeTestCase(unittest.TestCase):
    def test_ReactionParticipantAttribute(self):
        compart1 = core.Compartment(id='c')
        compart2 = core.Compartment(id='m')

        prot1 = core.ProteinSpeciesType(id='prot1')
        prot2 = core.ProteinSpeciesType(id='prot2')
        complex1 = core.ComplexSpeciesType(id='complex1')

        species1 = core.Species(species_type=prot1, compartment=compart1)
        species2 = core.Species(species_type=prot2, compartment=compart1)
        species3 = core.Species(species_type=complex1, compartment=compart1)
        species4 = core.Species(species_type=complex1, compartment=compart2)

        species_coeff1 = core.SpeciesCoefficient(
            species=species1, coefficient=-2)
        species_coeff2 = core.SpeciesCoefficient(
            species=species2, coefficient=-3)
        species_coeff3 = core.SpeciesCoefficient(
            species=species3, coefficient=5)
        species_coeff4 = core.SpeciesCoefficient(
            species=species4, coefficient=7)

        self.assertEqual(
            core.ReactionParticipantAttribute().serialize(participants=[]), '')
        self.assertEqual(core.ReactionParticipantAttribute().serialize(participants=[species_coeff1]),
                         '[c]: (2) prot1 ==> ')
        self.assertEqual(core.ReactionParticipantAttribute().serialize(participants=[species_coeff1, species_coeff2]),
                         '[c]: (2) prot1 + (3) prot2 ==> ')
        self.assertEqual(core.ReactionParticipantAttribute().serialize(participants=[species_coeff1, species_coeff2, species_coeff3]),
                         '[c]: (2) prot1 + (3) prot2 ==> (5) complex1')
        self.assertEqual(core.ReactionParticipantAttribute().serialize(participants=[species_coeff1, species_coeff2, species_coeff4]),
                         '(2) prot1[c] + (3) prot2[c] ==> (7) complex1[m]')

        objects = {
            core.DnaSpeciesType: {},
            core.RnaSpeciesType: {},
            core.MetaboliteSpeciesType: {},
            core.ComplexSpeciesType: {
                'complex1': complex1
            },
            core.Compartment: {
                'c': compart1, 'm': compart2
            },
            core.ProteinSpeciesType: {
                'prot1': prot1, 'prot2': prot2
            },
            core.Species: {
                'prot1[c]': species1, 'prot2[c]': species2, 'complex1[c]': species3, 'complex[m]': species4

            },
        }

        result = core.ReactionParticipantAttribute().deserialize(
            value='[c]: prot1 ==> prot2', objects=objects)
        self.assertEqual(result[0][0].species.species_type, prot1)
        self.assertEqual(result[0][1].species.species_type, prot2)
        self.assertEqual(result[0][0].coefficient, -1)
        self.assertEqual(result[0][1].coefficient,  1)
        self.assertEqual(result[0][0].species.compartment, compart1)
        self.assertEqual(result[0][1].species.compartment, compart1)
        self.assertEqual(result[0][0].species.id(), 'prot1[c]')
        self.assertEqual(result[0][1].species.id(), 'prot2[c]')
        self.assertEqual(result[1], None)

        result = core.ReactionParticipantAttribute().deserialize(
            value='(2) prot1[c] + (3) prot2[c] ==> (7) complex1[m]', objects=objects)
        self.assertEqual(result[0][0].species.species_type, prot1)
        self.assertEqual(result[0][1].species.species_type, prot2)
        self.assertEqual(result[0][2].species.species_type, complex1)
        self.assertEqual(result[0][0].coefficient, -2)
        self.assertEqual(result[0][1].coefficient, -3)
        self.assertEqual(result[0][2].coefficient, 7.0)
        self.assertEqual(result[0][0].species.compartment, compart1)
        self.assertEqual(result[0][1].species.compartment, compart1)
        self.assertEqual(result[0][2].species.compartment, compart2)
        self.assertEqual(result[0][0].species.id(), 'prot1[c]')
        self.assertEqual(result[0][1].species.id(), 'prot2[c]')
        self.assertEqual(result[0][2].species.id(), 'complex1[m]')
        self.assertEqual(result[1], None)

        result = core.ReactionParticipantAttribute().deserialize(
            value='prot2[c] ==>', objects=objects)
        self.assertEqual(result[0], None)
        self.assertEqual(
            result[1].messages[0], 'Incorrectly formatted participants: prot2[c] ==>')

        result = core.ReactionParticipantAttribute().deserialize(
            value='==> prot1[c]', objects=objects)
        self.assertEqual(result[0], None)
        self.assertEqual(
            result[1].messages[0], 'Incorrectly formatted participants: ==> prot1[c]')


class ObservableCoefficientTestCase(unittest.TestCase):
    def test_observable_coefficient(self):
        cell = core.Cell()
        comp1 = core.Compartment(id='c')
        prot1 = core.ProteinSpeciesType(id='prot1')
        rna1 = core.RnaSpeciesType(id="rna1")

        species1 = core.Species(species_type=prot1, compartment=comp1)
        species2 = core.Species(species_type=rna1, compartment=comp1)

        speciesCoefficient1 = core.SpeciesCoefficient(
            species=species1, coefficient=2.5)
        speciesCoefficient2 = core.SpeciesCoefficient(
            species=species2, coefficient=3)

        observable1 = core.Observable(
            id='test', cell=cell, species=[speciesCoefficient1, speciesCoefficient2])
        observableCoefficinet1 = core.ObservableCoefficient(
            observable=observable1, coefficient=2.3)

        self.assertIsInstance(observableCoefficinet1,
                              core.ObservableCoefficient)
        self.assertIsInstance(
            observableCoefficinet1.observable, core.Observable)
        self.assertIsInstance(observableCoefficinet1.coefficient, float)
        self.assertEqual(observableCoefficinet1.observable.id, 'test')


class ObservableTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def test_observables(self):
        cell = core.Cell()
        comp1 = core.Compartment(id='c')
        prot1 = core.ProteinSpeciesType(id='prot1')
        rna1 = core.RnaSpeciesType(id="rna1")

        species1 = core.Species(species_type=prot1, compartment=comp1)
        species2 = core.Species(species_type=rna1, compartment=comp1)

        speciesCoefficient1 = core.SpeciesCoefficient(
            species=species1, coefficient=2)
        speciesCoefficient2 = core.SpeciesCoefficient(
            species=species2, coefficient=3.3)

        observable1 = core.Observable(
            cell=cell, species=[speciesCoefficient1, speciesCoefficient2])
        observableCoefficinet1 = core.ObservableCoefficient(
            observable=observable1, coefficient=2)

        observable2 = core.Observable(cell=cell, species=[
                                      speciesCoefficient1, speciesCoefficient2], observables=[observableCoefficinet1])

        self.assertIsInstance(observable1, core.Observable)
        self.assertIsInstance(observable2, core.Observable)
        with self.assertRaisesRegex(AttributeError, ""):
            observable3 = core.Observable(
                cell=cell, species=[species1, species2])
        self.assertIsInstance(observable1.species[0], core.SpeciesCoefficient)
        self.assertIsInstance(
            observable2.observables[0], core.ObservableCoefficient)
        self.assertIsInstance(
            observable1.species[0].species, core.Species)
        self.assertEqual(observable1.species[0].species.id(), 'prot1[c]')
        self.assertEqual(
            observable1.species[1].species.species_type.id, 'rna1')
        self.assertIsInstance(
            observable2.observables[0].observable, core.Observable)
