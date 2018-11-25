""" Tests of the knowledge base core schema

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Author: Bilal Shaikh <bilal.shaikh@columbia.edu>
:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
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
import numpy
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


class DatabaseReferenceTestCase(unittest.TestCase):
    def test_constructor(self):
        db_ref1 = core.DatabaseReference(database='Sabio-Rk', id='123456')
        self.assertEqual(db_ref1.database, 'Sabio-Rk')
        self.assertEqual(db_ref1.id, '123456')


class ReferenceTestCase(unittest.TestCase):
    def test_constructor(self):
        ref1 = core.Reference(id='ref1', standard_id='10.1000/xyz123')
        self.assertEqual(ref1.id, 'ref1')
        self.assertEqual(ref1.standard_id, '10.1000/xyz123')


class KnowledgeBaseTestCase(unittest.TestCase):
    def test_constructor(self):
        kb = core.KnowledgeBase()
        self.assertEqual(kb.cell, None)
        

class CellTestCase(unittest.TestCase):
    def test_constructor(self):
        cell = core.Cell()

        self.assertEqual(cell.knowledge_base, None)
        self.assertEqual(cell.taxon, None)
        self.assertEqual(cell.observables, [])
        self.assertEqual(cell.species_types, [])
        self.assertEqual(cell.compartments, [])
        self.assertEqual(cell.reactions, [])
        self.assertEqual(cell.loci, [])

        self.assertEqual(cell.species_types.get(
            __type=core.DnaSpeciesType), [])
        

class CompartmentTestCase(unittest.TestCase):
    def test_constructor(self):
        comp = core.Compartment(volumetric_fraction=0.5)
        self.assertEqual(comp.volumetric_fraction, 0.5)


class SpeciesTypeTestCase(unittest.TestCase):
    def test_constructor(self):
        with self.assertRaisesRegex(TypeError, 'Can\'t instantiate abstract class'):
            core.SpeciesType()

        class ConcreteSpeciesType(core.SpeciesType):
            def get_empirical_formula(self): pass

            def get_charge(self): pass

            def get_mol_wt(self): pass

        species_type = ConcreteSpeciesType(
            id='species1', name='species1', half_life=3.)
        self.assertEqual(species_type.id, 'species1')
        self.assertEqual(species_type.name, 'species1')
        self.assertEqual(species_type.half_life, 3.)


class PolymerSpeciesTypeTestCase(unittest.TestCase):
    def test_constructor(self):
        with self.assertRaisesRegex(TypeError, 'Can\'t instantiate abstract class'):
            core.PolymerSpeciesType()

        class ConcretePolymerSpeciesType(core.PolymerSpeciesType):
            def get_empirical_formula(self): pass

            def get_charge(self): pass

            def get_mol_wt(self): pass

            def get_seq(self): return Bio.Seq.Seq(
                'AAATGCCC', alphabet=Bio.Alphabet.DNAAlphabet())

        pst1 = ConcretePolymerSpeciesType(
            id='pst1', name='pst1', half_life=2)

        # Test constructor
        self.assertEqual(pst1.id, 'pst1')
        self.assertEqual(pst1.name, 'pst1')
        self.assertEqual(pst1.half_life, 2)

        # Test methods: linear, single stranded case
        pst1.circular = False
        pst1.double_stranded = False

        self.assertEqual(pst1.get_len(), 8)
        self.assertEqual(pst1.get_subseq(1, 3), 'AAA')
        self.assertEqual(pst1.get_subseq(2, 4), 'AAT')

        with self.assertRaisesRegex(ValueError, 'Start and end coordinates'):
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
                                  circular=False, double_stranded=False, ploidy=2)

        self.assertEqual(dna.id, 'dna1')
        self.assertEqual(dna.name, 'dna1')
        self.assertEqual(dna.circular, False)
        self.assertEqual(dna.double_stranded, False)
        self.assertEqual(dna.ploidy, 2)

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


class ReactionAndRelatedClassesTestCase(unittest.TestCase):

    def setUp(self):
        self.cell_1 = cell_1 = core.Cell()
        self.compartment_1 = compartment_1 = core.Compartment(cell=cell_1, id='c')
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
            expression=species_1.id()
        )

        self.rate_law_1 = rate_law_1 = core.RateLaw(
            reaction=reaction_1,
            direction=core.RateLawDirection.forward,
            k_m=0.1,
            k_cat=0.5,
            equation=rate_law_equation_1
        )

        self.reaction_2 = reaction_2 = core.Reaction(
            id='reaction_2',
            name='test_reaction2',
            cell=cell_1,
            participants=[participant_1, participant_2],
            reversible=False)

        self.parameter_1 = parameter_1 = core.Parameter(
            id='p1',
            name='test_parameter1',
            value=0.7,
            error=0.01,
            units='s^-1'
        )

        self.parameter_2 = parameter_2 = core.Parameter(
            id='p2',
            name='test_parameter2',
            value=2.,
            error=0.1,
            units='M'
        )

        self.parameter_3 = parameter_3 = core.Parameter(
            id='p3',
            name='test_parameter3',
            value=2.3,
            error=0.15,
            units='M'
        )

        self.rate_law_equation_2 = rate_law_equation_2 = core.RateLawEquation(
            expression='p1*species_type_1[c]*species_type_2[c]/(p2+species_type_2[c]+(species_type_2[c]^2/p3))', 
            modifiers=[species_1, species_2], 
            parameters=[parameter_1, parameter_2]
        )

        self.rate_law_2 = rate_law_2 = core.RateLaw(
            reaction=reaction_2,
            direction=core.RateLawDirection.backward,
            equation=rate_law_equation_2
        )

        self.objects = {
            core.Species: {
                species_1.id(): species_1,
                species_2.id(): species_2
            },
            core.Parameter: {
                parameter_1.id: parameter_1,
                parameter_2.id: parameter_2,
                parameter_3.id: parameter_3
            }
        }

    def test_Reaction(self):
        self.assertEqual(self.reaction_1.id, 'reaction_1')
        self.assertEqual(self.reaction_1.name, 'test_reaction')
        self.assertEqual(self.reaction_1.cell, self.cell_1)
        self.assertEqual(self.reaction_1.participants, [self.participant_1, self.participant_2])
        self.assertEqual(self.reaction_1.reversible, 0)

    def test_RateLawEquation(self):
        self.assertEqual(self.rate_law_equation_1.serialize(), self.species_1.id())

    def test_RateLaw(self):
        self.assertEqual(self.rate_law_1.direction, core.RateLawDirection.forward)
        self.assertEqual(self.rate_law_1.equation, self.rate_law_equation_1)
        self.assertEqual(self.rate_law_1.k_m, 0.1)
        self.assertEqual(self.rate_law_1.k_cat, 0.5)
        self.assertIn(self.reaction_1.id, self.rate_law_1.serialize())
        self.assertIn('forward', self.rate_law_1.serialize())

        self.assertEqual(self.rate_law_2.direction, core.RateLawDirection.backward)
        self.assertEqual(self.rate_law_2.equation, self.rate_law_equation_2)
        self.assertEqual(numpy.isnan(self.rate_law_2.k_m), True)
        self.assertEqual(numpy.isnan(self.rate_law_2.k_cat), True)
        self.assertIn(self.reaction_2.id, self.rate_law_2.serialize())
        self.assertIn('backward', self.rate_law_2.serialize())

    def test_Parameter(self):
        self.assertEqual(self.parameter_1.id, 'p1')
        self.assertEqual(self.parameter_2.name, 'test_parameter2')
        self.assertEqual(self.parameter_3.value, 2.3)
        self.assertEqual(self.parameter_3.error, 0.15)
        self.assertEqual(self.parameter_3.units, 'M')
        self.assertEqual(self.parameter_3.references, [])
        self.assertEqual(self.parameter_3.database_references, [])                                                

    def test_deserialize_RateLawEquation(self):
        attr = core.RateLawEquation.expression

        # add RateLawEquation to self.objects
        rle_1, error = core.RateLawEquation.deserialize(attr, self.species_1.id(), self.objects)
        self.assertTrue(error is None)
        self.assertEqual(rle_1.modifiers, [self.species_1])
        rle_2, error = core.RateLawEquation.deserialize(attr, 
            'p1*species_type_1[c]*species_type_2[c]/(p2+species_type_2[c]+(species_type_2[c]^2/p3))', self.objects)
        self.assertTrue(error is None)
        self.assertEqual(rle_2.modifiers, [self.species_1, self.species_2])
        self.assertEqual(rle_2.parameters, [self.parameter_1, self.parameter_2, self.parameter_3])

        # RateLawEquation is in self.objects, and value is in self.objects[RateLawEquation]
        rle_3, error = core.RateLawEquation.deserialize(attr, self.species_1.id(), self.objects)
        self.assertTrue(error is None)
        self.assertEqual(rle_3.modifiers, [self.species_1])
        self.assertEqual(rle_1, rle_3)

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
        self.assertEqual(complex1.composition_in_uniprot, '')
        self.assertEqual(complex1.formation_process, None)
        self.assertEqual(complex1.subunits, [])

        cofactor1 = core.MetaboliteSpeciesType(id='cofactor1', 
            structure='InChI=1S/C8H7NO3/c10-6-1-4-5(2-7(6)11)9-3-8(4)12/h1-2,8-9,12H,3H2')
        cofactor2 = core.MetaboliteSpeciesType(id='cofactor2',
            structure='InChI=1S/Zn/q+2')

        # Test adding subunit composition 
        # Add subunit composition: (2) cofactor1 + (3) cofactor2 ==> complex1        
        species_type_coeff1 = core.SpeciesTypeCoefficient(
            species_type=cofactor1, coefficient=2)
        species_type_coeff2 = core.SpeciesTypeCoefficient(
            species_type=cofactor2, coefficient=3)
        complex1.subunits = [species_type_coeff1, species_type_coeff2]

        self.assertEqual(complex1.get_charge(), 6)
        self.assertAlmostEqual(complex1.get_mol_wt(),
                               (2*cofactor1.get_mol_wt() + 3 * cofactor2.get_mol_wt()))
        self.assertEqual(complex1.get_empirical_formula(),
                         chem.EmpiricalFormula('C16H14N2O6Zn3'))


class SpeciesTestCase(unittest.TestCase):
    def test_SpeciesType(self):
        comp1 = core.Compartment(id='c')
        met1 = core.MetaboliteSpeciesType(id='met1')
        species1 = core.Species(species_type=met1, compartment=comp1)

        self.assertEqual(core.Species.gen_id(met1, comp1), 'met1[c]')
        self.assertEqual(core.Species.gen_id('met1', 'c'), 'met1[c]')
        with self.assertRaisesRegex(ValueError, 'incorrect species type'):
            core.Species.gen_id(None, 'c')
        with self.assertRaisesRegex(ValueError, 'incorrect compartment type'):
            core.Species.gen_id('met1', None)

        self.assertEqual(species1.id(), 'met1[c]')

    def test_serialize(self):
        comp1 = core.Compartment(id='c')
        met1 = core.MetaboliteSpeciesType(id='met1')
        species1 = core.Species(species_type=met1, compartment=comp1)

        self.assertEqual(species1.serialize(), 'met1[c]')

    def test_deserialize(self):
        comp1 = core.Compartment(id='c')
        met1 = core.MetaboliteSpeciesType(id='met1')
        dna1 = core.DnaSpeciesType(id='dna1')

        objects = {
            core.Compartment: {
                'c': comp1,
            },
            core.MetaboliteSpeciesType: {
                'met1': met1,
            },
            core.DnaSpeciesType: {
                'dna1': dna1,
            },
        }

        attr = core.SpeciesCoefficient.species
        result = core.Species.deserialize(attr, 'met1[c]', objects)
        self.assertEqual(result[0].species_type, met1)
        self.assertEqual(result[0].compartment, comp1)
        self.assertEqual(result[1], None)

        result2 = core.Species.deserialize(attr, 'met1[c]', objects)
        self.assertEqual(result2[0], result[0])
        self.assertIn(core.Species, objects)
        self.assertIn('met1[c]', objects[core.Species])

        self.assertNotIn('dna1[c]', objects[core.Species])
        self.assertEqual(core.Species.deserialize(
            attr, 'dna1[c]', objects)[1], None)
        self.assertIn('dna1[c]', objects[core.Species])

        self.assertNotEqual(core.Species.deserialize(
            attr, 'met2[c]', objects)[1], None)
        self.assertNotEqual(core.Species.deserialize(
            attr, 'met1[e]', objects)[1], None)
        self.assertNotEqual(core.Species.deserialize(
            attr, 'met1', objects)[1], None)


class ConcentrationTestCase(unittest.TestCase):
    def test_constructor(self):
        comp = core.Compartment(id='c')
        met = core.MetaboliteSpeciesType(id='met')
        spec = core.Species(species_type=met, compartment=comp)

        conc = core.Concentration(species=spec, value=0.2)

        self.assertEqual(conc.serialize(), 'met[c]')
        self.assertEqual(conc.value, 0.2)
        self.assertEqual(conc.units, 2)


class SpeciesCoefficientTestCase(unittest.TestCase):
    def test_constructor(self):
        comp = core.Compartment(id='c')
        met = core.MetaboliteSpeciesType(id='met')

        spec = core.Species(species_type=met, compartment=comp)
        spec_coeff = core.SpeciesCoefficient(species=spec, coefficient=3)

        self.assertEqual(spec_coeff.species.species_type.id, 'met')
        self.assertEqual(spec_coeff.species.compartment.id, 'c')
        self.assertEqual(spec_coeff.coefficient, 3)

    def test_serialize(self):
        comp = core.Compartment(id='c')
        met = core.MetaboliteSpeciesType(id='met')

        spec = core.Species(species_type=met, compartment=comp)
        spec_coeff = core.SpeciesCoefficient(species=spec, coefficient=3)

        self.assertEqual(spec_coeff.serialize(), '(3) met[c]')
        self.assertEqual(core.SpeciesCoefficient._serialize(
            species=spec, coefficient=3), '(3) met[c]')
        self.assertEqual(core.SpeciesCoefficient._serialize(
            species=spec, coefficient=1), 'met[c]')
        self.assertEqual(core.SpeciesCoefficient._serialize(
            species=spec, coefficient=2000), '(2.000000e+03) met[c]')
        self.assertEqual(core.SpeciesCoefficient._serialize(
            species=spec, coefficient=3, show_compartment=False), '(3) met')
        self.assertEqual(core.SpeciesCoefficient._serialize(
            species=spec, coefficient=-1), '(-1) met[c]')
        self.assertEqual(core.SpeciesCoefficient._serialize(
            species=spec, coefficient=-1, show_coefficient_sign=False), 'met[c]')

    def test_deserialize(self):
        comp = core.Compartment(id='c')
        dna = core.DnaSpeciesType(id='dna')
        met = core.MetaboliteSpeciesType(id='met')

        spec = core.Species(species_type=dna, compartment=comp)

        objects = {
            core.Compartment: {
                'c': comp,
            },
            core.DnaSpeciesType: {
                'dna': dna,
            },
            core.MetaboliteSpeciesType: {
                'met': met,
            },
            core.Species: {
                'dna[c]': spec,
            },
        }

        attr = core.Reaction.participants

        result = core.SpeciesCoefficient.deserialize(
            attr, '(3) dna[c]', objects)
        self.assertEqual(result[0].species, spec)
        self.assertEqual(result[0].coefficient, 3)
        self.assertEqual(result[1], None)
        self.assertEqual(len(objects[core.Species]), 1)
        self.assertEqual(len(objects[core.SpeciesCoefficient]), 1)

        result2 = core.SpeciesCoefficient.deserialize(
            attr, '(3) dna[c]', objects)
        self.assertEqual(result2[0], result[0])
        self.assertEqual(result2[1], None)
        self.assertEqual(len(objects[core.Species]), 1)
        self.assertEqual(len(objects[core.SpeciesCoefficient]), 1)

        result = core.SpeciesCoefficient.deserialize(
            attr, '(-2) dna[c]', objects)
        self.assertEqual(result[0].species, spec)
        self.assertEqual(result[0].coefficient, -2)
        self.assertEqual(result[1], None)
        self.assertEqual(len(objects[core.Species]), 1)
        self.assertEqual(len(objects[core.SpeciesCoefficient]), 2)

        result = core.SpeciesCoefficient.deserialize(attr, 'dna[c]', objects)
        self.assertEqual(result[0].species, spec)
        self.assertEqual(result[0].coefficient, 1)
        self.assertEqual(result[1], None)
        self.assertEqual(len(objects[core.Species]), 1)
        self.assertEqual(len(objects[core.SpeciesCoefficient]), 3)

        result = core.SpeciesCoefficient.deserialize(attr, 'met[c]', objects)
        self.assertEqual(result[0].species.species_type, met)
        self.assertEqual(result[0].species.compartment, comp)
        self.assertEqual(result[0].coefficient, 1)
        self.assertEqual(result[1], None)
        self.assertEqual(len(objects[core.Species]), 2)
        self.assertEqual(len(objects[core.SpeciesCoefficient]), 4)

        result = core.SpeciesCoefficient.deserialize(attr, 'dna2[c]', objects)
        self.assertEqual(result[0], None)
        self.assertNotEqual(result[1], None)

        result = core.SpeciesCoefficient.deserialize(attr, 'dna', objects)
        self.assertEqual(result[0], None)
        self.assertNotEqual(result[1], None)

        result = core.SpeciesCoefficient.deserialize(
            attr, 'dna', objects, compartment=comp)
        self.assertEqual(result[0].species.species_type, dna)
        self.assertEqual(result[0].species.compartment, comp)
        self.assertEqual(result[0].coefficient, 1)
        self.assertEqual(result[1], None)

        result = core.SpeciesCoefficient.deserialize(
            attr, '(2) dna', objects, compartment=comp)
        self.assertEqual(result[0].species.species_type, dna)
        self.assertEqual(result[0].species.compartment, comp)
        self.assertEqual(result[0].coefficient, 2)
        self.assertEqual(result[1], None)

        result = core.SpeciesCoefficient.deserialize(
            attr, '(-3) dna', objects, compartment=comp)
        self.assertEqual(result[0].species.species_type, dna)
        self.assertEqual(result[0].species.compartment, comp)
        self.assertEqual(result[0].coefficient, -3)
        self.assertEqual(result[1], None)


class SpeciesTypeCoefficientTestCase(unittest.TestCase):
    def test_constructor(self):
        
        met = core.MetaboliteSpeciesType(id='met')
        species_type_coeff = core.SpeciesTypeCoefficient(species_type=met, coefficient=3)

        self.assertEqual(species_type_coeff.species_type.id, 'met')
        self.assertEqual(species_type_coeff.coefficient, 3)

    def test_serialize(self):
        
        met = core.MetaboliteSpeciesType(id='met')
        species_type_coeff = core.SpeciesTypeCoefficient(species_type=met, coefficient=3)

        self.assertEqual(species_type_coeff.serialize(), '(3) met')
        self.assertEqual(core.SpeciesTypeCoefficient._serialize(
            species_type=met, coefficient=3), '(3) met')
        self.assertEqual(core.SpeciesTypeCoefficient._serialize(
            species_type=met, coefficient=1), 'met')
        self.assertEqual(core.SpeciesTypeCoefficient._serialize(
            species_type=met, coefficient=2000), '(2.000000e+03) met')        
        
            
class SubunitAttributeTestCase(unittest.TestCase):
    def test_SubunitAttribute(self):
        compart1 = core.Compartment(id='c')
        
        met1 = core.MetaboliteSpeciesType(id='met1')
        dna1 = core.DnaSpeciesType(id='dna1')
        complex1 = core.ComplexSpeciesType(id='complex1')

        species1 = core.Species(species_type=complex1, compartment=compart1)
        
        species_type_coeff1 = core.SpeciesTypeCoefficient(
            species_type=met1, coefficient=2)
        species_type_coeff2 = core.SpeciesTypeCoefficient(
            species_type=dna1, coefficient=3)
        species_type_coeff3 = core.SpeciesTypeCoefficient(
            species_type=complex1, coefficient=5)
        
        self.assertEqual(
            core.SubunitAttribute().serialize(subunits=[]), '')
        self.assertEqual(core.SubunitAttribute().serialize(subunits=[species_type_coeff1]),
                         '(2) met1')
        self.assertEqual(core.SubunitAttribute().serialize(subunits=[species_type_coeff1, species_type_coeff2]),
                         '(3) dna1 + (2) met1')
        self.assertEqual(core.SubunitAttribute().serialize(subunits=[species_type_coeff1, species_type_coeff2, species_type_coeff3]),
                         '(5) complex1 + (3) dna1 + (2) met1')
        
        objects = {
            core.DnaSpeciesType: {
                'dna1': dna1
            },            
            core.MetaboliteSpeciesType: {
                'met1': met1
            },
            core.ComplexSpeciesType: {
                'complex1': complex1
            },
        }

        result = core.SubunitAttribute().deserialize(
            value='met1 + (2) dna1', objects=objects)
        self.assertEqual(result[0][0].species_type, met1)
        self.assertEqual(result[0][1].species_type, dna1)
        self.assertEqual(result[0][0].coefficient, 1)
        self.assertEqual(result[0][1].coefficient, 2)
        self.assertEqual(result[0][0].species_type.id, 'met1')
        self.assertEqual(result[0][1].species_type.id, 'dna1')
        self.assertEqual(result[1], None)

        result = core.SubunitAttribute().deserialize(
            value='(2) met1 + (7) complex1', objects=objects)
        self.assertEqual(result[0][0].species_type, met1)
        self.assertEqual(result[0][1].species_type, complex1)
        self.assertEqual(result[0][0].coefficient, 2)
        self.assertEqual(result[0][1].coefficient, 7)
        self.assertEqual(result[0][0].species_type.id, 'met1')
        self.assertEqual(result[0][1].species_type.id, 'complex1')
        self.assertEqual(result[1], None)

        result = core.SubunitAttribute().deserialize(
            value='[c]met1 + (2) dna1', objects=objects)
        self.assertEqual(result[0], None)
        self.assertEqual(
            result[1].messages[0], 'Incorrectly formatted participants: [c]met1 + (2) dna1')


class ReactionParticipantAttributeTestCase(unittest.TestCase):
    def test_ReactionParticipantAttribute(self):
        compart1 = core.Compartment(id='c')
        compart2 = core.Compartment(id='m')

        met1 = core.MetaboliteSpeciesType(id='met1')
        met2 = core.MetaboliteSpeciesType(id='met2')
        complex1 = core.ComplexSpeciesType(id='complex1')

        species1 = core.Species(species_type=met1, compartment=compart1)
        species2 = core.Species(species_type=met2, compartment=compart1)
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
                         '[c]: (2) met1 ==> ')
        self.assertEqual(core.ReactionParticipantAttribute().serialize(participants=[species_coeff1, species_coeff2]),
                         '[c]: (2) met1 + (3) met2 ==> ')
        self.assertEqual(core.ReactionParticipantAttribute().serialize(participants=[species_coeff1, species_coeff2, species_coeff3]),
                         '[c]: (2) met1 + (3) met2 ==> (5) complex1')
        self.assertEqual(core.ReactionParticipantAttribute().serialize(participants=[species_coeff1, species_coeff2, species_coeff4]),
                         '(2) met1[c] + (3) met2[c] ==> (7) complex1[m]')

        objects = {
            core.DnaSpeciesType: {},            
            core.MetaboliteSpeciesType: {
                'met1': met1, 'met2': met2
            },
            core.ComplexSpeciesType: {
                'complex1': complex1
            },
            core.Compartment: {
                'c': compart1, 'm': compart2
            },            
            core.Species: {
                'met1[c]': species1, 'met2[c]': species2, 'complex1[c]': species3, 'complex[m]': species4

            },
        }

        result = core.ReactionParticipantAttribute().deserialize(
            value='[c]: met1 ==> met2', objects=objects)
        self.assertEqual(result[0][0].species.species_type, met1)
        self.assertEqual(result[0][1].species.species_type, met2)
        self.assertEqual(result[0][0].coefficient, -1)
        self.assertEqual(result[0][1].coefficient,  1)
        self.assertEqual(result[0][0].species.compartment, compart1)
        self.assertEqual(result[0][1].species.compartment, compart1)
        self.assertEqual(result[0][0].species.id(), 'met1[c]')
        self.assertEqual(result[0][1].species.id(), 'met2[c]')
        self.assertEqual(result[1], None)

        result = core.ReactionParticipantAttribute().deserialize(
            value='(2) met1[c] + (3) met2[c] ==> (7) complex1[m]', objects=objects)
        self.assertEqual(result[0][0].species.species_type, met1)
        self.assertEqual(result[0][1].species.species_type, met2)
        self.assertEqual(result[0][2].species.species_type, complex1)
        self.assertEqual(result[0][0].coefficient, -2)
        self.assertEqual(result[0][1].coefficient, -3)
        self.assertEqual(result[0][2].coefficient, 7.0)
        self.assertEqual(result[0][0].species.compartment, compart1)
        self.assertEqual(result[0][1].species.compartment, compart1)
        self.assertEqual(result[0][2].species.compartment, compart2)
        self.assertEqual(result[0][0].species.id(), 'met1[c]')
        self.assertEqual(result[0][1].species.id(), 'met2[c]')
        self.assertEqual(result[0][2].species.id(), 'complex1[m]')
        self.assertEqual(result[1], None)

        result = core.ReactionParticipantAttribute().deserialize(
            value='met2[c] ==>', objects=objects)
        self.assertEqual(result[0], None)
        self.assertEqual(
            result[1].messages[0], 'Incorrectly formatted participants: met2[c] ==>')

        result = core.ReactionParticipantAttribute().deserialize(
            value='==> met1[c]', objects=objects)
        self.assertEqual(result[0], None)
        self.assertEqual(
            result[1].messages[0], 'Incorrectly formatted participants: ==> met1[c]')


class ObservableCoefficientTestCase(unittest.TestCase):
    def test_observable_coefficient(self):
        cell = core.Cell()
        comp1 = core.Compartment(id='c')
        met1 = core.MetaboliteSpeciesType(id='met1')
        dna1 = core.DnaSpeciesType(id="dna1")

        species1 = core.Species(species_type=met1, compartment=comp1)
        species2 = core.Species(species_type=dna1, compartment=comp1)

        speciesCoefficient1 = core.SpeciesCoefficient(
            species=species1, coefficient=2.5)
        speciesCoefficient2 = core.SpeciesCoefficient(
            species=species2, coefficient=3)

        observable1 = core.Observable(
            id='test', cell=cell, species=[speciesCoefficient1, speciesCoefficient2])
        observableCoefficient1 = core.ObservableCoefficient(
            observable=observable1, coefficient=2.3)

        self.assertIsInstance(observableCoefficient1,
                              core.ObservableCoefficient)
        self.assertIsInstance(
            observableCoefficient1.observable, core.Observable)
        self.assertIsInstance(observableCoefficient1.coefficient, float)
        self.assertEqual(observableCoefficient1.observable.id, 'test')


class ObservableTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def test_observables(self):
        cell = core.Cell()
        comp1 = core.Compartment(id='c')
        met1 = core.MetaboliteSpeciesType(id='met1')
        dna1 = core.DnaSpeciesType(id="dna1")

        species1 = core.Species(species_type=met1, compartment=comp1)
        species2 = core.Species(species_type=dna1, compartment=comp1)

        speciesCoefficient1 = core.SpeciesCoefficient(
            species=species1, coefficient=2)
        speciesCoefficient2 = core.SpeciesCoefficient(
            species=species2, coefficient=3.3)

        observable1 = core.Observable(
            cell=cell, species=[speciesCoefficient1, speciesCoefficient2])
        observableCoefficient1 = core.ObservableCoefficient(
            observable=observable1, coefficient=2)

        observable2 = core.Observable(cell=cell, species=[
                                      speciesCoefficient1, speciesCoefficient2], observables=[observableCoefficient1])

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
        self.assertEqual(observable1.species[0].species.id(), 'met1[c]')
        self.assertEqual(
            observable1.species[1].species.species_type.id, 'dna1')
        self.assertIsInstance(
            observable2.observables[0].observable, core.Observable)
