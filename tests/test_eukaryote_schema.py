""" Tests of the knowledge base schema for prokaryotes

:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2018-09-18
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_kb import core, eukaryote
from wc_utils.util import chem
import Bio.Alphabet
import Bio.Seq
import Bio.SeqIO
import Bio.SeqUtils
import mendeleev
import os
import shutil
import tempfile
import unittest


class CellTestCase(unittest.TestCase):
    def test_constructor(self):
        cell = core.Cell(taxon=9606)

        self.assertEqual(cell.knowledge_base, None)
        self.assertEqual(cell.taxon, 9606)
        self.assertEqual(cell.observables, [])
        self.assertEqual(cell.species_types, [])
        self.assertEqual(cell.compartments, [])
        self.assertEqual(cell.reactions, [])
        self.assertEqual(cell.loci, [])


class GenericLocusTestCase(unittest.TestCase):
    def test_serialize(self):
        gen_locus1 = eukaryote.GenericLocus(start=102, end=137)
        self.assertEqual(gen_locus1.serialize(), '102:137')


class LocusAttributeTestCase(unittest.TestCase):
    def test_LocusAttribute(self):
        gen_locus1 = eukaryote.GenericLocus(start=102, end=137)
        gen_locus2 = eukaryote.GenericLocus(start=285, end=379)

        self.assertEqual(eukaryote.LocusAttribute().serialize(
            coordinates=[gen_locus1, gen_locus2]), '102:137, 285:379')

        objects = {
            eukaryote.GenericLocus:
                {
                '102:137': gen_locus1,
                '285:379': gen_locus2,
                }
            }

        result = eukaryote.LocusAttribute().deserialize(
            value='102:137, 285:379', objects=objects)
        self.assertEqual(result[0][0].start, 102)
        self.assertEqual(result[0][0].end, 137)
        self.assertEqual(result[0][1].start, 285)
        self.assertEqual(result[0][1].end, 379)
        self.assertEqual(result[1], None)


class TranscriptSpeciesTypeTestCase(unittest.TestCase):

    def setUp(self):
        self.tmp_dirname = tempfile.mkdtemp()
        self.sequence_path = os.path.join(self.tmp_dirname, 'test_seq.fasta')
        with open(self.sequence_path, 'w') as f:
            f.write('>dna1\nACGTACGTACGTACGTTTT\n'
                    '>dna2\nAcTGAGTTACGTACGTTTT\n'
                    '>dna3\nACGT\n'
                    '>dna4\nATAT\n'
                    '>dna5\nAAAA\n'
                    '>dna6\nAACCGGTT\n')

    def tearDown(self):
        shutil.rmtree(self.tmp_dirname)

    def test_constructor(self):
        dna1 = core.DnaSpeciesType(id='dna1', sequence_path=self.sequence_path)
        gene1 = eukaryote.GeneLocus(polymer=dna1, start=1, end=15)

        transcript1 = eukaryote.TranscriptSpeciesType(
        	id='t1', name='transcript1', gene=gene1, type=eukaryote.TranscriptType.mRna)

        self.assertEqual(transcript1.id, 't1')
        self.assertEqual(transcript1.name, 'transcript1')
        self.assertEqual(transcript1.gene, gene1)
        self.assertEqual(transcript1.type.name, 'mRna')
        self.assertEqual(transcript1.exons, [])
        self.assertEqual(transcript1.comments, '')
        self.assertEqual(transcript1.references, [])
        self.assertEqual(transcript1.identifiers, [])

        exon1 = eukaryote.GenericLocus(start=1, end=1)
        exon2 = eukaryote.GenericLocus(start=2, end=2)
        transcript1.exons = [exon1, exon2]
        transcript2 = eukaryote.TranscriptSpeciesType(
            id='t2', name='transcript2', gene=gene1, exons=[exon2])

        self.assertEqual(transcript1.exons, [exon1, exon2])
        self.assertEqual(transcript2.exons, [exon2])

    def test_get_seq(self):
        dna1 = core.DnaSpeciesType(id='dna2', sequence_path=self.sequence_path)

        gene1 = eukaryote.GeneLocus(
            polymer=dna1, start=1, end=15, strand=core.PolymerStrand.positive)
        exon1_1 = eukaryote.GenericLocus(start=1, end=4)
        exon1_2 = eukaryote.GenericLocus(start=7, end=8)
        transcript1 = eukaryote.TranscriptSpeciesType(
            gene=gene1, exons=[exon1_1, exon1_2])

        gene2 = eukaryote.GeneLocus(
            polymer=dna1, start=4, end=18, strand=core.PolymerStrand.negative)
        exon2_1 = eukaryote.GenericLocus(start=4, end=10)
        exon2_2 = eukaryote.GenericLocus(start=14, end=16)
        transcript2 = eukaryote.TranscriptSpeciesType(
            gene=gene2, exons=[exon2_1, exon2_2])

        self.assertEqual(transcript1.get_seq(), 'AcUGUU')
        self.assertEqual(transcript2.get_seq(), 'ACGGUAACUC')

    def test_get_empirical_formula(self):
        dna1 = core.DnaSpeciesType(id='dna3', sequence_path=self.sequence_path)
        gene1 = eukaryote.GeneLocus(polymer=dna1, start=1, end=4)

        exon1 = eukaryote.GenericLocus(start=1, end=1)
        transcript1 = eukaryote.TranscriptSpeciesType(gene=gene1, exons=[exon1])
        self.assertEqual(transcript1.get_empirical_formula(),
                         chem.EmpiricalFormula('C10H12N5O7P'))

        exon2 = eukaryote.GenericLocus(start=2, end=2)
        transcript2 = eukaryote.TranscriptSpeciesType(gene=gene1, exons=[exon2])
        self.assertEqual(transcript2.get_empirical_formula(),
                         chem.EmpiricalFormula('C9H12N3O8P'))

        exon3 = eukaryote.GenericLocus(start=3, end=3)
        transcript3 = eukaryote.TranscriptSpeciesType(gene=gene1, exons=[exon3])
        self.assertEqual(transcript3.get_empirical_formula(),
                         chem.EmpiricalFormula('C10H12N5O8P'))

        exon4 = eukaryote.GenericLocus(start=4, end=4)
        transcript4 = eukaryote.TranscriptSpeciesType(gene=gene1, exons=[exon4])
        self.assertEqual(transcript4.get_empirical_formula(),
                         chem.EmpiricalFormula('C9H11N2O9P'))

        dna2 = core.DnaSpeciesType(id='dna4', sequence_path=self.sequence_path)
        gene2 = eukaryote.GeneLocus(polymer=dna2, start=1, end=4)
        exon5_1 = eukaryote.GenericLocus(start=1, end=1)
        exon5_2 = eukaryote.GenericLocus(start=3, end=3)
        transcript5 = eukaryote.TranscriptSpeciesType(gene=gene2, exons=[exon5_1, exon5_2])
        self.assertEqual(transcript5.get_empirical_formula(),
                         chem.EmpiricalFormula('C20H23N10O13P2'))

        # Test using input sequence
        test_trans = eukaryote.TranscriptSpeciesType()
        self.assertEqual(test_trans.get_empirical_formula(seq_input=Bio.Seq.Seq('AA')),
                         chem.EmpiricalFormula('C20H23N10O13P2'))

    def test_get_charge(self):
        dna1 = core.DnaSpeciesType(id='dna5', sequence_path=self.sequence_path)

        gene1 = eukaryote.GeneLocus(polymer=dna1, start=1, end=1)
        exon1 = eukaryote.GenericLocus(start=1, end=1)
        transcript1 = eukaryote.TranscriptSpeciesType(gene=gene1, exons=[exon1])
        self.assertEqual(transcript1.get_charge(), -2)

        gene2 = eukaryote.GeneLocus(polymer=dna1, start=2, end=4)
        exon2_1 = eukaryote.GenericLocus(start=2, end=2)
        exon2_2 = eukaryote.GenericLocus(start=4, end=4)
        transcript2 = eukaryote.TranscriptSpeciesType(gene=gene2, exons=[exon2_1, exon2_2])
        self.assertEqual(transcript2.get_charge(), -3)

        # Test using input sequence
        test_trans = eukaryote.TranscriptSpeciesType()
        self.assertEqual(test_trans.get_charge(seq_input=Bio.Seq.Seq('CG')), -3)

    def test_get_mol_wt(self):
        dna1 = core.DnaSpeciesType(id='dna6', sequence_path=self.sequence_path)
        gene1 = eukaryote.GeneLocus(polymer=dna1, start=1, end=6)

        exon1 = eukaryote.GenericLocus(start=1, end=1)
        transcript1 = eukaryote.TranscriptSpeciesType(gene=gene1, exons=[exon1])
        exp_mol_wt = \
            + Bio.SeqUtils.molecular_weight(transcript1.get_seq()) \
            - (transcript1.get_len() + 1) * mendeleev.element('H').atomic_weight
        self.assertAlmostEqual(transcript1.get_mol_wt(), exp_mol_wt, places=1)

        exon2 = eukaryote.GenericLocus(start=3, end=3)
        transcript2 = eukaryote.TranscriptSpeciesType(gene=gene1, exons=[exon2])
        exp_mol_wt = \
            + Bio.SeqUtils.molecular_weight(transcript2.get_seq()) \
            - (transcript2.get_len() + 1) * mendeleev.element('H').atomic_weight
        self.assertAlmostEqual(transcript2.get_mol_wt(), exp_mol_wt, places=1)

        exon3 = eukaryote.GenericLocus(start=5, end=5)
        transcript3 = eukaryote.TranscriptSpeciesType(gene=gene1, exons=[exon3])
        exp_mol_wt = \
            + Bio.SeqUtils.molecular_weight(transcript3.get_seq()) \
            - (transcript3.get_len() + 1) * mendeleev.element('H').atomic_weight
        self.assertAlmostEqual(transcript3.get_mol_wt(), exp_mol_wt, places=1)

        # Test using input sequence
        test_trans = eukaryote.TranscriptSpeciesType()
        exp_mol_wt = \
            + Bio.SeqUtils.molecular_weight(transcript1.get_seq()) \
            - (transcript1.get_len() + 1) * mendeleev.element('H').atomic_weight
        self.assertAlmostEqual(test_trans.get_mol_wt(seq_input=Bio.Seq.Seq('A')), exp_mol_wt, places=1)


class ProteinSpeciesTypeTestCase(unittest.TestCase):

    def setUp(self):
        self.tmp_dirname = tempfile.mkdtemp()
        sequence_path = os.path.join(self.tmp_dirname, 'test_seq.fasta')
        with open(sequence_path, 'w') as f:
            f.write('>dna1\nTTTATGAARGTNCTCATHAAYAARAAYGARCTCTAGTTTATGAARTTYAARTTYCTCCTCACNCCNCTCTAATTT\n')

        dna1 = core.DnaSpeciesType(id='dna1', sequence_path=sequence_path)

        cell1 = dna1.cell = core.Cell()

        gene1 = eukaryote.GeneLocus(polymer=dna1, start=1, end=36)
        exon1 = eukaryote.GenericLocus(start=4, end=36)
        transcript1 = eukaryote.TranscriptSpeciesType(gene=gene1, exons=[exon1])
        cds1 = eukaryote.GenericLocus(start=4, end=36)
        self.prot1 = eukaryote.ProteinSpeciesType(id='prot1', name='protein1',
            uniprot='Q12X34', transcript=transcript1, coding_regions=[cds1])

        gene2 = eukaryote.GeneLocus(polymer=dna1,
            start=30, end=75, strand=core.PolymerStrand.positive)
        exon2_1 = eukaryote.GenericLocus(start=32, end=35)
        exon2_2 = eukaryote.GenericLocus(start=38, end=45)
        exon2_3 = eukaryote.GenericLocus(start=49, end=54)
        exon2_4 = eukaryote.GenericLocus(start=55, end=72)
        exon2_5 = eukaryote.GenericLocus(start=73, end=74)
        transcript2 = eukaryote.TranscriptSpeciesType(
            gene=gene2, exons=[exon2_1, exon2_2, exon2_3, exon2_4, exon2_5])
        cds2_2 = eukaryote.GenericLocus(start=40, end=45)
        cds2_3 = eukaryote.GenericLocus(start=49, end=54)
        cds2_4 = eukaryote.GenericLocus(start=55, end=72)
        self.prot2 = eukaryote.ProteinSpeciesType(id='prot2', name='protein2',
            uniprot='P12345', cell=cell1, transcript=transcript2,
            coding_regions=[cds2_4, cds2_2, cds2_3])

    def tearDown(self):
        shutil.rmtree(self.tmp_dirname)

    def test_constructor(self):

        self.assertEqual(self.prot1.id, 'prot1')
        self.assertEqual(self.prot1.name, 'protein1')
        self.assertEqual(self.prot1.uniprot, 'Q12X34')
        self.assertEqual(self.prot1.coding_regions[0].serialize(), '4:36')
        self.assertEqual(self.prot1.comments, '')
        self.assertEqual(self.prot1.references, [])
        self.assertEqual(self.prot1.identifiers, [])
        self.assertEqual(self.prot1.cell, None)

    def test_get_seq(self):
        # Default translation table used is 1 (standard)
        self.assertEqual(self.prot1.get_seq(), 'MKVLINKNEL')
        self.assertEqual(self.prot2.get_seq(), 'MKKFLLTPL')

    def test_get_seq_and_start_codon(self):
        # Default translation table used is 1 (standard)
        coding_rna_seq, aa_seq, start_codon = self.prot1.get_seq_and_start_codon()
        self.assertEqual(coding_rna_seq, 'AUGAARGUNCUCAUHAAYAARAAYGARCUCUAG')
        self.assertEqual(aa_seq, 'MKVLINKNEL')
        self.assertEqual(start_codon, 'AUG')

        coding_rna_seq, aa_seq, start_codon = self.prot2.get_seq_and_start_codon()
        self.assertEqual(coding_rna_seq, 'AUGAARAARUUYCUCCUCACNCCNCUCUAA')
        self.assertEqual(aa_seq, 'MKKFLLTPL')
        self.assertEqual(start_codon, 'AUG')
        
    def test_get_empirical_formula(self):
        # Default translation table used is 1 (standard)
        self.assertEqual(self.prot1.get_empirical_formula(),
                         chem.EmpiricalFormula('C53H96N14O15S1'))
        self.assertEqual(self.prot2.get_empirical_formula(),
                         chem.EmpiricalFormula('C53H91N11O11S1'))

        # Test using input sequence
        test_prot = eukaryote.ProteinSpeciesType()
        self.assertEqual(test_prot.get_empirical_formula(seq_input=Bio.Seq.Seq('MKVLINKNEL')),
                         chem.EmpiricalFormula('C53H96N14O15S1'))
        self.assertEqual(test_prot.get_empirical_formula(seq_input=Bio.Seq.Seq('MKKFLLTPL')),
                         chem.EmpiricalFormula('C53H91N11O11S1'))

    def test_get_mol_wt(self):
        # Default translation table used is 1 (standard)
        self.assertAlmostEqual(self.prot1.get_mol_wt(), 1201.49, delta=0.3)
        self.assertAlmostEqual(self.prot2.get_mol_wt(), 1090.43, delta=0.3)

        # Test using input sequence
        test_prot = eukaryote.ProteinSpeciesType()
        self.assertAlmostEqual(test_prot.get_mol_wt(seq_input=Bio.Seq.Seq('MKVLINKNEL')), 1201.49, delta=0.3)
        self.assertAlmostEqual(test_prot.get_mol_wt(seq_input=Bio.Seq.Seq('MKKFLLTPL')), 1090.43, delta=0.3)

    def test_get_charge(self):
        # Default translation table used is 1 (standard)
        self.assertEqual(self.prot1.get_charge(), 1)
        self.assertEqual(self.prot2.get_charge(), 2)

        # Test using input sequence
        test_prot = eukaryote.ProteinSpeciesType()
        self.assertEqual(test_prot.get_charge(seq_input=Bio.Seq.Seq('MKVLINKNEL')), 1)
        self.assertEqual(test_prot.get_charge(seq_input=Bio.Seq.Seq('MKKFLLTPL')), 2)

class ComplexSpeciesTypeTestCase(unittest.TestCase):
    def test_ComplexSpeciesType(self):

        self.tmp_dirname = tempfile.mkdtemp()
        sequence_path = os.path.join(self.tmp_dirname, 'test_seq.fasta')
        with open(sequence_path, 'w') as f:
            f.write('>dna1\nTTTATGAARGTNCTCATHAAYAARAAYGARCTCTAGTTTATGAARTTYAARTTYCTCCTCACNCCNCTCTAATTT\n')

        dna1 = core.DnaSpeciesType(id='dna1', sequence_path=sequence_path)

        # Protein subunit 1
        gene1 = eukaryote.GeneLocus(polymer=dna1, start=1, end=36)
        exon1 = eukaryote.GenericLocus(start=4, end=36)
        transcript1 = eukaryote.TranscriptSpeciesType(gene=gene1, exons=[exon1])
        cds1 = eukaryote.GenericLocus(start=4, end=36)
        prot1 = eukaryote.ProteinSpeciesType(transcript=transcript1, coding_regions=[cds1])

        # Protein subunit 2
        gene2 = eukaryote.GeneLocus(polymer=dna1, start=37, end=75)
        exon2 = eukaryote.GenericLocus(start=40, end=72)
        transcript2 = eukaryote.TranscriptSpeciesType(gene=gene2, exons=[exon2])
        cds2 = eukaryote.GenericLocus(start=40, end=72)
        prot2 = eukaryote.ProteinSpeciesType(transcript=transcript2, coding_regions=[cds2])

        species_coeff1 = core.SpeciesTypeCoefficient(
            species_type=prot1, coefficient=2)
        species_coeff2 = core.SpeciesTypeCoefficient(
            species_type=prot2, coefficient=3)
        complex1 = core.ComplexSpeciesType(
            subunits = [species_coeff1, species_coeff2])

        self.assertEqual(complex1.get_charge(), 8)
        self.assertAlmostEqual(complex1.get_mol_wt(),
                               (2*prot1.get_mol_wt() + 3 * prot2.get_mol_wt()))
        self.assertEqual(complex1.get_empirical_formula(),
                         chem.EmpiricalFormula('C292H492N64O66S5'))

        shutil.rmtree(self.tmp_dirname)


class GeneLocusTestCase(unittest.TestCase):
    def test_constructor(self):
        gene = eukaryote.GeneLocus(id='gene1', name='gene1', symbol='gene_1',
            strand=core.PolymerStrand.negative, start=1, end=2)
        self.assertEqual(gene.id, 'gene1')
        self.assertEqual(gene.name, 'gene1')
        self.assertEqual(gene.symbol, 'gene_1')
        self.assertEqual(gene.strand, core.PolymerStrand.negative)
        self.assertEqual(gene.start, 1)
        self.assertEqual(gene.end, 2)
        self.assertEqual(gene.comments, '')
        self.assertEqual(gene.references, [])
        self.assertEqual(gene.identifiers, [])


class TranscriptionFactorRegulationTestCase(unittest.TestCase):
    def test_constructor(self):

        tf1 = eukaryote.ProteinSpeciesType(id='tf1')
        tf2 = eukaryote.ProteinSpeciesType(id='tf2')
        tf_reg1 = eukaryote.TranscriptionFactorRegulation(
                    transcription_factor=tf1,
                    direction=eukaryote.RegulatoryDirection.repression)
        tf_reg2 = eukaryote.TranscriptionFactorRegulation(
                    transcription_factor=tf2,
                    direction=eukaryote.RegulatoryDirection.activation)
        tf_reg3 = eukaryote.TranscriptionFactorRegulation(
                    transcription_factor=tf2,
                    direction=eukaryote.RegulatoryDirection.repression)

        self.assertEqual(tf_reg1.transcription_factor, tf1)
        self.assertEqual(tf_reg1.direction.name, 'repression')
        self.assertEqual(tf_reg2.direction.name, 'activation')
        self.assertEqual(tf_reg3.transcription_factor, tf2)

    def test_serialize(self):

        tf1 = eukaryote.ProteinSpeciesType(id='tf1')
        tf_reg1 = eukaryote.TranscriptionFactorRegulation(
                    transcription_factor=tf1,
                    direction=eukaryote.RegulatoryDirection.repression)
        tf_reg2 = eukaryote.TranscriptionFactorRegulation(
                    transcription_factor=tf1,
                    direction=eukaryote.RegulatoryDirection.activation)

        self.assertEqual(tf_reg1.serialize(), 'tf1:repression')
        self.assertEqual(tf_reg2.serialize(), 'tf1:activation')

    def test_deserialize(self):
        
        tf1 = eukaryote.ProteinSpeciesType(id='tf1')
        
        objects = {
            eukaryote.ProteinSpeciesType: {
                'tf1': tf1,
            },
        }    
                
        result = eukaryote.TranscriptionFactorRegulation.deserialize('tf1:activation', objects)
        self.assertEqual(result[0].transcription_factor, tf1)
        self.assertEqual(result[0].direction, eukaryote.RegulatoryDirection.activation)
        self.assertEqual(result[1], None)

        result = eukaryote.TranscriptionFactorRegulation.deserialize('tf1:repression', objects)
        self.assertEqual(result[0].transcription_factor, tf1)
        self.assertEqual(result[0].direction, eukaryote.RegulatoryDirection.repression)
        self.assertEqual(result[1], None)

        result = eukaryote.TranscriptionFactorRegulation.deserialize('tf1:unknown', objects)
        self.assertEqual(result[0].transcription_factor, tf1)
        self.assertEqual(result[0].direction, eukaryote.RegulatoryDirection.unknown)
        self.assertEqual(result[1], None)

        result = eukaryote.TranscriptionFactorRegulation.deserialize('tf:activation', objects)
        self.assertEqual(result[0], None)

        result = eukaryote.TranscriptionFactorRegulation.deserialize('tf1:1', objects)
        self.assertEqual(result[0], None)

        result = eukaryote.TranscriptionFactorRegulation.deserialize('tf1:error', objects)
        self.assertEqual(result[0], None)

        result = eukaryote.TranscriptionFactorRegulation.deserialize('tf1:3.6', objects)
        self.assertEqual(result[0], None)


class RegulatoryModuleTestCase(unittest.TestCase):
    def test_constructor(self):
        dna1 = core.DnaSpeciesType(circular=False, double_stranded=False)

        gene1 = eukaryote.GeneLocus(polymer=dna1, start=9, end=15)
        gene2 = eukaryote.GeneLocus(polymer=dna1, start=17, end=18)

        promoter1 = 'ENSR00000172399'
        promoter2 = 'ENSR00000309980'

        tf1 = eukaryote.ProteinSpeciesType(id='tf1')
        tf2 = eukaryote.ProteinSpeciesType(id='tf2')

        reg_module1 = eukaryote.RegulatoryModule(
            gene=gene1,
            promoter=promoter1,
            activity=eukaryote.ActivityLevel.active,
            type=eukaryote.RegulationType.proximal,
            transcription_factor_regulation=[
                eukaryote.TranscriptionFactorRegulation(
                    transcription_factor=tf1,
                    direction=eukaryote.RegulatoryDirection.activation),
                eukaryote.TranscriptionFactorRegulation(
                    transcription_factor=tf2,
                    direction=eukaryote.RegulatoryDirection.repression)
            ]) 

        reg_module2 = eukaryote.RegulatoryModule(
            gene=gene1,
            promoter=promoter1,
            activity=eukaryote.ActivityLevel.active,
            type=eukaryote.RegulationType.distal,
            transcription_factor_regulation=[eukaryote.TranscriptionFactorRegulation(
                transcription_factor=tf1,
                direction=eukaryote.RegulatoryDirection.repression)]) 

        reg_module3 = eukaryote.RegulatoryModule(
            id='rm3',
            name='reg_module3',
            activity=eukaryote.ActivityLevel.inactive,
            gene=gene2,
            promoter=promoter2)

        self.assertEqual(reg_module1.gene, gene1)
        self.assertEqual(reg_module1.promoter, promoter1)
        self.assertEqual(reg_module1.activity.name, 'active')
        self.assertEqual(reg_module1.type.value, 1)
        self.assertEqual(sorted([i.transcription_factor.id for i in reg_module1.transcription_factor_regulation]), 
            ['tf1', 'tf2'])
        self.assertEqual(sorted([i.direction.name for i in reg_module1.transcription_factor_regulation]), 
            ['activation', 'repression'])
        self.assertEqual(reg_module2.gene, gene1)
        self.assertEqual(reg_module2.promoter, promoter1)
        self.assertEqual(reg_module2.activity.name, 'active')
        self.assertEqual(reg_module2.type.value, 2)
        self.assertEqual(reg_module2.transcription_factor_regulation[0].transcription_factor, tf1)
        self.assertEqual(reg_module2.transcription_factor_regulation[0].direction.name, 'repression')
        self.assertEqual(reg_module3.id, 'rm3')
        self.assertEqual(reg_module3.activity.name, 'inactive')
        self.assertEqual(reg_module3.name, 'reg_module3')
        self.assertEqual(reg_module3.promoter, promoter2)


class PtmSiteTestCase(unittest.TestCase):
    def test_constructor(self):
        # defining modified protein name
        mp = eukaryote.ProteinSpeciesType(id='mp')

        # testing example of modification: one protein one modified site
        ptm1 = eukaryote.PtmSite(id='protptm1', name='ptm1', modified_protein=mp,
             type='phosphorylation', modified_residue='s145', fractional_abundance='0.5', comments='oneprot_onesite')

        self.assertEqual(ptm1.id, 'protptm1')
        self.assertEqual(ptm1.name, 'ptm1')
        self.assertEqual(ptm1.modified_protein, mp)
        self.assertEqual(ptm1.type, 'phosphorylation')
        self.assertEqual(ptm1.modified_residue, 's145')
        self.assertEqual(ptm1.fractional_abundance, '0.5')
        self.assertEqual(ptm1.comments, 'oneprot_onesite')
        self.assertEqual(ptm1.references, [])
        self.assertEqual(ptm1.identifiers, [])
