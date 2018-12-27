""" Tests of the knowledge base schema for prokaryotes

:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2018-09-18
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_kb import core, eukaryote_schema
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


class TestCore(unittest.TestCase):
    
    def test_RegulatoryElementType(self):
        self.assertEqual(eukaryote_schema.RegulatoryElementType.promoter.value, 1)
        self.assertEqual(eukaryote_schema.RegulatoryElementType.promoter_flanking_region.value, 2)
        self.assertEqual(eukaryote_schema.RegulatoryElementType.enhancer.value, 3)
        self.assertEqual(eukaryote_schema.RegulatoryElementType.CTCF_binding_site.value, 4)
        self.assertEqual(eukaryote_schema.RegulatoryElementType.TF_binding_site.value, 5)
        self.assertEqual(eukaryote_schema.RegulatoryElementType.open_chromatin_region.value, 6)

    def test_ActivityLevel(self):
        self.assertEqual(eukaryote_schema.ActivityLevel.active.value, 1)
        self.assertEqual(eukaryote_schema.ActivityLevel.poised.value, 2)
        self.assertEqual(eukaryote_schema.ActivityLevel.repressed.value, 3)
        self.assertEqual(eukaryote_schema.ActivityLevel.inactive.value, 4)
        self.assertEqual(eukaryote_schema.ActivityLevel.na.value, 5)


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

        self.assertEqual(cell.species_types.get(
            __type=eukaryote_schema.PreRnaSpeciesType), [])
        self.assertEqual(cell.loci.get(__type=eukaryote_schema.RegulatoryElementLocus), [])


class PreRnaSpeciesTypeTestCase(unittest.TestCase):

    def setUp(self):
        self.tmp_dirname = tempfile.mkdtemp()
        self.sequence_path = os.path.join(self.tmp_dirname, 'test_seq.fasta')
        with open(self.sequence_path, 'w') as f:
            f.write('>dna1\nACTGAGTTACGTACGTTTT\n'
                    '>dna2\nACGT\n'
                    '>dna3\nAAAA\n'
                    '>dna4\nAACCGGTT\n')

    def tearDown(self):    
        shutil.rmtree(self.tmp_dirname)  

    def test_constructor(self):
        dna1 = core.DnaSpeciesType(id='dna1')
        gene1 = eukaryote_schema.GeneLocus(polymer=dna1, start=1, end=15)
        rna1 = eukaryote_schema.PreRnaSpeciesType(id='rna1', name='rna1', gene=gene1, 
        	type=1, half_life=2)

        self.assertEqual(rna1.id, 'rna1')
        self.assertEqual(rna1.name, 'rna1')
        self.assertEqual(rna1.gene, gene1)
        self.assertEqual(rna1.type, 1)
        self.assertEqual(rna1.half_life, 2)
        self.assertEqual(rna1.comments, '')
        self.assertEqual(rna1.references, [])
        self.assertEqual(rna1.database_references, [])

    def test_get_seq(self):
        dna1 = core.DnaSpeciesType(id='dna1', sequence_path=self.sequence_path)
        gene1 = eukaryote_schema.GeneLocus(
            polymer=dna1, start=1, end=15, strand=core.PolymerStrand.positive)
        rna1 = eukaryote_schema.PreRnaSpeciesType(gene=gene1)

        gene2 = eukaryote_schema.GeneLocus(
            polymer=dna1, start=4, end=18, strand=core.PolymerStrand.negative)
        rna2 = eukaryote_schema.PreRnaSpeciesType(gene=gene2)            

        self.assertEqual(rna1.get_seq(), 'ACUGAGUUACGUACG')
        self.assertEqual(rna2.get_seq(), 'AAACGUACGUAACUC')                

    def test_get_empirical_formula(self):
        dna1 = core.DnaSpeciesType(id='dna2', sequence_path=self.sequence_path)
        
        gene1 = eukaryote_schema.GeneLocus(polymer=dna1, start=1, end=1)
        rna1 = eukaryote_schema.PreRnaSpeciesType(gene=gene1)
        self.assertEqual(rna1.get_empirical_formula(),
                         chem.EmpiricalFormula('C10H12N5O7P'))

        gene2 = eukaryote_schema.GeneLocus(polymer=dna1, start=2, end=2)
        rna2 = eukaryote_schema.PreRnaSpeciesType(gene=gene2)
        self.assertEqual(rna2.get_empirical_formula(),
                         chem.EmpiricalFormula('C9H12N3O8P'))

        gene3 = eukaryote_schema.GeneLocus(polymer=dna1, start=3, end=3)
        rna3 = eukaryote_schema.PreRnaSpeciesType(gene=gene3)
        self.assertEqual(rna3.get_empirical_formula(),
                         chem.EmpiricalFormula('C10H12N5O8P'))

        gene4 = eukaryote_schema.GeneLocus(polymer=dna1, start=4, end=4)
        rna4 = eukaryote_schema.PreRnaSpeciesType(gene=gene4)
        self.assertEqual(rna4.get_empirical_formula(),
                         chem.EmpiricalFormula('C9H11N2O9P'))

        dna2 = core.DnaSpeciesType(id='dna3', sequence_path=self.sequence_path)
        gene5 = eukaryote_schema.GeneLocus(polymer=dna2, start=1, end=2)
        rna5 = eukaryote_schema.PreRnaSpeciesType(gene=gene5)
        self.assertEqual(rna5.get_empirical_formula(),
                         chem.EmpiricalFormula('C20H23N10O13P2'))

    def test_get_charge(self):
        dna1 = core.DnaSpeciesType(id='dna3', sequence_path=self.sequence_path)
        
        gene1 = eukaryote_schema.GeneLocus(polymer=dna1, start=1, end=1)
        rna1 = eukaryote_schema.PreRnaSpeciesType(gene=gene1)
        self.assertEqual(rna1.get_charge(), -2)

        gene2 = eukaryote_schema.GeneLocus(polymer=dna1, start=2, end=3)
        rna2 = eukaryote_schema.PreRnaSpeciesType(gene=gene2)
        self.assertEqual(rna2.get_charge(), -3)

    def test_get_mol_wt(self):
        dna1 = core.DnaSpeciesType(id='dna4', sequence_path=self.sequence_path)
        
        gene1 = eukaryote_schema.GeneLocus(polymer=dna1, start=1, end=1)
        rna1 = eukaryote_schema.PreRnaSpeciesType(gene=gene1)
        exp_mol_wt = \
            + Bio.SeqUtils.molecular_weight(rna1.get_seq()) \
            - (rna1.get_len() + 1) * mendeleev.element('H').atomic_weight
        self.assertAlmostEqual(rna1.get_mol_wt(), exp_mol_wt, places=1)

        gene2 = eukaryote_schema.GeneLocus(polymer=dna1, start=3, end=3)
        rna2 = eukaryote_schema.PreRnaSpeciesType(gene=gene2)
        exp_mol_wt = \
            + Bio.SeqUtils.molecular_weight(rna2.get_seq()) \
            - (rna2.get_len() + 1) * mendeleev.element('H').atomic_weight
        self.assertAlmostEqual(rna2.get_mol_wt(), exp_mol_wt, places=1)

        gene3 = eukaryote_schema.GeneLocus(polymer=dna1, start=5, end=5)
        rna3 = eukaryote_schema.PreRnaSpeciesType(gene=gene3)
        exp_mol_wt = \
            + Bio.SeqUtils.molecular_weight(rna3.get_seq()) \
            - (rna3.get_len() + 1) * mendeleev.element('H').atomic_weight
        self.assertAlmostEqual(rna3.get_mol_wt(), exp_mol_wt, places=1)


class TranscriptSpeciesTypeTestCase(unittest.TestCase):

    def setUp(self):
        self.tmp_dirname = tempfile.mkdtemp()
        self.sequence_path = os.path.join(self.tmp_dirname, 'test_seq.fasta')
        with open(self.sequence_path, 'w') as f:
            f.write('>dna1\nACGTACGTACGTACGTTTT\n'
                    '>dna2\nACTGAGTTACGTACGTTTT\n'
                    '>dna3\nACGT\n'
                    '>dna4\nATAT\n'
                    '>dna5\nAAAA\n'
                    '>dna6\nAACCGGTT\n')

    def tearDown(self):    
        shutil.rmtree(self.tmp_dirname)

    def test_constructor(self):
        dna1 = core.DnaSpeciesType(id='dna1', sequence_path=self.sequence_path)
        gene1 = eukaryote_schema.GeneLocus(polymer=dna1, start=1, end=15)
        rna1 = eukaryote_schema.PreRnaSpeciesType(gene=gene1)

        transcript1 = eukaryote_schema.TranscriptSpeciesType(
        	id='t1', name='transcript1', rna=rna1, half_life=2)

        self.assertEqual(transcript1.id, 't1')
        self.assertEqual(transcript1.name, 'transcript1')
        self.assertEqual(transcript1.rna, rna1)
        self.assertEqual(transcript1.half_life, 2)
        self.assertEqual(transcript1.exons, [])
        self.assertEqual(transcript1.comments, '')
        self.assertEqual(transcript1.references, [])
        self.assertEqual(transcript1.database_references, [])

        exon1 = eukaryote_schema.ExonLocus(start=1, end=1)
        exon2 = eukaryote_schema.ExonLocus(start=2, end=2)
        transcript1.exons = [exon1, exon2]
        transcript2 = eukaryote_schema.TranscriptSpeciesType(
            id='t2', name='transcript2', rna=rna1, exons=[exon2])

        self.assertEqual(transcript1.exons, [exon1, exon2])
        self.assertEqual(transcript2.exons, [exon2])

    def test_get_seq(self):
        dna1 = core.DnaSpeciesType(id='dna2', sequence_path=self.sequence_path)
        
        gene1 = eukaryote_schema.GeneLocus(
            polymer=dna1, start=1, end=15, strand=core.PolymerStrand.positive)
        rna1 = eukaryote_schema.PreRnaSpeciesType(gene=gene1)        
        exon1_1 = eukaryote_schema.ExonLocus(start=1, end=4)
        exon1_2 = eukaryote_schema.ExonLocus(start=7, end=8)
        transcript1 = eukaryote_schema.TranscriptSpeciesType(
            rna=rna1, exons=[exon1_1, exon1_2])
        
        gene2 = eukaryote_schema.GeneLocus(
            polymer=dna1, start=4, end=18, strand=core.PolymerStrand.negative)
        rna2 = eukaryote_schema.PreRnaSpeciesType(gene=gene2)      
        exon2_1 = eukaryote_schema.ExonLocus(start=4, end=10)
        exon2_2 = eukaryote_schema.ExonLocus(start=14, end=16)
        transcript2 = eukaryote_schema.TranscriptSpeciesType(
            rna=rna2, exons=[exon2_1, exon2_2])

        self.assertEqual(transcript1.get_seq(), 'ACUGUU')
        self.assertEqual(transcript2.get_seq(), 'ACGGUAACUC') 

    def test_get_empirical_formula(self):
        dna1 = core.DnaSpeciesType(id='dna3', sequence_path=self.sequence_path)        
        gene1 = eukaryote_schema.GeneLocus(polymer=dna1, start=1, end=4)
        rna1 = eukaryote_schema.PreRnaSpeciesType(gene=gene1)
        
        exon1 = eukaryote_schema.ExonLocus(start=1, end=1)
        transcript1 = eukaryote_schema.TranscriptSpeciesType(rna=rna1, exons=[exon1])        
        self.assertEqual(transcript1.get_empirical_formula(),
                         chem.EmpiricalFormula('C10H12N5O7P'))

        exon2 = eukaryote_schema.ExonLocus(start=2, end=2)
        transcript2 = eukaryote_schema.TranscriptSpeciesType(rna=rna1, exons=[exon2])        
        self.assertEqual(transcript2.get_empirical_formula(),
                         chem.EmpiricalFormula('C9H12N3O8P'))

        exon3 = eukaryote_schema.ExonLocus(start=3, end=3)
        transcript3 = eukaryote_schema.TranscriptSpeciesType(rna=rna1, exons=[exon3])        
        self.assertEqual(transcript3.get_empirical_formula(),
                         chem.EmpiricalFormula('C10H12N5O8P'))

        exon4 = eukaryote_schema.ExonLocus(start=4, end=4)
        transcript4 = eukaryote_schema.TranscriptSpeciesType(rna=rna1, exons=[exon4])
        self.assertEqual(transcript4.get_empirical_formula(),
                         chem.EmpiricalFormula('C9H11N2O9P'))

        dna2 = core.DnaSpeciesType(id='dna4', sequence_path=self.sequence_path)
        gene2 = eukaryote_schema.GeneLocus(polymer=dna2, start=1, end=4)
        rna2 = eukaryote_schema.PreRnaSpeciesType(gene=gene2)
        exon5_1 = eukaryote_schema.ExonLocus(start=1, end=1)
        exon5_2 = eukaryote_schema.ExonLocus(start=3, end=3)
        transcript5 = eukaryote_schema.TranscriptSpeciesType(rna=rna2, exons=[exon5_1, exon5_2])        
        self.assertEqual(transcript5.get_empirical_formula(),
                         chem.EmpiricalFormula('C20H23N10O13P2'))

    def test_get_charge(self):
        dna1 = core.DnaSpeciesType(id='dna5', sequence_path=self.sequence_path)        
        
        gene1 = eukaryote_schema.GeneLocus(polymer=dna1, start=1, end=1)
        rna1 = eukaryote_schema.PreRnaSpeciesType(gene=gene1)
        exon1 = eukaryote_schema.ExonLocus(start=1, end=1)
        transcript1 = eukaryote_schema.TranscriptSpeciesType(rna=rna1, exons=[exon1])
        self.assertEqual(transcript1.get_charge(), -2)

        gene2 = eukaryote_schema.GeneLocus(polymer=dna1, start=2, end=4)
        rna2 = eukaryote_schema.PreRnaSpeciesType(gene=gene2)
        exon2_1 = eukaryote_schema.ExonLocus(start=2, end=2)
        exon2_2 = eukaryote_schema.ExonLocus(start=4, end=4)
        transcript2 = eukaryote_schema.TranscriptSpeciesType(rna=rna2, exons=[exon2_1, exon2_2])        
        self.assertEqual(transcript2.get_charge(), -3)

    def test_get_mol_wt(self):
        dna1 = core.DnaSpeciesType(id='dna6', sequence_path=self.sequence_path)
        gene1 = eukaryote_schema.GeneLocus(polymer=dna1, start=1, end=6)
        rna1 = eukaryote_schema.PreRnaSpeciesType(gene=gene1)

        exon1 = eukaryote_schema.ExonLocus(start=1, end=1)
        transcript1 = eukaryote_schema.TranscriptSpeciesType(rna=rna1, exons=[exon1])        
        exp_mol_wt = \
            + Bio.SeqUtils.molecular_weight(transcript1.get_seq()) \
            - (transcript1.get_len() + 1) * mendeleev.element('H').atomic_weight
        self.assertAlmostEqual(transcript1.get_mol_wt(), exp_mol_wt, places=1)

        exon2 = eukaryote_schema.ExonLocus(start=3, end=3)
        transcript2 = eukaryote_schema.TranscriptSpeciesType(rna=rna1, exons=[exon2])        
        exp_mol_wt = \
            + Bio.SeqUtils.molecular_weight(transcript2.get_seq()) \
            - (transcript2.get_len() + 1) * mendeleev.element('H').atomic_weight
        self.assertAlmostEqual(transcript2.get_mol_wt(), exp_mol_wt, places=1)

        exon3 = eukaryote_schema.ExonLocus(start=5, end=5)
        transcript3 = eukaryote_schema.TranscriptSpeciesType(rna=rna1, exons=[exon3])        
        exp_mol_wt = \
            + Bio.SeqUtils.molecular_weight(transcript3.get_seq()) \
            - (transcript3.get_len() + 1) * mendeleev.element('H').atomic_weight
        self.assertAlmostEqual(transcript3.get_mol_wt(), exp_mol_wt, places=1)


class ProteinSpeciesTypeTestCase(unittest.TestCase):

    def setUp(self):
        self.tmp_dirname = tempfile.mkdtemp()
        sequence_path = os.path.join(self.tmp_dirname, 'test_seq.fasta')
        with open(sequence_path, 'w') as f:
            f.write('>dna1\nTTTATGAARGTNCTCATHAAYAARAAYGARCTCTAGTTTATGAARTTYAARTTYCTCCTCACNCCNCTCTAATTT\n')
    
        dna1 = core.DnaSpeciesType(id='dna1', sequence_path=sequence_path)

        cell1 = dna1.cell = core.Cell()           

        gene1 = eukaryote_schema.GeneLocus(polymer=dna1, start=1, end=36)
        rna1 = eukaryote_schema.PreRnaSpeciesType(gene=gene1)
        exon1 = eukaryote_schema.ExonLocus(start=4, end=36)
        transcript1 = eukaryote_schema.TranscriptSpeciesType(rna=rna1, exons=[exon1])
        cds1 = eukaryote_schema.CdsLocus(id='cds1', exon=exon1, start=4, end=36)        
        self.prot1 = eukaryote_schema.ProteinSpeciesType(id='prot1', name='protein1', 
            uniprot='Q12X34', transcript=transcript1, coding_regions=[cds1], half_life=0.35)

        gene2 = eukaryote_schema.GeneLocus(polymer=dna1,
            start=30, end=75, strand=core.PolymerStrand.positive)
        rna2 = eukaryote_schema.PreRnaSpeciesType(gene=gene2)
        exon2_1 = eukaryote_schema.ExonLocus(start=32, end=35)
        exon2_2 = eukaryote_schema.ExonLocus(start=38, end=45)
        exon2_3 = eukaryote_schema.ExonLocus(start=49, end=54)
        exon2_4 = eukaryote_schema.ExonLocus(start=55, end=72)
        exon2_5 = eukaryote_schema.ExonLocus(start=73, end=74)
        transcript2 = eukaryote_schema.TranscriptSpeciesType(
            rna=rna2, exons=[exon2_1, exon2_2, exon2_3, exon2_4, exon2_5])
        cds2_2 = eukaryote_schema.CdsLocus(id='cds2', exon=exon2_2, start=40, end=45)
        cds2_3 = eukaryote_schema.CdsLocus(id='cds2', exon=exon2_3, start=49, end=54)
        cds2_4 = eukaryote_schema.CdsLocus(id='cds2', exon=exon2_4, start=55, end=72)        
        self.prot2 = eukaryote_schema.ProteinSpeciesType(id='prot2', name='protein2', 
            uniprot='P12345', cell=cell1, transcript=transcript2, 
            coding_regions=[cds2_4, cds2_2, cds2_3])

    def tearDown(self):    
        shutil.rmtree(self.tmp_dirname)    

    def test_constructor(self):
        
        self.assertEqual(self.prot1.id, 'prot1')
        self.assertEqual(self.prot1.name, 'protein1')
        self.assertEqual(self.prot1.uniprot, 'Q12X34')
        self.assertEqual(self.prot1.coding_regions[0].id, 'cds1')
        self.assertEqual(self.prot1.half_life, 0.35)
        self.assertEqual(self.prot1.comments, '')
        self.assertEqual(self.prot1.references, [])
        self.assertEqual(self.prot1.database_references, [])
        self.assertEqual(self.prot1.cell, None)

    def test_get_seq(self):
        # Default translation table used is 1 (standard)
        self.assertEqual(self.prot1.get_seq(), 'MKVLINKNEL')
        self.assertEqual(self.prot2.get_seq(), 'MKKFLLTPL')

    def test_get_empirical_formula(self):
        # Default translation table used is 1 (standard)
        self.assertEqual(self.prot1.get_empirical_formula(),
                         chem.EmpiricalFormula('C53H96N14O15S1'))
        self.assertEqual(self.prot2.get_empirical_formula(),
                         chem.EmpiricalFormula('C53H91N11O11S1'))

    def test_get_mol_wt(self):
        # Default translation table used is 1 (standard)
        self.assertAlmostEqual(self.prot1.get_mol_wt(), 1201.49, delta=0.3)
        self.assertAlmostEqual(self.prot2.get_mol_wt(), 1090.43, delta=0.3)

    def test_get_charge(self):
        # Default translation table used is 1 (standard)
        self.assertEqual(self.prot1.get_charge(), 1)
        self.assertEqual(self.prot2.get_charge(), 2)


class ComplexSpeciesTypeTestCase(unittest.TestCase):
    def test_ComplexSpeciesType(self):

        self.tmp_dirname = tempfile.mkdtemp()
        sequence_path = os.path.join(self.tmp_dirname, 'test_seq.fasta')
        with open(sequence_path, 'w') as f:
            f.write('>dna1\nTTTATGAARGTNCTCATHAAYAARAAYGARCTCTAGTTTATGAARTTYAARTTYCTCCTCACNCCNCTCTAATTT\n')        
        
        dna1 = core.DnaSpeciesType(id='dna1', sequence_path=sequence_path)

        # Protein subunit 1
        gene1 = eukaryote_schema.GeneLocus(polymer=dna1, start=1, end=36)
        rna1 = eukaryote_schema.PreRnaSpeciesType(gene=gene1)
        exon1 = eukaryote_schema.ExonLocus(start=4, end=36)
        transcript1 = eukaryote_schema.TranscriptSpeciesType(rna=rna1, exons=[exon1])
        cds1 = eukaryote_schema.CdsLocus(exon=exon1, start=4, end=36)        
        prot1 = eukaryote_schema.ProteinSpeciesType(transcript=transcript1, coding_regions=[cds1])

        # Protein subunit 2
        gene2 = eukaryote_schema.GeneLocus(polymer=dna1, start=37, end=75)
        rna2 = eukaryote_schema.PreRnaSpeciesType(gene=gene2)
        exon2 = eukaryote_schema.ExonLocus(start=40, end=72)
        transcript2 = eukaryote_schema.TranscriptSpeciesType(rna=rna2, exons=[exon2])
        cds2 = eukaryote_schema.CdsLocus(exon=exon2, start=40, end=72)        
        prot2 = eukaryote_schema.ProteinSpeciesType(transcript=transcript2, coding_regions=[cds2])

        # Complex formation: (2) prot1 + (3) prot2 ==> complex1  
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
        gene = eukaryote_schema.GeneLocus(id='gene1', name='gene1', symbol='gene_1', 
            strand=core.PolymerStrand.negative, start=1, end=2)
        self.assertEqual(gene.id, 'gene1')
        self.assertEqual(gene.name, 'gene1')
        self.assertEqual(gene.symbol, 'gene_1')
        self.assertEqual(gene.strand, core.PolymerStrand.negative)
        self.assertEqual(gene.start, 1)
        self.assertEqual(gene.end, 2)
        self.assertEqual(gene.comments, '')
        self.assertEqual(gene.references, [])
        self.assertEqual(gene.database_references, [])


class ExonLocusTestCase(unittest.TestCase):
    def test_constructor(self):
        exon = eukaryote_schema.ExonLocus(id='e1', name='exon1', start=21, end=30, 
            comments='No comment')
        self.assertEqual(exon.id, 'e1')
        self.assertEqual(exon.name, 'exon1')
        self.assertEqual(exon.start, 21)
        self.assertEqual(exon.end, 30)
        self.assertEqual(exon.comments, 'No comment')
        self.assertEqual(exon.references, [])
        self.assertEqual(exon.database_references, [])


class CdsLocusTestCase(unittest.TestCase):
    def test_constructor(self):
        cds1 = eukaryote_schema.CdsLocus(id='c1', name='cds1', start=21, end=30, 
            comments='No comment')
        self.assertEqual(cds1.id, 'c1')
        self.assertEqual(cds1.name, 'cds1')
        self.assertEqual(cds1.start, 21)
        self.assertEqual(cds1.end, 30)
        self.assertEqual(cds1.comments, 'No comment')
        self.assertEqual(cds1.references, [])
        self.assertEqual(cds1.database_references, [])


class RegulatoryElementLocusTestCase(unittest.TestCase):
    def test_constructor(self):        
        
        promoter = eukaryote_schema.RegulatoryElementLocus(id='p1', name='promoter1', 
            type=eukaryote_schema.RegulatoryElementType.promoter, 
            activity=eukaryote_schema.ActivityLevel.active, 
            strand=core.PolymerStrand.positive, start=2, end=10)

        self.assertEqual(promoter.id, 'p1')
        self.assertEqual(promoter.name, 'promoter1')
        self.assertEqual(promoter.type, eukaryote_schema.RegulatoryElementType.promoter)
        self.assertEqual(promoter.activity, eukaryote_schema.ActivityLevel.active)
        self.assertEqual(promoter.strand, core.PolymerStrand.positive)
        self.assertEqual(promoter.start, 2)
        self.assertEqual(promoter.end, 10)

        tf1 = eukaryote_schema.ProteinSpeciesType(id='TF1')
        
        TF_binding_site = eukaryote_schema.RegulatoryElementLocus( 
            type=eukaryote_schema.RegulatoryElementType.TF_binding_site, 
            start=2, end=10, bound_start=3, bound_end=8, motif_features=[tf1])
        
        self.assertEqual(TF_binding_site.type.value, 5)
        self.assertEqual(TF_binding_site.start, 2)
        self.assertEqual(TF_binding_site.end, 10)
        self.assertEqual(TF_binding_site.bound_start, 3)
        self.assertEqual(TF_binding_site.bound_end, 8)
        self.assertEqual(TF_binding_site.motif_features[0].id, 'TF1')


class RegulatoryModuleTestCase(unittest.TestCase):
    def test_constructor(self):
        dna1 = core.DnaSpeciesType(circular=False, double_stranded=False)
        
        gene1 = eukaryote_schema.GeneLocus(polymer=dna1, start=9, end=15)
        gene2 = eukaryote_schema.GeneLocus(polymer=dna1, start=17, end=18)

        promoter = eukaryote_schema.RegulatoryElementLocus(polymer=dna1, start=6, end=8)
        enhancer = eukaryote_schema.RegulatoryElementLocus(polymer=dna1, start=2, end=4)

        tf = eukaryote_schema.ProteinSpeciesType(id='tf')

        reg_module1 = eukaryote_schema.RegulatoryModule(gene=gene1, 
            regulatory_element=promoter, binding_factor=tf, 
            type=eukaryote_schema.RegulationType.proximal, 
            direction=eukaryote_schema.RegulatoryDirection.positive)
        reg_module2 = eukaryote_schema.RegulatoryModule(gene=gene1, 
            regulatory_element=enhancer, binding_factor=tf, 
            type=eukaryote_schema.RegulationType.distal, 
            direction=eukaryote_schema.RegulatoryDirection.negative)
        reg_module3 = eukaryote_schema.RegulatoryModule(id='rm3', name='reg_module3', 
            gene=gene2, regulatory_element=enhancer)
         
        self.assertEqual(reg_module1.gene, gene1)
        self.assertEqual(reg_module1.regulatory_element, promoter)
        self.assertEqual(reg_module1.binding_factor, tf)
        self.assertEqual(reg_module1.type.value, 1)
        self.assertEqual(reg_module1.direction.value, 1)
        self.assertEqual(reg_module2.gene, gene1)
        self.assertEqual(reg_module2.regulatory_element, enhancer)
        self.assertEqual(reg_module2.binding_factor, tf)
        self.assertEqual(reg_module2.type.value, 2)
        self.assertEqual(reg_module2.direction.value, -1)
        self.assertEqual(reg_module3.id, 'rm3')
        self.assertEqual(reg_module3.name, 'reg_module3')
        self.assertEqual(set([i.gene for i in enhancer.regulatory_modules]), set([gene1, gene2]))
