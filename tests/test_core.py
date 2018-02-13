""" Tests of the knowledge base

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-02-07
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_kb import core
from wc_utils.util import chem
import Bio.Seq
import Bio.SeqIO
import mendeleev
import unittest


class TestCore(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        records = Bio.SeqIO.parse('tests/fixtures/seq.fna', 'fasta')
        cls.chr_seq = next(records).seq

    def test_ChromosomeStrand(self):
        self.assertEqual(core.ChromosomeStrand.negative.value, -1)
        self.assertEqual(core.ChromosomeStrand.positive.value, 1)

    def test_KnowledgeBase(self):
        kb = core.KnowledgeBase()
        self.assertEqual(kb.cell, None)

    def test_Cell(self):
        cell = core.Cell()
        self.assertEqual(cell.knowledge_base, None)
        self.assertEqual(cell.chromosomes, [])

        chromosomes = [
            core.Chromosome(id='chr1', seq='AAA'),
            core.Chromosome(id='chr2', seq='CCC'),
        ]
        cell = core.Cell(chromosomes=chromosomes)
        self.assertEqual(cell.chromosomes, chromosomes)

    def test_Chromosome(self):
        chr = core.Chromosome(seq=Bio.Seq.Seq('AAA'))

    def test_Chromosome_get_len(self):
        chr = core.Chromosome(seq=Bio.Seq.Seq('AAA'))
        self.assertEqual(chr.get_len(), 3)

    def test_Chromosome_get_subseq(self):
        chr = core.Chromosome(seq=Bio.Seq.Seq('AAATGCCC'))
        self.assertEqual(chr.get_subseq(3, 6), 'ATGC')
        self.assertEqual(chr.get_subseq(3, 6, strand=core.ChromosomeStrand.positive), 'ATGC')
        self.assertEqual(chr.get_subseq(3, 6, strand=core.ChromosomeStrand.negative), 'GCAT')

        self.assertEqual(chr.get_subseq(0, 1), 'CA')
        self.assertEqual(chr.get_subseq(-1, 1), 'CCA')
        self.assertEqual(chr.get_subseq(-3, 1), 'GCCCA')

        self.assertEqual(chr.get_subseq(6, 10), 'CCCAA')
        self.assertEqual(chr.get_subseq(6, 10, strand=core.ChromosomeStrand.negative), 'TTGGG')

        self.assertEqual(chr.get_subseq(6, 26), 'CCCAAATGCCCAAATGCCCAA')
        self.assertEqual(chr.get_subseq(-2, 18), 'CCCAAATGCCCAAATGCCCAA')
        self.assertEqual(chr.get_subseq(-10, 10), 'CCCAAATGCCCAAATGCCCAA')

    def test_Chromosome_get_subseq_real_genes(self):
        chr = core.Chromosome(seq=self.chr_seq)

        # MPN001
        gene_seq = chr.get_subseq(692, 1834)
        self.assertEqual(gene_seq[0:10], 'ATGAAAGTTT')
        self.assertEqual(gene_seq[-10:], 'TTCCAAGTAA')
        self.assertEqual(gene_seq.transcribe()[0:-3].translate(table=4)[0:10], 'MKVLINKNEL')

        # MPN011
        gene_seq = chr.get_subseq(12838, 13533, strand=core.ChromosomeStrand.negative)
        self.assertEqual(gene_seq[0:10], 'ATGAAATTTA')
        self.assertEqual(gene_seq[-10:], 'AATTGAGTAA')
        self.assertEqual(gene_seq.transcribe().translate(table=4)[0:10], 'MKFKFLLTPL')

    def test_TranscriptionUnit(self):
        tu = core.TranscriptionUnit(id='tu_1',
                                          start=1, end=2, strand=core.ChromosomeStrand.positive,
                                          pribnow_start=3, pribnow_end=4)

        self.assertEqual(tu.id, 'tu_1')
        self.assertEqual(tu.start, 1)
        self.assertEqual(tu.end, 2)
        self.assertEqual(tu.strand, core.ChromosomeStrand.positive)
        self.assertEqual(tu.pribnow_start, 3)
        self.assertEqual(tu.pribnow_end, 4)

    def test_TranscriptionUnit_get_3_prime(self):
        tu = core.TranscriptionUnit(start=100, end=200, strand=core.ChromosomeStrand.positive)
        self.assertEqual(tu.get_3_prime(), 200)

        tu = core.TranscriptionUnit(start=100, end=200, strand=core.ChromosomeStrand.negative)
        self.assertEqual(tu.get_3_prime(), 100)

    def test_TranscriptionUnit_get_5_prime(self):
        tu = core.TranscriptionUnit(start=100, end=200, strand=core.ChromosomeStrand.positive)
        self.assertEqual(tu.get_5_prime(), 100)

        tu = core.TranscriptionUnit(start=100, end=200, strand=core.ChromosomeStrand.negative)
        self.assertEqual(tu.get_5_prime(), 200)

    def test_TranscriptionUnit_get_len(self):
        tu = core.TranscriptionUnit(start=100, end=200)
        self.assertEqual(tu.get_len(), 101)

    def test_TranscriptionUnit_get_seq(self):
        chr = core.Chromosome(seq=Bio.Seq.Seq('AAATGCCC'))

        tu = chr.transcription_units.create(start=3, end=6)
        self.assertEqual(tu.get_seq(), 'AUGC')

        tu = chr.transcription_units.create(start=3, end=6, strand=core.ChromosomeStrand.positive)
        self.assertEqual(tu.get_seq(), 'AUGC')

        tu = chr.transcription_units.create(start=3, end=6, strand=core.ChromosomeStrand.negative)
        self.assertEqual(tu.get_seq(), 'GCAU')

    def test_TranscriptionUnit_get_empirical_formula(self):
        chr = core.Chromosome()

        tu = core.TranscriptionUnit(chromosome=chr, start=1, end=1)
        chr.seq = Bio.Seq.Seq('A')
        self.assertEqual(tu.get_empirical_formula(), chem.EmpiricalFormula('C10H12N5O7P'))

        tu = core.TranscriptionUnit(chromosome=chr, start=1, end=1)
        chr.seq = Bio.Seq.Seq('C')
        self.assertEqual(tu.get_empirical_formula(), chem.EmpiricalFormula('C9H12N3O8P'))

        tu = core.TranscriptionUnit(chromosome=chr, start=1, end=1)
        chr.seq = Bio.Seq.Seq('G')
        self.assertEqual(tu.get_empirical_formula(), chem.EmpiricalFormula('C10H12N5O8P'))

        tu = core.TranscriptionUnit(chromosome=chr, start=1, end=1)
        chr.seq = Bio.Seq.Seq('T')
        self.assertEqual(tu.get_empirical_formula(), chem.EmpiricalFormula('C9H11N2O9P'))

        tu = core.TranscriptionUnit(chromosome=chr, start=1, end=2)
        chr.seq = Bio.Seq.Seq('AA')
        self.assertEqual(tu.get_empirical_formula(), chem.EmpiricalFormula('C20H23N10O13P2'))

    def test_TranscriptionUnit_get_charge(self):
        tu = core.TranscriptionUnit(start=1, end=1)
        self.assertEqual(tu.get_charge(), -2)

        tu = core.TranscriptionUnit(start=1, end=2)
        self.assertEqual(tu.get_charge(), -3)

    def test_TranscriptionUnit_get_molecular_weight(self):
        chr = core.Chromosome(seq=Bio.Seq.Seq('AACCGGTT'))

        tu = core.TranscriptionUnit(chromosome=chr, start=1, end=1)
        self.assertAlmostEqual(
            tu.get_molecular_weight() + tu.get_charge() * mendeleev.element('H').atomic_weight,
            tu.get_empirical_formula().get_molecular_weight(), places=1)

        tu = core.TranscriptionUnit(chromosome=chr, start=3, end=3)
        self.assertAlmostEqual(
            tu.get_molecular_weight() + tu.get_charge() * mendeleev.element('H').atomic_weight,
            tu.get_empirical_formula().get_molecular_weight(), places=1)

        tu = core.TranscriptionUnit(chromosome=chr, start=5, end=5)
        self.assertAlmostEqual(
            tu.get_molecular_weight() + tu.get_charge() * mendeleev.element('H').atomic_weight,
            tu.get_empirical_formula().get_molecular_weight(), places=1)

        tu = core.TranscriptionUnit(chromosome=chr, start=7, end=7)
        self.assertAlmostEqual(
            tu.get_molecular_weight() + tu.get_charge() * mendeleev.element('H').atomic_weight,
            tu.get_empirical_formula().get_molecular_weight(), places=1)

        tu = core.TranscriptionUnit(chromosome=chr, start=1, end=2)
        self.assertAlmostEqual(
            tu.get_molecular_weight() + tu.get_charge() * mendeleev.element('H').atomic_weight,
            tu.get_empirical_formula().get_molecular_weight(), places=1)

    def test_TranscriptionUnit_get_pribnow_len(self):
        tu = core.TranscriptionUnit(pribnow_start=-40, pribnow_end=-50)
        self.assertEqual(tu.get_pribnow_len(), 11)

    def test_TranscriptionUnit_get_pribnow_seq(self):
        chr = core.Chromosome(seq=Bio.Seq.Seq('ACGTACGTACGTACGT'))
        tu = chr.transcription_units.create(
            start=7, end=8, strand=core.ChromosomeStrand.positive,
            pribnow_start=-2, pribnow_end=-4)
        self.assertEqual(tu.get_pribnow_seq(), 'ATG')

        chr = core.Chromosome(seq=Bio.Seq.Seq('ACGTACGTACGTACGT'))
        tu = chr.transcription_units.create(
            start=7, end=8, strand=core.ChromosomeStrand.negative,
            pribnow_start=-2, pribnow_end=-4)
        self.assertEqual(tu.get_pribnow_seq(), 'GCA')

        chr = core.Chromosome(seq=self.chr_seq)
        tu = chr.transcription_units.create(
            start=652, end=1900, strand=core.ChromosomeStrand.positive,
            pribnow_start=-35, pribnow_end=-40)
        self.assertEqual(tu.get_pribnow_seq(), 'TAAAAC')

        tu = chr.transcription_units.create(
            start=1252, end=1377, strand=core.ChromosomeStrand.negative,
            pribnow_start=-125, pribnow_end=-130)
        self.assertEqual(tu.get_pribnow_seq(), 'TAAGTT')
