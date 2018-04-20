""" Tests of the knowledge base IO

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-02-07
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_kb import core
from wc_kb import io
import Bio.Seq
import obj_model.io
import os
import random
import shutil
import tempfile
import unittest
import wc_utils.workbook.io


class TestIO(unittest.TestCase):

    def setUp(self):
        self.kb = kb = core.KnowledgeBase(id='genus_species', name='Genus species', version='0.0.1')

        cell = kb.cell = core.Cell(id='genus_species_cell')

        for i_chr in range(5):
            dna = core.DnaSpeciesType(id='chr_{}'.format(i_chr + 1))
            cell.species_types.append(dna)

            seq_len = random.randint(100, 200)
            bases = 'ACGT'
            seq = ''
            for i_nt in range(seq_len):
                seq += bases[random.randint(0, 3)]
            dna.seq = Bio.Seq.Seq(seq)

            for i_trn in range(5):
                trn = core.TranscriptionUnitLocus(id='tu_{}_{}'.format(i_chr + 1, i_trn + 1))
                dna.loci.append(trn)
                trn.start = random.randint(100, 200)
                trn.end = ((trn.start + random.randint(1, 200) - 1) % seq_len) + 1
                trn.strand = core.PolymerStrand.positive

        self.dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dir)

    def test_write_read(self):
        core_path = os.path.join(self.dir, 'core.xlsx')
        seq_path = os.path.join(self.dir, 'seq.fna')

        writer = io.Writer()
        writer.run(self.kb, core_path, seq_path)

        reader = io.Reader()
        kb = reader.run(core_path, seq_path)

        self.assertTrue(kb.is_equal(self.kb))

    def test_write_read_sloppy(self):
        core_path = os.path.join(self.dir, 'core.xlsx')
        seq_path = os.path.join(self.dir, 'seq.fna')

        writer = io.Writer()
        writer.run(self.kb, core_path, seq_path)

        wb = wc_utils.workbook.io.read(core_path)
        row = wb['Knowledge base'].pop(0)
        wb['Knowledge base'].insert(1, row)
        wc_utils.workbook.io.write(core_path, wb)

        reader = io.Reader()
        with self.assertRaisesRegexp(ValueError, 'The attributes must be defined in this order'):
            kb = reader.run(core_path, seq_path)
        kb = reader.run(core_path, seq_path, strict=False)

        self.assertTrue(kb.is_equal(self.kb))

    def test_reader_no_kb(self):
        core_path = os.path.join(self.dir, 'core.xlsx')
        obj_model.io.WorkbookWriter().run(core_path, [], io.Writer.model_order)

        seq_path = os.path.join(self.dir, 'seq.fna')
        with open(seq_path, 'w') as file:
            pass

        kb = io.Reader().run(core_path, seq_path)
        self.assertEqual(kb, None)

    def test_reader_error_multiple_kbs(self):
        kb1 = core.KnowledgeBase(id='kb1', name='kb1', version='0.0.1')
        kb2 = core.KnowledgeBase(id='kb2', name='kb2', version='0.0.1')

        core_path = os.path.join(self.dir, 'core.xlsx')
        obj_model.io.WorkbookWriter().run(core_path, [kb1, kb2], io.Writer.model_order)

        seq_path = os.path.join(self.dir, 'seq.fna')
        with open(seq_path, 'w') as file:
            pass

        with self.assertRaisesRegexp(ValueError, ' should only define one knowledge base'):
            io.Reader().run(core_path, seq_path)

    def test_convert(self):
        path_core_1 = os.path.join(self.dir, 'core_1.xlsx')
        path_core_2 = os.path.join(self.dir, 'core_2-*.csv')
        path_core_3 = os.path.join(self.dir, 'core_3.xlsx')
        path_seq_1 = os.path.join(self.dir, 'seq_1.fna')
        path_seq_2 = os.path.join(self.dir, 'seq_2.fna')
        path_seq_3 = os.path.join(self.dir, 'seq_3.fna')

        io.Writer().run(self.kb, path_core_1, path_seq_1)

        io.convert(path_core_1, path_seq_1, path_core_2, path_seq_2)
        kb = io.Reader().run(path_core_2, path_seq_2)
        self.assertTrue(kb.is_equal(self.kb))

        io.convert(path_core_2, path_seq_2, path_core_3, path_seq_3)
        kb = io.Reader().run(path_core_3, path_seq_3)
        self.assertTrue(kb.is_equal(self.kb))

    def test_convert_sloppy(self):
        path_core_1 = os.path.join(self.dir, 'core_1.xlsx')
        path_core_2 = os.path.join(self.dir, 'core_2-*.csv')
        path_core_3 = os.path.join(self.dir, 'core_3.xlsx')
        path_seq_1 = os.path.join(self.dir, 'seq_1.fna')
        path_seq_2 = os.path.join(self.dir, 'seq_2.fna')
        path_seq_3 = os.path.join(self.dir, 'seq_3.fna')

        io.Writer().run(self.kb, path_core_1, path_seq_1)

        wb = wc_utils.workbook.io.read(path_core_1)
        row = wb['Knowledge base'].pop(0)
        wb['Knowledge base'].insert(1, row)
        wc_utils.workbook.io.write(path_core_1, wb)

        with self.assertRaisesRegexp(ValueError, 'The attributes must be defined in this order'):
            io.convert(path_core_1, path_seq_1, path_core_2, path_seq_2)
        io.convert(path_core_1, path_seq_1, path_core_2, path_seq_2, strict=False)
        kb = io.Reader().run(path_core_2, path_seq_2)
        self.assertTrue(kb.is_equal(self.kb))

        io.convert(path_core_2, path_seq_2, path_core_3, path_seq_3)
        kb = io.Reader().run(path_core_3, path_seq_3)
        self.assertTrue(kb.is_equal(self.kb))

    def test_create_template(self):
        path_core = os.path.join(self.dir, 'core.xlsx')
        path_seq = os.path.join(self.dir, 'seq.fna')
        io.create_template(path_core, path_seq)
        kb = io.Reader().run(path_core, path_seq)
