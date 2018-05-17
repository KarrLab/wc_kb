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
                trn.cell = cell
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
        writer.run(self.kb, core_path, seq_path, set_repo_metadata_from_path=False)

        reader = io.Reader()
        kb = reader.run(core_path, seq_path)

        core_path = os.path.join(self.dir, 'core2.xlsx')
        seq_path = os.path.join(self.dir, 'seq2.fna')
        writer.run(kb, core_path, seq_path, set_repo_metadata_from_path=False)

        self.assertTrue(self.kb.is_equal(kb))

    def test_write_with_repo_md(self):
        _, core_path = tempfile.mkstemp(suffix='.xlsx', dir='.')
        _, seq_path = tempfile.mkstemp(suffix='.fna', dir='.')

        self.assertEqual(self.kb.url, '')

        writer = io.Writer()
        writer.run(self.kb, core_path, seq_path, set_repo_metadata_from_path=True)

        self.assertIn(self.kb.url, [
            'https://github.com/KarrLab/wc_kb.git',
            'ssh://git@github.com/KarrLab/wc_kb.git',
            'git@github.com:KarrLab/wc_kb.git',
        ])

        os.remove(core_path)
        os.remove(seq_path)

    def test_write_without_cell_relationships(self):
        core_path = os.path.join(self.dir, 'core.xlsx')
        seq_path = os.path.join(self.dir, 'seq.fna')

        dna = core.DnaSpeciesType(id='chr_x')
        dna.seq = Bio.Seq.Seq('ACGT')
        self.kb.cell.species_types.append(dna)

        trn = core.TranscriptionUnitLocus(id='tu_x_0')
        dna.loci.append(trn)
        trn.cell = None

        writer = io.Writer()
        with self.assertRaisesRegexp(ValueError, 'must be set to the instance of `Cell`'):
            writer.run(self.kb, core_path, seq_path, set_repo_metadata_from_path=False)

    def test_write_read_sloppy(self):
        core_path = os.path.join(self.dir, 'core.xlsx')
        seq_path = os.path.join(self.dir, 'seq.fna')

        writer = io.Writer()
        writer.run(self.kb, core_path, seq_path, set_repo_metadata_from_path=False)

        wb = wc_utils.workbook.io.read(core_path)
        row = wb['Knowledge base'].pop(0)
        wb['Knowledge base'].insert(1, row)
        wc_utils.workbook.io.write(core_path, wb)

        reader = io.Reader()
        with self.assertRaisesRegexp(ValueError, "The columns of worksheet 'Knowledge base' must be defined in this order"):
            kb = reader.run(core_path, seq_path)
        kb = reader.run(core_path, seq_path, strict=False)

        self.assertTrue(kb.is_equal(self.kb))

    def test_reader_no_kb(self):
        core_path = os.path.join(self.dir, 'core.xlsx')
        obj_model.io.WorkbookWriter().run(core_path, [], io.Writer.model_order, include_all_attributes=False)

        seq_path = os.path.join(self.dir, 'seq.fna')
        with open(seq_path, 'w') as file:
            pass

        kb = io.Reader().run(core_path, seq_path)
        self.assertEqual(kb, None)

        obj_model.io.WorkbookWriter().run(core_path, [core.Cell(id='cell')], io.Writer.model_order, include_all_attributes=False)
        with self.assertRaisesRegexp(ValueError, 'cannot contain instances'):
            io.Reader().run(core_path, seq_path)

    def test_reader_error_multiple_kbs(self):
        kb1 = core.KnowledgeBase(id='kb1', name='kb1', version='0.0.1')
        kb2 = core.KnowledgeBase(id='kb2', name='kb2', version='0.0.1')

        core_path = os.path.join(self.dir, 'core.xlsx')
        obj_model.io.WorkbookWriter().run(core_path, [kb1, kb2], io.Writer.model_order, include_all_attributes=False)

        seq_path = os.path.join(self.dir, 'seq.fna')
        with open(seq_path, 'w') as file:
            pass

        with self.assertRaisesRegexp(ValueError, ' should define one knowledge base'):
            io.Reader().run(core_path, seq_path)

    def test_reader_error_no_cell(self):
        kb = core.KnowledgeBase(id='kb', name='kb1', version='0.0.1')
        dna = core.DnaSpeciesType(id='chr')

        core_path = os.path.join(self.dir, 'core.xlsx')
        obj_model.io.WorkbookWriter().run(core_path, [kb, dna], io.Writer.model_order, include_all_attributes=False)

        seq_path = os.path.join(self.dir, 'seq.fna')
        with open(seq_path, 'w') as file:
            pass

        with self.assertRaisesRegexp(ValueError, 'cannot contain instances'):
            io.Reader().run(core_path, seq_path)

    def test_reader_error_multiple_cells(self):
        kb = core.KnowledgeBase(id='kb', name='kb1', version='0.0.1')
        cell1 = core.Cell(id='cell1', name='cell1')
        cell2 = core.Cell(id='cell2', name='cell2')

        core_path = os.path.join(self.dir, 'core.xlsx')
        obj_model.io.WorkbookWriter().run(core_path, [kb, cell1, cell2], io.Writer.model_order, include_all_attributes=False)

        seq_path = os.path.join(self.dir, 'seq.fna')
        with open(seq_path, 'w') as file:
            pass

        with self.assertRaisesRegexp(ValueError, ' should define one cell'):
            io.Reader().run(core_path, seq_path)

    def test_convert(self):
        path_core_1 = os.path.join(self.dir, 'core_1.xlsx')
        path_core_2 = os.path.join(self.dir, 'core_2-*.csv')
        path_core_3 = os.path.join(self.dir, 'core_3.xlsx')
        path_seq_1 = os.path.join(self.dir, 'seq_1.fna')
        path_seq_2 = os.path.join(self.dir, 'seq_2.fna')
        path_seq_3 = os.path.join(self.dir, 'seq_3.fna')

        io.Writer().run(self.kb, path_core_1, path_seq_1, set_repo_metadata_from_path=False)

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

        io.Writer().run(self.kb, path_core_1, path_seq_1, set_repo_metadata_from_path=False)

        wb = wc_utils.workbook.io.read(path_core_1)
        row = wb['Knowledge base'].pop(0)
        wb['Knowledge base'].insert(1, row)
        wc_utils.workbook.io.write(path_core_1, wb)

        with self.assertRaisesRegexp(ValueError, "The columns of worksheet 'Knowledge base' must be defined in this order"):
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
        io.create_template(path_core, path_seq, set_repo_metadata_from_path=False)
        kb = io.Reader().run(path_core, path_seq)

    def test_validate_implicit_relationships(self):
        class TestModel(obj_model.Model):
            id = obj_model.StringAttribute(primary=True, unique=True)

        core.KnowledgeBase.Meta.attributes['test'] = obj_model.OneToOneAttribute(TestModel, related_name='a')
        with self.assertRaisesRegexp(Exception, 'Relationships from `KnowledgeBase` not supported'):
            io.Writer.validate_implicit_relationships()
        core.KnowledgeBase.Meta.attributes.pop('test')

        core.KnowledgeBase.Meta.related_attributes['test'] = obj_model.OneToOneAttribute(TestModel, related_name='b')
        with self.assertRaisesRegexp(Exception, 'Only one-to-one relationships to `KnowledgeBase` from `Cell` are supported'):
            io.Writer.validate_implicit_relationships()
        core.KnowledgeBase.Meta.related_attributes.pop('test')

        core.Cell.Meta.attributes['test'] = obj_model.OneToManyAttribute(TestModel, related_name='c')
        with self.assertRaisesRegexp(Exception, 'Only one-to-one relationships from `Cell` to `KnowledgeBase` are supported'):
            io.Writer.validate_implicit_relationships()
        core.Cell.Meta.attributes.pop('test')

        core.Cell.Meta.related_attributes['test'] = obj_model.OneToManyAttribute(TestModel, related_name='d')
        with self.assertRaisesRegexp(Exception, 'Only one-to-one and many-to-one relationships are supported to `Cell`'):
            io.Writer.validate_implicit_relationships()
        core.Cell.Meta.related_attributes.pop('test')
