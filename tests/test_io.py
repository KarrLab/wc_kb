""" Tests of the knowledge base IO

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2018-02-07
:Copyright: 2018, Karr Lab
:License: MIT
"""

from test.support import EnvironmentVarGuard
from wc_kb import core, prokaryote_schema
from wc_kb import io
import Bio.Seq
import Bio.SeqRecord
import filecmp
import obj_model.io
import os
import random
import shutil
import tempfile
import unittest
import wc_utils.workbook.io

class TestIO(unittest.TestCase):

    def setUp(self):

        self.dir = tempfile.mkdtemp()
        self.seq_path = os.path.join(self.dir, 'seq.fna')

        self.kb = kb = core.KnowledgeBase(id='genus_species', name='Genus species', version='0.0.1')

        cell = kb.cell = core.Cell(id='genus_species_cell')

        dna_seqs = []
        for i_chr in range(5):
            dna = core.DnaSpeciesType(id='chr_{}'.format(i_chr + 1), sequence_path=self.seq_path)
            cell.species_types.append(dna)

            seq_len = random.randint(100, 200)
            bases = 'ACGT'
            seq = ''
            for i_nt in range(seq_len):
                seq += bases[random.randint(0, 3)]
            dna_seqs.append(Bio.SeqRecord.SeqRecord(
                Bio.Seq.Seq(seq), dna.id))

            for i_trn in range(5):
                trn = prokaryote_schema.TranscriptionUnitLocus(id='tu_{}_{}'.format(i_chr + 1, i_trn + 1))
                trn.type = random.choice(['mRna', 'sRna', 'tRna', 'rRna', 'intergenic', 'mixed'])
                trn.cell = cell
                dna.loci.append(trn)
                trn.start = random.randint(100, 200)
                trn.end = ((trn.start + random.randint(1, 200) - 1) % seq_len) + 1
                trn.strand = core.PolymerStrand.positive

        with open(self.seq_path, 'w') as file:
            writer = Bio.SeqIO.FastaIO.FastaWriter(
                file, wrap=70, record2title=lambda record: record.id)
            writer.write_file(dna_seqs)

    def tearDown(self):
        shutil.rmtree(self.dir)

    def test_write_read(self):
        core_path = os.path.join(self.dir, 'core.xlsx')

        writer = io.Writer()
        writer.run(core_path, self.kb, set_repo_metadata_from_path=False)

        reader = io.Reader()
        kb = reader.run(core_path, seq_path=self.seq_path)[core.KnowledgeBase][0]

        core_path = os.path.join(self.dir, 'core2.xlsx')
        seq_path = os.path.join(self.dir, 'seq2.fna')
        writer.run(core_path, kb, seq_path, set_repo_metadata_from_path=False)

        self.assertTrue(self.kb.is_equal(kb))
        self.assertTrue(filecmp.cmp(self.seq_path, seq_path, shallow=False))

    def test_read_write_prokaryote(self):
        fixtures = os.path.join(os.path.dirname(__file__), 'fixtures')
        core_path = os.path.join(fixtures, 'prokaryote_core.xlsx')
        seq_path = os.path.join(fixtures, 'seq.fna')

        reader = io.Reader()
        kb = reader.run(core_path, seq_path=seq_path)[core.KnowledgeBase][0]

        tmp_core_path = os.path.join(self.dir, 'tmp_core.xlsx')
        tmp_seq_path = os.path.join(self.dir, 'tmp_seq.fna')

        writer = io.Writer()
        writer.run(tmp_core_path, kb, seq_path=tmp_seq_path, set_repo_metadata_from_path=False)

        tmp_kb = reader.run(tmp_core_path, seq_path)[core.KnowledgeBase][0]

        self.assertTrue(kb.is_equal(tmp_kb))
        self.assertTrue(filecmp.cmp(tmp_seq_path, seq_path, shallow=False))

    def test_read_write_eukaryote(self):
        fixtures = os.path.join(os.path.dirname(__file__), 'fixtures')
        core_path = os.path.join(fixtures, 'eukaryote_core.xlsx')
        seq_path = os.path.join(fixtures, 'eukaryote_seq.fna')

        reader = io.Reader()
        kb = reader.run(core_path, seq_path=seq_path, taxon='eukaryote')[core.KnowledgeBase][0]

        tmp_core_path = os.path.join(self.dir, 'tmp_eukaryote_core.xlsx')
        tmp_seq_path = os.path.join(self.dir, 'tmp_eukaryote_seq.fna')

        writer = io.Writer()
        writer.run(tmp_core_path, kb, seq_path=tmp_seq_path, taxon='eukaryote', set_repo_metadata_from_path=False)

        tmp_kb = reader.run(tmp_core_path, seq_path, taxon='eukaryote')[core.KnowledgeBase][0]

        self.assertTrue(kb.is_equal(tmp_kb))
        self.assertTrue(filecmp.cmp(tmp_seq_path, seq_path, shallow=False))

    def test_rewrite_seq_path_in_read_write(self):
        path_core_1 = os.path.join(self.dir, 'core_1.xlsx')
        path_core_2 = os.path.join(self.dir, 'core_2.xlsx')
        path_seq_1 = os.path.join(self.dir, 'seq_1.fna')
        path_seq_2 = os.path.join(self.dir, 'seq_2.fna')

        io.Writer().run(path_core_1, self.kb, seq_path=path_seq_1, set_repo_metadata_from_path=False)
        kb1 = io.Reader().run(path_core_1, seq_path=path_seq_1)[core.KnowledgeBase][0]
        kb2 = io.Reader().run(path_core_1, seq_path=path_seq_1, rewrite_seq_path=False)[core.KnowledgeBase][0]
        kb3 = io.Reader().run(path_core_1, seq_path=self.seq_path)[core.KnowledgeBase][0]
        kb4 = io.Reader().run(path_core_1, seq_path=self.seq_path, rewrite_seq_path=False)[core.KnowledgeBase][0]
        self.assertFalse(kb1.is_equal(self.kb))
        self.assertFalse(kb2.is_equal(self.kb))
        self.assertTrue(kb3.is_equal(self.kb))
        self.assertFalse(kb4.is_equal(self.kb))
        self.assertTrue(filecmp.cmp(path_seq_1, self.seq_path, shallow=False))

        io.Writer().run(path_core_2, self.kb, seq_path=path_seq_2, rewrite_seq_path=False, set_repo_metadata_from_path=False)
        kb5 = io.Reader().run(path_core_2, seq_path=path_seq_2)[core.KnowledgeBase][0]
        kb6 = io.Reader().run(path_core_2, seq_path=path_seq_2, rewrite_seq_path=False)[core.KnowledgeBase][0]
        kb7 = io.Reader().run(path_core_2, seq_path=self.seq_path)[core.KnowledgeBase][0]
        kb8 = io.Reader().run(path_core_2, seq_path=self.seq_path, rewrite_seq_path=False)[core.KnowledgeBase][0]
        self.assertFalse(kb5.is_equal(self.kb))
        self.assertTrue(kb6.is_equal(self.kb))
        self.assertTrue(kb7.is_equal(self.kb))
        self.assertTrue(kb8.is_equal(self.kb))
        self.assertTrue(filecmp.cmp(path_seq_2, self.seq_path, shallow=False))

    def test_write_with_repo_md(self):

        _, core_path = tempfile.mkstemp(suffix='.xlsx', dir='.')
        _, seq_path = tempfile.mkstemp(suffix='.fna', dir='.')

        self.assertEqual(self.kb.url, '')

        writer = io.Writer()
        writer.run(core_path, self.kb, seq_path=seq_path, set_repo_metadata_from_path=True)

        self.assertIn(self.kb.url, [
            'https://github.com/KarrLab/wc_kb.git',
            'ssh://git@github.com/KarrLab/wc_kb.git',
            'git@github.com:KarrLab/wc_kb.git',
        ])

        os.remove(core_path)
        os.remove(seq_path)

    def test_write_without_cell_relationships(self):
        core_path = os.path.join(self.dir, 'core.xlsx')
        seq_path = os.path.join(self.dir, 'test_seq.fna')

        with open(seq_path, 'w') as file:
            file.write('>chr_x\nACGT\n')

        dna = core.DnaSpeciesType(id='chr_x', sequence_path=seq_path)
        self.kb.cell.species_types.append(dna)

        trn = prokaryote_schema.TranscriptionUnitLocus(id='tu_x_0')
        dna.loci.append(trn)
        trn.cell = None

        writer = io.Writer()
        with self.assertRaisesRegex(ValueError, 'must be set to the instance of `Cell`'):
            writer.run(core_path, self.kb, seq_path=seq_path, set_repo_metadata_from_path=False)

    def test_write_read_sloppy(self):
        core_path = os.path.join(self.dir, 'core.xlsx')
        seq_path = os.path.join(self.dir, 'test_seq.fna')

        writer = io.Writer()
        writer.run(core_path, self.kb, seq_path=seq_path, set_repo_metadata_from_path=False)

        wb = wc_utils.workbook.io.read(core_path)

        row = wb['KB'].pop(0)
        wb['KB'].insert(1, row)
        wc_utils.workbook.io.write(core_path, wb)

        reader = io.Reader()
        with self.assertRaisesRegex(ValueError, "The model cannot be loaded because"):
            reader.run(core_path, seq_path=self.seq_path)
        env = EnvironmentVarGuard()
        env.set('CONFIG__DOT__wc_kb__DOT__io__DOT__strict', '0')
        with env:
            kb = reader.run(core_path, self.seq_path)[core.KnowledgeBase][0]

        self.assertTrue(kb.is_equal(self.kb))
        self.assertTrue(filecmp.cmp(self.seq_path, seq_path, shallow=False))

    def test_reader_no_kb(self):
        core_path = os.path.join(self.dir, 'core.xlsx')
        obj_model.io.WorkbookWriter().run(core_path, [], io.PROKARYOTE_MODELS, include_all_attributes=False)

        seq_path = os.path.join(self.dir, 'test_seq.fna')
        with open(seq_path, 'w') as file:
            pass

        with self.assertRaisesRegex(ValueError, 'should define one knowledge base'):
            io.Reader().run(core_path, seq_path=seq_path)

        obj_model.io.WorkbookWriter().run(core_path, [core.Cell(id='cell')], io.PROKARYOTE_MODELS, include_all_attributes=False)
        with self.assertRaisesRegex(ValueError, 'should define one knowledge base'):
            io.Reader().run(core_path, seq_path=seq_path)

    def test_reader_error_multiple_kbs(self):
        kb1 = core.KnowledgeBase(id='kb1', name='kb1', version='0.0.1')
        kb2 = core.KnowledgeBase(id='kb2', name='kb2', version='0.0.1')

        core_path = os.path.join(self.dir, 'core.xlsx')
        obj_model.io.WorkbookWriter().run(core_path, [kb1, kb2], io.PROKARYOTE_MODELS, include_all_attributes=False)

        seq_path = os.path.join(self.dir, 'test_seq.fna')
        with open(seq_path, 'w') as file:
            pass

        with self.assertRaisesRegex(ValueError, ' should define one knowledge base'):
            io.Reader().run(core_path, seq_path=seq_path)

    def test_reader_no_cell(self):
        kb = core.KnowledgeBase(id='kb', name='kb1', version='0.0.1')
        dna = core.DnaSpeciesType(id='chr')

        core_path = os.path.join(self.dir, 'core.xlsx')
        obj_model.io.WorkbookWriter().run(core_path, [kb, dna], io.PROKARYOTE_MODELS, include_all_attributes=False)

        seq_path = os.path.join(self.dir, 'test_seq.fna')
        with open(seq_path, 'w') as file:
            pass

        io.Reader().run(core_path, seq_path=seq_path)

    def test_reader_error_multiple_cells(self):
        kb = core.KnowledgeBase(id='kb', name='kb1', version='0.0.1')
        cell1 = core.Cell(id='cell1', name='cell1')
        cell2 = core.Cell(id='cell2', name='cell2')

        core_path = os.path.join(self.dir, 'core.xlsx')
        obj_model.io.WorkbookWriter().run(core_path, [kb, cell1, cell2], io.PROKARYOTE_MODELS, include_all_attributes=False)

        seq_path = os.path.join(self.dir, 'test_seq.fna')
        with open(seq_path, 'w') as file:
            pass

        with self.assertRaisesRegex(ValueError, ' should define zero or one cells'):
            io.Reader().run(core_path, seq_path=seq_path)

    def test_convert(self):
        path_core_1 = os.path.join(self.dir, 'core_1.xlsx')
        path_core_2 = os.path.join(self.dir, 'core_2-*.csv')
        path_core_3 = os.path.join(self.dir, 'core_3.xlsx')
        path_seq_1 = os.path.join(self.dir, 'seq_1.fna')
        path_seq_2 = os.path.join(self.dir, 'seq_2.fna')
        path_seq_3 = os.path.join(self.dir, 'seq_3.fna')

        io.Writer().run(path_core_1, self.kb, seq_path=path_seq_1, set_repo_metadata_from_path=False)
        self.assertTrue(filecmp.cmp(path_seq_1, self.seq_path, shallow=False))

        io.convert(path_core_1, path_seq_1, path_core_2, path_seq_2)
        kb = io.Reader().run(path_core_2, seq_path=self.seq_path)[core.KnowledgeBase][0]
        self.assertTrue(kb.is_equal(self.kb))
        self.assertTrue(filecmp.cmp(path_seq_1, path_seq_2, shallow=False))

        io.convert(path_core_2, path_seq_2, path_core_3, path_seq_3)
        kb = io.Reader().run(path_core_3, seq_path=self.seq_path)[core.KnowledgeBase][0]
        self.assertTrue(kb.is_equal(self.kb))
        self.assertTrue(filecmp.cmp(path_seq_2, path_seq_3, shallow=False))

    def test_convert_sloppy(self):
        path_core_1 = os.path.join(self.dir, 'core_1.xlsx')
        path_core_2 = os.path.join(self.dir, 'core_2-*.csv')
        path_core_3 = os.path.join(self.dir, 'core_3.xlsx')
        path_seq_1 = os.path.join(self.dir, 'seq_1.fna')
        path_seq_2 = os.path.join(self.dir, 'seq_2.fna')
        path_seq_3 = os.path.join(self.dir, 'seq_3.fna')

        io.Writer().run(path_core_1, self.kb, seq_path=path_seq_1, set_repo_metadata_from_path=False)
        self.assertTrue(filecmp.cmp(path_seq_1, self.seq_path, shallow=False))

        wb = wc_utils.workbook.io.read(path_core_1)
        row = wb['KB'].pop(0)
        wb['KB'].insert(1, row)
        wc_utils.workbook.io.write(path_core_1, wb)

        with self.assertRaisesRegex(ValueError, "The model cannot be loaded because"):
            io.convert(path_core_1, path_seq_1, path_core_2, path_seq_2)
        env = EnvironmentVarGuard()
        env.set('CONFIG__DOT__wc_kb__DOT__io__DOT__strict', '0')
        with env:
            io.convert(path_core_1, path_seq_1, path_core_2, path_seq_2)
        kb = io.Reader().run(path_core_2, seq_path=self.seq_path)[core.KnowledgeBase][0]
        self.assertTrue(kb.is_equal(self.kb))
        self.assertTrue(filecmp.cmp(path_seq_1, path_seq_2, shallow=False))

        io.convert(path_core_2, path_seq_2, path_core_3, path_seq_3)
        kb = io.Reader().run(path_core_3, seq_path=self.seq_path)[core.KnowledgeBase][0]
        self.assertTrue(kb.is_equal(self.kb))
        self.assertTrue(filecmp.cmp(path_seq_2, path_seq_3, shallow=False))

    def test_create_template(self):
        path_core = os.path.join(self.dir, 'template.xlsx')
        path_seq = os.path.join(self.dir, 'template_seq.fna')
        io.create_template(path_core, path_seq, set_repo_metadata_from_path=False)
        kb = io.Reader().run(path_core, seq_path=path_seq)[core.KnowledgeBase][0]

    def test_validate_implicit_relationships(self):
        class TestModel(obj_model.Model):
            id = obj_model.StringAttribute(primary=True, unique=True)

        try:
            core.KnowledgeBase.Meta.attributes['test'] = obj_model.OneToOneAttribute(TestModel, related_name='a')
            with self.assertRaisesRegex(Exception, 'Relationships from `KnowledgeBase` not supported:'):
                io.Writer.validate_implicit_relationships()
        finally:
            core.KnowledgeBase.Meta.attributes.pop('test')

        try:
            core.KnowledgeBase.Meta.related_attributes['test'] = obj_model.OneToManyAttribute(core.Cell, related_name='c')
            with self.assertRaisesRegex(Exception,
                                        'Relationships to `KnowledgeBase` that are not one-to-one are prohibited'):
                io.Writer.validate_implicit_relationships()
        finally:
            core.KnowledgeBase.Meta.related_attributes.pop('test')

        try:
            core.Cell.Meta.attributes['test'] = obj_model.OneToManyAttribute(TestModel, related_name='c')
            with self.assertRaisesRegex(Exception,
                                        'Relationships from `Cell` to `KnowledgeBase` that are not one-to-one are prohibited:'):
                io.Writer.validate_implicit_relationships()
        finally:
            core.Cell.Meta.attributes.pop('test')

        try:
            core.Cell.Meta.attributes['test'] = obj_model.OneToOneAttribute(TestModel, related_name='d')
            with self.assertRaisesRegex(Exception,
                                        'Relationships from `Cell` to classes other than `KnowledgeBase` are prohibited:'):
                io.Writer.validate_implicit_relationships()
        finally:
            core.Cell.Meta.attributes.pop('test')

        try:
            core.Cell.Meta.related_attributes['test'] = obj_model.OneToManyAttribute(TestModel, related_name='d')
            with self.assertRaisesRegex(Exception,
                                        'Relationships to `Cell` that are not one-to-one or many-to-one are prohibited: '):
                io.Writer.validate_implicit_relationships()
        finally:
            core.Cell.Meta.related_attributes.pop('test')

        try:
            core.KnowledgeBase.Meta.related_attributes['test'] = obj_model.OneToOneAttribute(TestModel, related_name='b')
            with self.assertRaisesRegex(Exception,
                                        'Relationships to `KnowledgeBase` from classes other than `Cell` are prohibited'):
                io.Writer.validate_implicit_relationships()
        finally:
            core.KnowledgeBase.Meta.related_attributes.pop('test')
