""" Tests of the knowledge base IO

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2018-02-07
:Copyright: 2018, Karr Lab
:License: MIT
"""

from obj_tables import utils
from test.support import EnvironmentVarGuard
from wc_kb import core, prokaryote
from wc_kb import io
from wc_utils.util.git import GitHubRepoForTests
import Bio.Seq
import Bio.SeqRecord
import filecmp
import obj_tables.io
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
                trn = prokaryote.TranscriptionUnitLocus(id='tu_{}_{}'.format(i_chr + 1, i_trn + 1))
                trn.type = random.choice(['WC:mRNA', 'WC:sRNA', 'WC:tRNA', 'WC:rRNA', 'WC:intergenic', 'WC:mixed'])
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
        writer.run(core_path, self.kb, data_repo_metadata=False)

        reader = io.Reader()
        kb = reader.run(core_path, seq_path=self.seq_path)[core.KnowledgeBase][0]

        core_path = os.path.join(self.dir, 'core2.xlsx')
        seq_path = os.path.join(self.dir, 'seq2.fna')
        writer.run(core_path, kb, seq_path, data_repo_metadata=False)

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
        writer.run(tmp_core_path, kb, seq_path=tmp_seq_path, data_repo_metadata=False)

        tmp_kb = reader.run(tmp_core_path, seq_path)[core.KnowledgeBase][0]

        self.assertTrue(kb.is_equal(tmp_kb))
        self.assertTrue(filecmp.cmp(tmp_seq_path, seq_path, shallow=False))

    def test_read_write_eukaryote(self):
        fixtures = os.path.join(os.path.dirname(__file__), 'fixtures')
        core_path = os.path.join(fixtures, 'eukaryote_core.xlsx')
        seq_path = os.path.join(fixtures, 'eukaryote_seq.fna')

        reader = io.Reader()
        kb = reader.run(core_path, seq_path=seq_path, taxon='eukaryote', rewrite_seq_path=False)[core.KnowledgeBase][0]

        tmp_core_path = os.path.join(self.dir, 'tmp_eukaryote_core.xlsx')
        tmp_seq_path = os.path.join(self.dir, 'tmp_eukaryote_seq.fna')

        writer = io.Writer()
        writer.run(tmp_core_path, kb, seq_path=tmp_seq_path, taxon='eukaryote', data_repo_metadata=False)

        tmp_kb = reader.run(tmp_core_path, seq_path, taxon='eukaryote')[core.KnowledgeBase][0]

        self.assertTrue(kb.is_equal(tmp_kb))
        self.assertTrue(filecmp.cmp(tmp_seq_path, seq_path, shallow=False))

    def test_rewrite_seq_path_in_read_write(self):
        path_core_1 = os.path.join(self.dir, 'core_1.xlsx')
        path_core_2 = os.path.join(self.dir, 'core_2.xlsx')
        path_seq_1 = os.path.join(self.dir, 'seq_1.fna')
        path_seq_2 = os.path.join(self.dir, 'seq_2.fna')

        io.Writer().run(path_core_1, self.kb, seq_path=path_seq_1, data_repo_metadata=False)
        kb1 = io.Reader().run(path_core_1, seq_path=path_seq_1)[core.KnowledgeBase][0]
        kb2 = io.Reader().run(path_core_1, seq_path=path_seq_1, rewrite_seq_path=False)[core.KnowledgeBase][0]
        kb3 = io.Reader().run(path_core_1, seq_path=self.seq_path)[core.KnowledgeBase][0]
        kb4 = io.Reader().run(path_core_1, seq_path=self.seq_path, rewrite_seq_path=False)[core.KnowledgeBase][0]
        self.assertFalse(kb1.is_equal(self.kb))
        self.assertFalse(kb2.is_equal(self.kb))
        self.assertTrue(kb3.is_equal(self.kb))
        self.assertFalse(kb4.is_equal(self.kb))
        self.assertTrue(filecmp.cmp(path_seq_1, self.seq_path, shallow=False))

        io.Writer().run(path_core_2, self.kb, seq_path=path_seq_2, rewrite_seq_path=False, data_repo_metadata=False)
        kb5 = io.Reader().run(path_core_2, seq_path=path_seq_2)[core.KnowledgeBase][0]
        kb6 = io.Reader().run(path_core_2, seq_path=path_seq_2, rewrite_seq_path=False)[core.KnowledgeBase][0]
        kb7 = io.Reader().run(path_core_2, seq_path=self.seq_path)[core.KnowledgeBase][0]
        kb8 = io.Reader().run(path_core_2, seq_path=self.seq_path, rewrite_seq_path=False)[core.KnowledgeBase][0]
        self.assertFalse(kb5.is_equal(self.kb))
        self.assertTrue(kb6.is_equal(self.kb))
        self.assertTrue(kb7.is_equal(self.kb))
        self.assertTrue(kb8.is_equal(self.kb))
        self.assertTrue(filecmp.cmp(path_seq_2, self.seq_path, shallow=False))

    def test_write_with_repo_metadata(self):

        with tempfile.TemporaryDirectory() as temp_dir:

            # create temp git repo & write file into it
            test_repo_name = 'test_wc_kb_test_io'
            test_github_repo = GitHubRepoForTests(test_repo_name)
            repo = test_github_repo.make_test_repo(temp_dir)

            _, core_path = tempfile.mkstemp(dir=temp_dir, suffix='.xlsx')
            _, seq_path = tempfile.mkstemp(dir=temp_dir, suffix='.fna')

            # write data repo metadata in data_file
            writer = io.Writer()
            writer.run(core_path, self.kb, seq_path=seq_path, data_repo_metadata=True)

            # deliberately read metadata
            reader = io.Reader()
            objs_read = reader.run(core_path, seq_path=seq_path, read_metadata=True)
            data_repo_metadata = objs_read[utils.DataRepoMetadata][0]
            self.assertTrue(data_repo_metadata.url.startswith('https://github.com/'))
            self.assertEqual(data_repo_metadata.branch, 'master')
            self.assertEqual(len(data_repo_metadata.revision), 40)

    def test_write_without_cell_relationships(self):
        core_path = os.path.join(self.dir, 'core.xlsx')
        seq_path = os.path.join(self.dir, 'test_seq.fna')

        with open(seq_path, 'w') as file:
            file.write('>chr_x\nACGT\n')

        dna = core.DnaSpeciesType(id='chr_x', sequence_path=seq_path)
        self.kb.cell.species_types.append(dna)

        trn = prokaryote.TranscriptionUnitLocus(id='tu_x_0')
        dna.loci.append(trn)
        trn.cell = None

        writer = io.Writer()
        with self.assertRaisesRegex(ValueError, 'must be set to the instance of `Cell`'):
            writer.run(core_path, self.kb, seq_path=seq_path, data_repo_metadata=False)

    def test_write_read_sloppy(self):
        core_path = os.path.join(self.dir, 'core.xlsx')
        seq_path = os.path.join(self.dir, 'test_seq.fna')

        writer = io.Writer()
        writer.run(core_path, self.kb, seq_path=seq_path, data_repo_metadata=False)

        wb = wc_utils.workbook.io.read(core_path)

        row = wb['!!KB'].pop(4)
        wb['!!KB'].insert(5, row)
        wc_utils.workbook.io.write(core_path, wb)

        reader = io.Reader()
        with self.assertRaisesRegex(ValueError, "cannot be loaded because"):
            reader.run(core_path, seq_path=self.seq_path)
        env = EnvironmentVarGuard()
        env.set('CONFIG__DOT__wc_kb__DOT__io__DOT__strict', '0')
        with env:
            kb = reader.run(core_path, self.seq_path)[core.KnowledgeBase][0]

        self.assertTrue(kb.is_equal(self.kb))
        self.assertTrue(filecmp.cmp(self.seq_path, seq_path, shallow=False))

    def test_reader_no_kb(self):
        core_path = os.path.join(self.dir, 'core.xlsx')
        obj_tables.io.WorkbookWriter().run(core_path, [], models=io.PROKARYOTE_MODELS, include_all_attributes=False)

        seq_path = os.path.join(self.dir, 'test_seq.fna')
        with open(seq_path, 'w') as file:
            pass

        with self.assertRaisesRegex(ValueError, 'should define one knowledge base'):
            io.Reader().run(core_path, seq_path=seq_path)

        obj_tables.io.WorkbookWriter().run(core_path, [core.Cell(id='cell')], models=io.PROKARYOTE_MODELS, include_all_attributes=False)
        with self.assertRaisesRegex(ValueError, 'should define one knowledge base'):
            io.Reader().run(core_path, seq_path=seq_path)

    def test_reader_error_multiple_kbs(self):
        kb1 = core.KnowledgeBase(id='kb1', name='kb1', version='0.0.1')
        kb2 = core.KnowledgeBase(id='kb2', name='kb2', version='0.0.1')

        core_path = os.path.join(self.dir, 'core.xlsx')
        obj_tables.io.WorkbookWriter().run(core_path, [kb1, kb2], models=io.PROKARYOTE_MODELS, include_all_attributes=False)

        seq_path = os.path.join(self.dir, 'test_seq.fna')
        with open(seq_path, 'w') as file:
            pass

        with self.assertRaisesRegex(ValueError, ' should define one knowledge base'):
            io.Reader().run(core_path, seq_path=seq_path)

    def test_reader_no_cell(self):
        kb = core.KnowledgeBase(id='kb', name='kb1', version='0.0.1')
        dna = core.DnaSpeciesType(id='chr')

        core_path = os.path.join(self.dir, 'core.xlsx')
        obj_tables.io.WorkbookWriter().run(core_path, [kb, dna], models=io.PROKARYOTE_MODELS, include_all_attributes=False)

        seq_path = os.path.join(self.dir, 'test_seq.fna')
        with open(seq_path, 'w') as file:
            pass

        io.Reader().run(core_path, seq_path=seq_path)

    def test_reader_error_multiple_cells(self):
        kb = core.KnowledgeBase(id='kb', name='kb1', version='0.0.1')
        cell1 = core.Cell(id='cell1', name='cell1')
        cell2 = core.Cell(id='cell2', name='cell2')

        core_path = os.path.join(self.dir, 'core.xlsx')
        obj_tables.io.WorkbookWriter().run(core_path, [kb, cell1, cell2], models=io.PROKARYOTE_MODELS, include_all_attributes=False)

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

        io.Writer().run(path_core_1, self.kb, seq_path=path_seq_1, data_repo_metadata=False)
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

        io.Writer().run(path_core_1, self.kb, seq_path=path_seq_1, data_repo_metadata=False)
        self.assertTrue(filecmp.cmp(path_seq_1, self.seq_path, shallow=False))

        wb = wc_utils.workbook.io.read(path_core_1)
        row = wb['!!KB'].pop(4)
        wb['!!KB'].insert(5, row)
        wc_utils.workbook.io.write(path_core_1, wb)

        with self.assertRaisesRegex(ValueError, "cannot be loaded because"):
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
        io.create_template(path_core, path_seq, data_repo_metadata=False)
        kb = io.Reader().run(path_core, seq_path=path_seq)[core.KnowledgeBase][0]

    def test_validate_implicit_relationships(self):
        class TestModel(obj_tables.Model):
            id = obj_tables.StringAttribute(primary=True, unique=True)

        try:
            core.KnowledgeBase.Meta.attributes['test'] = obj_tables.OneToOneAttribute(TestModel, related_name='a')
            with self.assertRaisesRegex(Exception, 'Relationships from `KnowledgeBase` not supported:'):
                io.Writer.validate_implicit_relationships()
        finally:
            core.KnowledgeBase.Meta.attributes.pop('test')

        try:
            core.KnowledgeBase.Meta.related_attributes['test'] = obj_tables.OneToManyAttribute(core.Cell, related_name='c')
            with self.assertRaisesRegex(Exception,
                                        'Relationships to `KnowledgeBase` that are not one-to-one are prohibited'):
                io.Writer.validate_implicit_relationships()
        finally:
            core.KnowledgeBase.Meta.related_attributes.pop('test')

        try:
            core.Cell.Meta.attributes['test'] = obj_tables.OneToManyAttribute(TestModel, related_name='c')
            with self.assertRaisesRegex(Exception,
                                        'Relationships from `Cell` to `KnowledgeBase` that are not one-to-one are prohibited:'):
                io.Writer.validate_implicit_relationships()
        finally:
            core.Cell.Meta.attributes.pop('test')

        try:
            core.Cell.Meta.attributes['test'] = obj_tables.OneToOneAttribute(TestModel, related_name='d')
            with self.assertRaisesRegex(Exception,
                                        'Relationships from `Cell` to classes other than `KnowledgeBase` are prohibited:'):
                io.Writer.validate_implicit_relationships()
        finally:
            core.Cell.Meta.attributes.pop('test')

        try:
            core.Cell.Meta.related_attributes['test'] = obj_tables.OneToManyAttribute(TestModel, related_name='d')
            with self.assertRaisesRegex(Exception,
                                        'Relationships to `Cell` that are not one-to-one or many-to-one are prohibited: '):
                io.Writer.validate_implicit_relationships()
        finally:
            core.Cell.Meta.related_attributes.pop('test')

        try:
            core.KnowledgeBase.Meta.related_attributes['test'] = obj_tables.OneToOneAttribute(TestModel, related_name='b')
            with self.assertRaisesRegex(Exception,
                                        'Relationships to `KnowledgeBase` from classes other than `Cell` are prohibited'):
                io.Writer.validate_implicit_relationships()
        finally:
            core.KnowledgeBase.Meta.related_attributes.pop('test')

    @unittest.skip('Assertions need to be uncommented')
    def test_seq_path_consistency(self):
        pass

        #core_path = os.path.join(self.dir, 'core.xlsx')
        #import pdb; pdb.set_trace()
        #kb = io.Reader().run(core_path, seq_path=self.seq_path)[core.KnowledgeBase][0]
        #core_path = os.path.join(self.dir, 'core2.xlsx')
        #seq_path = os.path.join(self.dir, 'seq2.fna')

    def test_read_flat_list_of_objects(self):
        core_path = os.path.join(self.dir, 'core.xlsx')

        writer = io.Writer()
        writer.run(core_path, self.kb, data_repo_metadata=False)

        reader = io.Reader()

        objs = reader.run(core_path, seq_path=self.seq_path)
        self.assertIsInstance(objs, dict)

        objs = reader.run(core_path, seq_path=self.seq_path, 
                          group_objects_by_model=False)
        self.assertIsInstance(objs, list)
        kb = next(obj for obj in objs if isinstance(obj, core.KnowledgeBase))
        self.assertTrue(kb.is_equal(self.kb))
