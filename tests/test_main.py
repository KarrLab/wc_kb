""" Tests of command line program

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-04-20
:Copyright: 2018, Karr Lab
:License: MIT
"""

from capturer import CaptureOutput
from obj_tables import Validator
from os import path
from shutil import rmtree
from tempfile import mkdtemp
from wc_kb import __main__
from wc_kb import io
import mock
import unittest
import wc_kb
import wc_kb.core


class TestCli(unittest.TestCase):

    def setUp(self):
        self.tempdir = mkdtemp()

    def tearDown(self):
        rmtree(self.tempdir)

    def test_get_version(self):
        with CaptureOutput() as capturer:
            with __main__.App(argv=['-v']) as app:
                with self.assertRaises(SystemExit):
                    app.run()
                self.assertEqual(capturer.get_text(), wc_kb.__version__)

        with CaptureOutput() as capturer:
            with __main__.App(argv=['--version']) as app:
                with self.assertRaises(SystemExit):
                    app.run()
                self.assertEqual(capturer.get_text(), wc_kb.__version__)

    def test_validate(self):
        kb = wc_kb.core.KnowledgeBase(id='kb', name='KB', version='0.0.1a', wc_kb_version='0.0.1')
        self.assertEqual(Validator().run(kb, get_related=True), None)
        filename_core = path.join(self.tempdir, 'core.xlsx')
        filename_seq = path.join(self.tempdir, 'seq.fna')
        io.Writer().run(filename_core, kb, seq_path=filename_seq, data_repo_metadata=False)

        with CaptureOutput() as capturer:
            with __main__.App(argv=['validate', filename_core, filename_seq]) as app:
                app.run()
            self.assertEqual(capturer.get_text(), 'Knowledge base is valid')

    def test_validate_exception(self):
        kb = wc_kb.core.KnowledgeBase(id='kb', name='KB', version='0.0.1a', wc_kb_version='0.0.1')
        kb.cell = wc_kb.core.Cell(id='cell')
        kb.cell.compartments.create(id='c')
        kb.cell.compartments.create(id='c')

        self.assertNotEqual(Validator().run(kb, get_related=True), None)
        filename_core = path.join(self.tempdir, 'core.xlsx')
        filename_seq = path.join(self.tempdir, 'seq.fna')
        io.Writer().run(filename_core, kb, seq_path=filename_seq, data_repo_metadata=False)

        with self.assertRaisesRegex(SystemExit, '^Knowledge base is invalid: '):
            with __main__.App(argv=['validate', filename_core, filename_seq]) as app:
                app.run()

    def test_difference(self):
        kb1 = wc_kb.core.KnowledgeBase(id='kb', name='KB', version='0.0.1a', wc_kb_version='0.0.0')
        filename_core_1 = path.join(self.tempdir, 'core1.xlsx')
        filename_seq_1 = path.join(self.tempdir, 'seq1.fna')
        io.Writer().run(filename_core_1, kb1, seq_path=filename_seq_1, data_repo_metadata=False)

        kb2 = wc_kb.core.KnowledgeBase(id='kb', name='KB', version='0.0.1a', wc_kb_version='0.0.0')
        filename_core_2 = path.join(self.tempdir, 'core2.xlsx')
        filename_seq_2 = path.join(self.tempdir, 'seq2.fna')
        io.Writer().run(filename_core_2, kb2, seq_path=filename_seq_2, data_repo_metadata=False)

        kb3 = wc_kb.core.KnowledgeBase(id='kb', name='KB', version='0.0.1a', wc_kb_version='0.0.1')
        filename_core_3 = path.join(self.tempdir, 'core3.xlsx')
        filename_seq_3 = path.join(self.tempdir, 'seq3.fna')
        io.Writer().run(filename_core_3, kb3, seq_path=filename_seq_3, data_repo_metadata=False)

        with CaptureOutput(termination_delay=0.1) as capturer:
            with __main__.App(argv=['difference',
                                    filename_core_1, filename_seq_1,
                                    filename_core_2, filename_seq_2,
                                    ]) as app:
                app.run()
            self.assertEqual(capturer.get_text(), 'Knowledge bases are identical')

        with CaptureOutput() as capturer:
            with __main__.App(argv=['difference',
                                    filename_core_1, filename_seq_1,
                                    filename_core_2, filename_seq_2,
                                    '--compare-files']) as app:
                app.run()
            self.assertEqual(capturer.get_text(), 'Knowledge bases are identical')

        with CaptureOutput() as capturer:
            with __main__.App(argv=['difference',
                                    filename_core_1, filename_seq_1,
                                    filename_core_3, filename_seq_3,
                                    ]) as app:
                app.run()
            diff = ('Objects (KnowledgeBase: "kb", KnowledgeBase: "kb") have different attribute values:\n  '
                    '`wc_kb_version` are not equal:\n    0.0.0 != 0.0.1')
            self.assertEqual(capturer.get_text(), diff)

        with CaptureOutput() as capturer:
            with __main__.App(argv=['difference',
                                    filename_core_1, filename_seq_1,
                                    filename_core_3, filename_seq_3,
                                    '--compare-files']) as app:
                app.run()
            diff = 'Sheet KB:\n  Row 8:\n    Cell B: 0.0.0 != 0.0.1'
            #diff = 'Sheet KBnowledge base:\n  Row 11:\n    Cell B: 0.0.0 != 0.0.1'
            self.assertEqual(capturer.get_text(), diff)

    def test_normalize(self):
        filename_core_1 = path.join(self.tempdir, 'model-1.xlsx')
        filename_seq_1 = path.join(self.tempdir, 'seq-1.fna')
        filename_core_2 = path.join(self.tempdir, 'model-2.xlsx')
        filename_seq_2 = path.join(self.tempdir, 'seq-2.fna')

        kb = wc_kb.core.KnowledgeBase(id='kb', name='KB', version='0.0.1a', wc_kb_version='0.0.0')
        io.Writer().run(filename_core_1, kb, seq_path=filename_seq_1, data_repo_metadata=False)

        # with same dest
        with __main__.App(argv=['normalize', filename_core_1, filename_seq_1]) as app:
            app.run()

        kb2 = io.Reader().run(filename_core_1, seq_path=filename_seq_1)[wc_kb.core.KnowledgeBase][0]
        self.assertTrue(kb2.is_equal(kb))

        # with different dest
        with __main__.App(argv=['normalize', filename_core_1, filename_seq_1,
                                '--dest-core', filename_core_2, '--dest-seq', filename_seq_2]) as app:
            app.run()

        kb2 = io.Reader().run(filename_core_2, seq_path=filename_seq_2)[wc_kb.core.KnowledgeBase][0]
        self.assertTrue(kb2.is_equal(kb))

    def test_convert(self):
        filename_in_core = path.join(self.tempdir, 'in.core.xlsx')
        filename_in_seq = path.join(self.tempdir, 'in.seq.fna')
        filename_out_core = path.join(self.tempdir, 'out.core-*.csv')
        filename_out_seq = path.join(self.tempdir, 'out.seq.fna')

        kb = wc_kb.core.KnowledgeBase(id='kb', name='KB', version='0.0.1a', wc_kb_version='0.0.0')
        io.Writer().run(filename_in_core, kb, seq_path=filename_in_seq, data_repo_metadata=False)

        with __main__.App(argv=['convert',
                                filename_in_core, filename_in_seq,
                                filename_out_core, filename_out_seq,
                                ]) as app:
            app.run()

        self.assertTrue(path.isfile(path.join(self.tempdir, 'out.core-KB.csv')))

    def test_create_template(self):
        filename_core = path.join(self.tempdir, 'core.xlsx')
        filename_seq = path.join(self.tempdir, 'seq.fna')

        with __main__.App(argv=['create-template', filename_core, filename_seq, '--ignore-repo-metadata']) as app:
            app.run()

        self.assertTrue(path.isfile(filename_core))
        self.assertTrue(path.isfile(filename_seq))

    def test_update_version_metadata(self):
        filename_core = path.join(self.tempdir, 'core.xlsx')
        filename_seq = path.join(self.tempdir, 'seq.fna')

        kb = wc_kb.core.KnowledgeBase(id='kb', name='KB', version='0.0.1a', wc_kb_version='0.0.0')
        self.assertNotEqual(kb.wc_kb_version, wc_kb.__version__)
        io.Writer().run(filename_core, kb, seq_path=filename_seq, data_repo_metadata=False)

        with __main__.App(argv=['update-version-metadata', filename_core, filename_seq, '--ignore-repo-metadata']) as app:
            app.run()

        kb = io.Reader().run(filename_core, seq_path=filename_seq)[wc_kb.core.KnowledgeBase][0]
        self.assertEqual(kb.wc_kb_version, wc_kb.__version__)

    def test_raw_cli(self):
        with mock.patch('sys.argv', ['wc-kb', '--help']):
            with self.assertRaises(SystemExit) as context:
                __main__.main()
                self.assertRegex(context.Exception, 'usage: wc-kb')

        with mock.patch('sys.argv', ['wc-kb']):
            with CaptureOutput() as capturer:
                __main__.main()
                self.assertRegex(capturer.get_text(), 'usage: wc-kb')
