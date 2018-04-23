""" Tests of utilities.

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-02-12
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_kb import core
from wc_kb import util
import shutil
import tempfile
import unittest


class TestUtil(unittest.TestCase):
    """ Test utilities """

    def test_get_models(self):
        self.assertIn(core.KnowledgeBase, util.get_models())
        self.assertIn(core.Cell, util.get_models())

        self.assertIn(core.KnowledgeBase, util.get_models(inline=False))
        self.assertIn(core.Cell, util.get_models(inline=False))

    def test_set_git_repo_metadata_from_path(self):
        kb = core.KnowledgeBase()
        self.assertEqual(kb.url, '')

        util.set_git_repo_metadata_from_path(kb, path='.')
        self.assertIn(kb.url, [
            'https://github.com/KarrLab/wc_kb.git',
            'ssh://git@github.com/KarrLab/wc_kb.git',
            'git@github.com:KarrLab/wc_kb.git',
        ])

    def test_set_git_repo_metadata_from_path_error(self):
        tempdir = tempfile.mkdtemp()

        kb = core.KnowledgeBase()
        self.assertEqual(kb.url, '')

        with self.assertRaisesRegexp(ValueError, 'is not a Git repository'):
            util.set_git_repo_metadata_from_path(kb, path=tempdir)
        self.assertEqual(kb.url, '')

        shutil.rmtree(tempdir)
