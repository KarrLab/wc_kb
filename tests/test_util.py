""" Tests of utilities.

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-02-12
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_kb import core
from wc_kb import util
import unittest


class TestUtil(unittest.TestCase):
    """ Test utilities """

    def test_get_models(self):
        self.assertIn(core.KnowledgeBase, util.get_models())
        self.assertIn(core.Cell, util.get_models())

        self.assertIn(core.KnowledgeBase, util.get_models(inline=False))
        self.assertIn(core.Cell, util.get_models(inline=False))
