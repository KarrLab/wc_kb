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
        non_inline_models = set([
            core.KnowledgeBase,
            core.Cell,
            core.SpeciesType,
            core.Chromosome,
            core.TranscriptionUnit,
            core.Reaction,
        ])
        inline_models = set()
        self.assertEqual(set(util.get_models()), non_inline_models | inline_models)
        self.assertEqual(set(util.get_models(inline=False)), non_inline_models)
