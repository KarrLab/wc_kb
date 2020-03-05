""" Utilities

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-02-12
:Copyright: 2018, Karr Lab
:License: MIT
"""

from . import core
from wc_utils.util import git
import obj_tables


def get_models(inline=True):
    """ Get list of models

    Args:
        inline (:obj:`bool`, optional): if true, return inline models

    Returns:
        :obj:`list` of :obj:`class`: list of models
    """

    return obj_tables.get_models(module=core, inline=inline)
