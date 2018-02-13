""" Utilities

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-02-12
:Copyright: 2018, Karr Lab
:License: MIT
"""

from . import core
import obj_model.core

def get_models(inline=True):
    """ Get list of models

    Args:
        inline (:obj:`bool`, optional): if true, return inline models

    Returns:
        :obj:`list` of `class`: list of models
    """

    return obj_model.core.get_models(module=core, inline=inline)
