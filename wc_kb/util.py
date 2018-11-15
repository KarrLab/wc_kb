""" Utilities

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-02-12
:Copyright: 2018, Karr Lab
:License: MIT
"""

from . import core
from wc_utils.util import git
import obj_model.core


def get_models(inline=True):
    """ Get list of models

    Args:
        inline (:obj:`bool`, optional): if true, return inline models

    Returns:
        :obj:`list` of :obj:`class`: list of models
    """

    return obj_model.core.get_models(module=core, inline=inline)


def set_git_repo_metadata_from_path(kb, path='.'):
    """ Use Git to set the Git repository URL, branch, and revision metadata for a knowledge base

    Args:
        kb (:obj:`core.KnowledgeBase`): knowledge base
        path (:obj:`str`, optional): path to the Git repository for the knowledge base
    """
    md = git.get_repo_metadata(dirname=path)
    kb.url = md.url
    kb.branch = md.branch
    kb.revision = md.revision
