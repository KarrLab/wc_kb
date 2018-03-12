import pkg_resources

with open(pkg_resources.resource_filename('wc_kb', 'VERSION'), 'r') as file:
    __version__ = file.read().strip()
# :obj:`str`: version

# API
from .core import (PolymerStrand,
                   KnowledgeBaseObject, KnowledgeBase, Cell, Compartment,
                   SpeciesType, MetaboliteSpeciesType, PolymerSpeciesType,
                   DnaSpeciesType, RnaType, RnaSpeciesType, ProteinSpeciesType,
                   PolymerLocus, GeneType, GeneLocus, PromoterLocus, OpenReadingFrameLocus,
                   ReactionParticipant, Reaction)
from . import io
from . import util
