import pkg_resources

with open(pkg_resources.resource_filename('wc_kb', 'VERSION'), 'r') as file:
    __version__ = file.read().strip()
# :obj:`str`: version

# API
from .core import (KnowledgeBaseObject,
                   KnowledgeBase,
                   Cell,
                   Compartment,
                   SpeciesType,
                   MetaboliteSpeciesType,
                   PolymerSpeciesType,
                   ComplexSpeciesType,
                   DnaSpeciesType,
                   RnaType,
                   RnaSpeciesType,
                   ProteinSpeciesType,
                   PolymerLocus,
                   GeneType,
                   GeneLocus,
                   PromoterLocus,
                   TranscriptionUnitLocus,
                   Reaction,
                   Species,
                   SpeciesCoefficient,
                   ReactionParticipantAttribute,
                   SubunitAttribute,
                   Property,
                   Observable,
                   ObservableCoefficient,
                   ObservableObservableParticipantAttribute,
                   ObservableSpeciesParticipantAttribute,
                   PolymerStrand,
                   RateLawDirection,
                   RateLawEquationAttribute,
                   RateLaw,
                   RateLawEquation)

from . import io
from . import util
