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
                   SpeciesTypeCoefficient,
                   MetaboliteSpeciesType,
                   PolymerSpeciesType,
                   ComplexSpeciesType,
                   DnaSpeciesType,
                   RnaType,
                   PolymerLocus,
                   GeneType,
                   Reaction,
                   Species,
                   SpeciesCoefficient,
                   Concentration,
                   OneToOneSpeciesAttribute,
                   DatabaseReferenceAttribute,
                   ReactionParticipantAttribute,
                   SubunitAttribute,
                   Property,
                   Observable,
                   ObservableCoefficient,
                   ObservableObservableParticipantAttribute,
                   ObservableSpeciesParticipantAttribute,
                   PolymerStrand,
                   Parameter,
                   RateLawDirection,
                   RateLawEquationAttribute,
                   RateLaw,
                   RateLawEquation,
                   Reference,
                   DatabaseReference)

from . import io
from . import util
