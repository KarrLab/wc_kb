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
                   ChromosomeFeature,
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
                   Observable,
                   ObservableExpression,
                   PolymerStrand,
                   Parameter,
                   RateLawDirection,
                   RateLaw,
                   RateLawExpression,
                   Reference,
                   DatabaseReference,
                   Evidence,
                   Experiment,
                   SpeciesTypeProperty)

from . import io
from . import util
from . import kb_transformer
