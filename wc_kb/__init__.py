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
                   ChromosomeFeature,
                   PolymerLocus,
                   Reaction,
                   Species,
                   SpeciesCoefficient,
                   Concentration,
                   OneToOneSpeciesAttribute,
                   IdentifierAttribute,
                   ReactionParticipantAttribute,
                   SubunitAttribute,
                   Observable,
                   ObservableExpression,
                   Parameter,
                   RateLawDirection,
                   RateLaw,
                   RateLawExpression,
                   Reference,
                   Evidence,
                   Experiment,
                   SpeciesTypeProperty)

from . import io
from . import util
from . import prokaryote
from . import eukaryote
