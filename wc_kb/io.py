""" Reading and writing knowledge bases to/from files.

Supported file types:

* Comma separated values (.csv)
* Excel (.xlsx)
* Tab separated values (.tsv)

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-02-12
:Copyright: 2018, Karr Lab
:License: MIT
"""

from . import core
from . import util
from obj_model import io
import Bio.SeqIO
import Bio.SeqRecord
import wc_kb


class Writer(object):
    """ Write knowledge base to file(s) """

    model_order = (
        core.KnowledgeBase, core.Cell,
        core.DnaSpeciesType, core.RnaSpeciesType, core.ProteinSpeciesType,
        core.GeneLocus, core.PromoterLocus, core.OpenReadingFrameLocus,
    )

    def run(self, knowledge_base, core_path, seq_path):
        """ Write knowledge base to file(s)

        Args:
            knowledge_base (:obj:`core.KnowledgeBase`): knowledge base
            core_path (:obj:`str`): path to save core knowledge base
            seq_path (:obj:`str`): path to save genome sequence
        """
        cell = knowledge_base.cell

        # gather DNA sequences
        dna_seqs = []
        if cell:
            dna_species_types = cell.get_dna_species_types()
            for species_type in dna_species_types:
                dna_seqs.append(Bio.SeqRecord.SeqRecord(species_type.seq, species_type.id))
                species_type.seq = None

        # export core
        kwargs = {
            'language': 'wc_kb',
            'creator': '{}.{}'.format(self.__class__.__module__, self.__class__.__name__),
        }
        objects = [knowledge_base]
        kwargs['title'] = knowledge_base.id
        kwargs['description'] = knowledge_base.name
        kwargs['version'] = knowledge_base.version

        io.Writer().run(core_path, objects, self.model_order, **kwargs)

        # export sequences
        with open(seq_path, 'w') as file:
            writer = Bio.SeqIO.FastaIO.FastaWriter(file, wrap=70, record2title=lambda record: record.id)
            writer.write_file(dna_seqs)

        # restore DNA sequences
        if cell:
            for species_type, seq in zip(dna_species_types, dna_seqs):
                species_type.seq = seq.seq


class Reader(object):
    """ Read knowledge base from file(s) """

    def run(self, core_path, seq_path):
        """ Read knowledge base from file(s)

        Args:
            core_path (:obj:`str`): path to core knowledge base
            seq_path (:obj:`str`): path to genome sequence

        Returns:
            :obj:`core.KnowledgeBase`: knowledge base

        Raises:
            :obj:`ValueError`: if :obj:`core_path` defines multiple knowledge bases
        """
        objects = io.Reader().run(core_path, util.get_models(inline=False))

        if not objects[core.KnowledgeBase]:
            return None

        if len(objects[core.KnowledgeBase]) > 1:
            raise ValueError('Knowledge base file "{}" should only define one knowledge base'.format(core_path))

        kb = objects[core.KnowledgeBase].pop()

        for dna in Bio.SeqIO.parse(seq_path, "fasta"):
            kb.cell.species_types.get_one(id=dna.id).seq = dna.seq

        return kb


def convert(core_source, core_destination):
    """ Convert among Excel (.xlsx), comma separated (.csv), and tab separated (.tsv) file formats

    Read a knowledge base from the `source` files(s) and write it to the `destination` files(s). A path to a
    delimiter separated set of knowledge base files must be represented by a Unix glob pattern (with a \*) that
    matches all delimiter separated files.

    Args:
        core_source (:obj:`str`): path to source core knowledge base
        core_destination (:obj:`str`): path to save converted core knowledge base
    """
    io.convert(core_source, core_destination, models=Writer.model_order)


def create_template(core_path, seq_path):
    """ Create file with knowledge base template, including row and column headings

    Args:
        core_path (:obj:`str`): path to save temploate of core knowledge base
        seq_path (:obj:`str`): path to save genome sequence
    """
    kb = core.KnowledgeBase(id='template', name='Template', version=wc_kb.__version__)
    Writer().run(kb, core_path, seq_path)
