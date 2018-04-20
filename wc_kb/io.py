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
import Bio.SeqIO
import Bio.SeqRecord
import obj_model
import os
import shutil
import wc_utils.cache
import wc_kb


class Writer(object):
    """ Write knowledge base to file(s) """

    model_order = (
        core.KnowledgeBase,
        core.Cell,
        core.Compartment,
        core.MetaboliteSpeciesType,
        core.DnaSpeciesType,
        core.PromoterLocus,
        core.TranscriptionUnitLocus,
        core.RnaSpeciesType,
        core.GeneLocus,
        core.ProteinSpeciesType,
        core.ComplexSpeciesType,
        core.Reaction,
    )

    def run(self, knowledge_base, core_path, seq_path):
        """ Write knowledge base to file(s)

        Args:
            knowledge_base (:obj:`core.KnowledgeBase`): knowledge base
            core_path (:obj:`str`): path to save core knowledge base
            seq_path (:obj:`str`): path to save genome sequence
        """
        self.validate_implicit_kb_and_cell_relationships()

        cell = knowledge_base.cell

        # check that there is only 1 :obj:`KnowledgeBase` and only 1 :obj:`Cell` and that each relationship
        # to :obj:`KnowledgeBase` and :obj:`Cell` is set. This is necessary to enable the :obj:`KnowledgeBase` and
        # :obj:`Cell` relationships to be implicit in the Excel output and added by :obj:`Reader.run`
        # todo

        # gather DNA sequences
        dna_seqs = []
        if cell:
            dna_species_types = cell.species_types.get(__type=core.DnaSpeciesType)
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

        _, ext = os.path.splitext(core_path)
        writer = obj_model.io.get_writer(ext)()
        writer.run(core_path, objects, models=self.model_order, include_all_attributes=False, **kwargs)

        # export sequences
        with open(seq_path, 'w') as file:
            writer = Bio.SeqIO.FastaIO.FastaWriter(file, wrap=70, record2title=lambda record: record.id)
            writer.write_file(dna_seqs)

        # restore DNA sequences
        if cell:
            for species_type, seq in zip(dna_species_types, dna_seqs):
                species_type.seq = seq.seq

    @classmethod
    def validate_implicit_kb_and_cell_relationships(cls):
        """ Check that relationships to :obj:`core.KnowledgeBase` and :obj:`core.Cell` do not need to be explicitly written to 
        workbooks because they can be inferred by :obj:`Reader.run`
        """
        for attr in core.KnowledgeBase.Meta.attributes.values():
            if isinstance(attr, obj_model.RelatedAttribute):
                raise Exception('Relationships from `KnowledgeBase` not supported')

        for attr in core.KnowledgeBase.Meta.related_attributes.values():
            if attr.primary_class != core.Cell or not isinstance(attr, obj_model.OneToOneAttribute):
                raise Exception('Only one-to-one relationships to `KnowledgeBase` from `Cell` are supported')

        for attr in core.Cell.Meta.attributes.values():
            if isinstance(attr, obj_model.RelatedAttribute) and \
                    (not isinstance(attr, obj_model.OneToOneAttribute) or attr.related_class != core.KnowledgeBase):
                raise Exception('Only one-to-one relationships from `Cell` to `KnowledgeBase` are supported')

        for attr in core.Cell.Meta.related_attributes.values():
            if not isinstance(attr, (obj_model.OneToOneAttribute, obj_model.ManyToOneAttribute)):
                raise Exception('Only one-to-one and many-to-one relationships are supported to `Cell`')


class Reader(object):
    """ Read knowledge base from file(s) """

    @wc_utils.cache.memoize(filename_args=[1, 2])
    def run(self, core_path, seq_path, strict=True):
        """ Read knowledge base from file(s)

        Args:
            core_path (:obj:`str`): path to core knowledge base
            seq_path (:obj:`str`): path to genome sequence
            strict (:obj:`str`, optional): if :obj:`True`, validate that the the model file(s) strictly follow the
                :obj:`obj_model` serialization format:

                * The worksheets are in the expected order
                * There are no missing worksheets
                * There are no extra worksheets
                * The columns are in the expected order
                * There are no missing columns
                * There are no extra columns

        Returns:
            :obj:`core.KnowledgeBase`: knowledge base

        Raises:
            :obj:`ValueError`: if :obj:`core_path` defines multiple knowledge bases
        """
        Writer.validate_implicit_kb_and_cell_relationships()

        # read core objects from file
        _, ext = os.path.splitext(core_path)
        reader = obj_model.io.get_reader(ext)()

        kwargs = {}
        if isinstance(reader, obj_model.io.WorkbookReader) and not strict:
            kwargs['ignore_missing_sheets'] = True
            kwargs['ignore_extra_sheets'] = True
            kwargs['ignore_sheet_order'] = True
            kwargs['ignore_missing_attributes'] = True
            kwargs['ignore_extra_attributes'] = True
            kwargs['ignore_attribute_order'] = True

        objects = reader.run(core_path, models=Writer.model_order, include_all_attributes=False, **kwargs)

        # check that file has 0 or 1 knowledge bases
        if not objects[core.KnowledgeBase]:
            for model, model_objects in objects.items():
                if model_objects:
                    raise ValueError('"{}" cannot contain instances of `{}` without an instance of `KnowledgeBase`'.format(
                        core_path, model.__name__))
            return None

        elif len(objects[core.KnowledgeBase]) > 1:
            raise ValueError('"{}" should define one knowledge base'.format(core_path))

        else:
            kb = objects[core.KnowledgeBase].pop()

        # check that file has 0 or 1 cells
        if not objects[core.Cell]:
            for model, model_objects in objects.items():
                if model_objects:
                    raise ValueError('"{}" cannot contain instances of `{}` without an instance of `Cell`'.format(
                        core_path, model.__name__))
            cell = None

        elif len(objects[core.Cell]) > 1:
            raise ValueError('"{}" should define one cell'.format(core_path))

        else:
            cell = objects[core.Cell].pop()

        # add implict relationships to `KnowledgeBase` and `Cell`
        kb.cell = cell

        for model, model_objects in objects.items():
            for attr in model.Meta.attributes.values():
                if isinstance(attr, obj_model.RelatedAttribute) and attr.related_class == core.Cell:
                    for model_obj in model_objects:
                        setattr(model_obj, attr.name, cell)

        # read genome sequence and link to the DNA species types
        for dna in Bio.SeqIO.parse(seq_path, "fasta"):
            kb.cell.species_types.get_one(id=dna.id).seq = dna.seq

        return kb


def convert(source_core, source_seq, dest_core, dest_seq, strict=True):
    """ Convert among Excel (.xlsx), comma separated (.csv), and tab separated (.tsv) file formats

    Read a knowledge base from the `source` files(s) and write it to the `destination` files(s). A path to a
    delimiter separated set of knowledge base files must be represented by a Unix glob pattern (with a \*) that
    matches all delimiter separated files.

    Args:
        source_core (:obj:`str`): path to the core of the source knowledge base
        source_seq (:obj:`str`): path to the genome sequence of the source knowledge base
        dest_core (:obj:`str`): path to save the converted core of the knowledge base
        dest_seq (:obj:`str`): path to save the converted genome sequence of the knowledge base
        strict (:obj:`str`, optional): if :obj:`True`, validate that the the model file(s) strictly follow the
                :obj:`obj_model` serialization format:

                * The worksheets are in the expected order
                * There are no missing worksheets
                * There are no extra worksheets
                * The columns are in the expected order
                * There are no missing columns
                * There are no extra columns
    """
    kb = Reader().run(source_core, source_seq, strict=strict)
    Writer().run(kb, dest_core, dest_seq)


def create_template(core_path, seq_path):
    """ Create file with knowledge base template, including row and column headings

    Args:
        core_path (:obj:`str`): path to save temploate of core knowledge base
        seq_path (:obj:`str`): path to save genome sequence
    """
    kb = core.KnowledgeBase(id='template', name='Template', version=wc_kb.__version__)
    Writer().run(kb, core_path, seq_path)
