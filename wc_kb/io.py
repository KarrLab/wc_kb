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
from . import eukaryote
from . import prokaryote
from . import util
from wc_utils.util.string import indent_forest
import Bio.SeqIO
import Bio.SeqRecord
import obj_tables
import os
import shutil
import wc_kb
import wc_kb.config.core
import wc_utils.cache
import warnings

PROKARYOTE_MODELS = (
    core.KnowledgeBase,
    core.Cell,
    core.Compartment,
    core.DnaSpeciesType,
    core.ChromosomeFeature,
    prokaryote.TranscriptionUnitLocus,
    prokaryote.GeneLocus,
    prokaryote.RnaSpeciesType,
    prokaryote.ProteinSpeciesType,
    core.ComplexSpeciesType,
    core.MetaboliteSpeciesType,
    core.SpeciesTypeProperty,
    core.Concentration,
    core.Observable,
    core.Reaction,
    core.RateLaw,
    core.Parameter,
    core.Evidence,
    core.Experiment,
    core.Reference)

EUKARYOTE_MODELS = (
    core.KnowledgeBase,
    core.Cell,
    core.Compartment,
    core.DnaSpeciesType,
    eukaryote.GeneLocus,
    eukaryote.RegulatoryModule,
    eukaryote.TranscriptSpeciesType,
    eukaryote.ProteinSpeciesType,
    eukaryote.PtmSite,
    core.ComplexSpeciesType,
    core.MetaboliteSpeciesType,
    core.SpeciesTypeProperty,
    core.Concentration,
    core.Observable,
    core.Reaction,
    core.RateLaw,
    core.Parameter,
    core.Evidence,
    core.Experiment,
    core.Reference)

class Writer(obj_tables.io.Writer):
    """ Write knowledge base to file(s) """

    def run(self, core_path, knowledge_base,
            seq_path=None, rewrite_seq_path=True, taxon='prokaryote',
            models=None, get_related=True, include_all_attributes=False, validate=True,
            title=None, description=None, keywords=None, version=None, language=None, creator=None,
            extra_entries=0, data_repo_metadata=False, schema_package=None):
        """ Write knowledge base to file(s)

        Args:
            knowledge_base (:obj:`core.KnowledgeBase`): knowledge base
            core_path (:obj:`str`): path to save core knowledge base
            seq_path (:obj:`str`, optional): path to save genome sequence
            rewrite_seq_path (:obj:`bool`, optional): if :obj:`True`, the path to genome sequence in the saved knowledge base
                will be updated to the newly saved seq_path
            taxon (:obj:`str`, optional): type of model order to use
            models (:obj:`list` of :obj:`Model`, optional): models in the order that they should
                appear as worksheets; all models which are not in `models` will
                follow in alphabetical order
            get_related (:obj:`bool`, optional): if :obj:`True`, write object and all related objects
            include_all_attributes (:obj:`bool`, optional): if :obj:`True`, export all attributes including those
                not explictly included in `Model.Meta.attribute_order`
            validate (:obj:`bool`, optional): if :obj:`True`, validate the data
            title (:obj:`str`, optional): title
            description (:obj:`str`, optional): description
            keywords (:obj:`str`, optional): keywords
            version (:obj:`str`, optional): version
            language (:obj:`str`, optional): language
            creator (:obj:`str`, optional): creator
            extra_entries (:obj:`int`, optional): additional entries to display
            data_repo_metadata (:obj:`bool`, optional): if :obj:`True`, try to write metadata information
                about the file's Git repo; the repo must be current with origin, except for the file
            schema_package (:obj:`str`, optional): the package which defines the `obj_tables` schema
                used by the file; if not :obj:`None`, try to write metadata information about the
                the schema's Git repository: the repo must be current with origin

        Raises:
            :obj:`ValueError`: if any of the relationships with knowledge bases and cells are not set
        """
        if issubclass(self.get_writer(core_path), obj_tables.io.WorkbookWriter):
            self.validate_implicit_relationships()
            self.validate_implicit_relationships_are_set(knowledge_base)

        if taxon == 'prokaryote':
            models = PROKARYOTE_MODELS
        elif taxon == 'eukaryote':
            models = EUKARYOTE_MODELS

        # default metadata for exported file
        if title is None:
            title = knowledge_base.id
        if description is None:
            description = knowledge_base.name
        if version is None:
            version = knowledge_base.version
        if language is None:
            language = 'wc_kb'
        if creator is None:
            creator = '{}.{}'.format(self.__class__.__module__, self.__class__.__name__)

        # export sequences, if a path is provided
        if seq_path:
            dna_seqs = []
            original_seq_paths = []
            if knowledge_base.cell:
                dna_species_types = knowledge_base.cell.species_types.get(
                    __type=core.DnaSpeciesType)
                for species_type in dna_species_types:
                    dna_seqs.append(Bio.SeqRecord.SeqRecord(
                        species_type.get_seq(), species_type.id))
                    if rewrite_seq_path:
                        original_seq_paths.append((species_type, species_type.sequence_path))
                        species_type.sequence_path = seq_path

            with open(seq_path, 'w') as file:
                writer = Bio.SeqIO.FastaIO.FastaWriter(
                    file, wrap=70, record2title=lambda record: record.id)
                writer.write_file(dna_seqs)

            file.close()

        # export core
        super(Writer, self).run(core_path, knowledge_base, models=models, get_related=get_related,
                                include_all_attributes=include_all_attributes, validate=validate,
                                title=title, description=description, version=version, language=language,
                                creator=creator, extra_entries=extra_entries,
                                data_repo_metadata=data_repo_metadata, schema_package=schema_package)

        # reset sequence paths
        if seq_path and rewrite_seq_path:
            for species_type, path in original_seq_paths:
                species_type.sequence_path = path

    @classmethod
    def validate_implicit_relationships(cls):
        """ Check that relationships to :obj:`core.KnowledgeBase` and :obj:`core.Cell` do not need to be explicitly written to
        workbooks because they can be inferred by :obj:`Reader.run`

        Raises:
            :obj:`Exception`: if the Excel serialization involves an unsupported implicit relationship
        """
        for name, attr in core.KnowledgeBase.Meta.attributes.items():
            if isinstance(attr, obj_tables.RelatedAttribute):
                raise Exception(
                    "Relationships from `KnowledgeBase` not supported: {}.{} to {}".format(
                        'KnowledgeBase', name, attr.related_class.__name__))

        for name, attr in core.KnowledgeBase.Meta.related_attributes.items():
            if not isinstance(attr, obj_tables.OneToOneAttribute):
                raise Exception(
                    "Relationships to `KnowledgeBase` that are not one-to-one are prohibited: {}.{} to {}".format(
                        attr.related_class.__name__, name, 'KnowledgeBase'))

        for name, attr in core.Cell.Meta.attributes.items():
            if isinstance(attr, obj_tables.RelatedAttribute):
                if not isinstance(attr, obj_tables.OneToOneAttribute):
                    raise Exception(
                        "Relationships from `Cell` to `KnowledgeBase` that are not one-to-one are prohibited: {}.{} to {}".format(
                            'Cell', name, 'KnowledgeBase'))
                if attr.related_class != core.KnowledgeBase:
                    raise Exception(
                        "Relationships from `Cell` to classes other than `KnowledgeBase` are prohibited: {}.{} to {}".format(
                            'Cell', name, attr.related_class.__name__))

        for attr in core.Cell.Meta.related_attributes.values():
            if not isinstance(attr, (obj_tables.OneToOneAttribute, obj_tables.ManyToOneAttribute)):
                raise Exception(
                    "Relationships to `Cell` that are not one-to-one or many-to-one are prohibited: {}.{} to {}".format(
                        attr.related_class.__name__, attr.related_name, 'Cell'))

        for name, attr in core.KnowledgeBase.Meta.related_attributes.items():
            if attr.primary_class != core.Cell:
                raise Exception(
                    "Relationships to `KnowledgeBase` from classes other than `Cell` are prohibited: {}.{} to {}".format(
                        attr.related_class.__name__, name, 'KnowledgeBase'))

        return None  # pragma: no cover; avoids missing branch coverage on previous for loop

    def validate_implicit_relationships_are_set(self, knowledge_base):
        """ Check that there is only 1 :obj:`KnowledgeBase` and <= 1 :obj:`Cell` and that each relationship
        to :obj:`KnowledgeBase` and :obj:`Cell` is set. This is necessary to enable the :obj:`KnowledgeBase` and
        :obj:`Cell` relationships to be implicit in the Excel output and added by :obj:`Reader.run`

        Args:
            knowledge_base (:obj:`core.KnowledgeBase`): knowledge base

        Raises:
            :obj:`ValueError`: if there are multiple instances of :obj:`core.KnowledgeBase` in the object graph
        """
        cell = knowledge_base.cell

        for obj in knowledge_base.get_related():
            for attr in obj.Meta.attributes.values():
                if isinstance(attr, obj_tables.RelatedAttribute) and attr.related_class == core.Cell:
                    val = getattr(obj, attr.name)
                    if val is None or val != cell:
                        raise ValueError('{}.{} must be set to the instance of `Cell`'.format(
                            obj.__class__.__name__, attr.name))


class Reader(obj_tables.io.Reader):
    """ Read knowledge base from file(s) """

    #@wc_utils.cache.memoize(filename_args=[1, 2])
    def run(self, core_path,
            seq_path='', rewrite_seq_path=True, taxon='prokaryote',
            models=None, ignore_missing_models=None, ignore_extra_models=None, ignore_sheet_order=None,
            include_all_attributes=False, ignore_missing_attributes=None, ignore_extra_attributes=None, ignore_attribute_order=None,
            group_objects_by_model=True, validate=True, read_metadata=False):


                """ Read knowledge base from file(s)

                Args:
                    core_path (:obj:`str`): path to core knowledge base
                    seq_path (:obj:`str`): path to genome sequence
                    rewrite_seq_path (:obj:`bool`, optional): if :obj:`True`, the path to genome sequence in the knowledge base
                        will be updated to the provided seq_path
                    taxon (:obj:`str`, optional): type of model order to use
                    models (:obj:`types.TypeType` or :obj:`list` of :obj:`types.TypeType`, optional): type
                        of object to read or list of types of objects to read
                    ignore_missing_models (:obj:`bool`, optional): if :obj:`False`, report an error if a worksheet/
                        file is missing for one or more models
                    ignore_extra_models (:obj:`bool`, optional): if :obj:`True` and all `models` are found, ignore
                        other worksheets or files
                    ignore_sheet_order (:obj:`bool`, optional): if :obj:`True`, do not require the sheets to be provided
                        in the canonical order
                    include_all_attributes (:obj:`bool`, optional): if :obj:`True`, export all attributes including those
                        not explictly included in `Model.Meta.attribute_order`
                    ignore_missing_attributes (:obj:`bool`, optional): if :obj:`False`, report an error if a
                        worksheet/file doesn't contain all of attributes in a model in `models`
                    ignore_extra_attributes (:obj:`bool`, optional): if :obj:`True`, do not report errors if
                        attributes in the data are not in the model
                    ignore_attribute_order (:obj:`bool`): if :obj:`True`, do not require the attributes to be provided
                        in the canonical order
                    group_objects_by_model (:obj:`bool`, optional): if :obj:`True`, group decoded objects by their
                        types
                    validate (:obj:`bool`, optional): if :obj:`True`, validate the data
                    read_metadata (:obj:`bool`, optional): if :obj:`True`, read metadata models

                Returns:
                    :obj:`dict`: model objects grouped by `obj_tables.Model` class

                Raises:
                    :obj:`ValueError`: if :obj:`core_path`

                        * Defines multiple knowledge bases or cells
                        * Represents objects that cannot be linked to a knowledge base and/or cell
                """
                if issubclass(self.get_reader(core_path), obj_tables.io.WorkbookReader):
                    Writer.validate_implicit_relationships()

                if taxon == 'prokaryote':
                    models = PROKARYOTE_MODELS
                elif taxon == 'eukaryote':
                    models = EUKARYOTE_MODELS
                else:
                    raise ValueError('Unsupported taxon "{}"'.format(taxon))

                if read_metadata:
                    models = list(models) + [obj_tables.utils.DataRepoMetadata, obj_tables.utils.SchemaRepoMetadata]
                    ignore_missing_models = True
                    ignore_sheet_order = True

                config = wc_kb.config.core.get_config()['wc_kb']['io']
                if ignore_missing_models is None:
                    ignore_missing_models = not config['strict']
                if ignore_extra_models is None:
                    ignore_extra_models = not config['strict']
                if ignore_sheet_order is None:
                    ignore_sheet_order = not config['strict']
                if ignore_missing_attributes is None:
                    ignore_missing_attributes = not config['strict']
                if ignore_extra_attributes is None:
                    ignore_extra_attributes = not config['strict']
                if ignore_attribute_order is None:
                    ignore_attribute_order = not config['strict']

                # read core objects from file
                objects = super(Reader, self).run(core_path, models=models,
                                                  ignore_missing_models=ignore_missing_models,
                                                  ignore_extra_models=ignore_extra_models,
                                                  ignore_sheet_order=ignore_sheet_order,
                                                  include_all_attributes=include_all_attributes,
                                                  ignore_missing_attributes=ignore_missing_attributes,
                                                  ignore_extra_attributes=ignore_extra_attributes,
                                                  ignore_attribute_order=ignore_attribute_order,
                                                  group_objects_by_model=group_objects_by_model,
                                                  validate=False)

                # Check if sequence pathes are consistent
                for idx, chromosome in enumerate(objects[wc_kb.core.DnaSpeciesType]):

                    if (chromosome.sequence_path is None) or (chromosome.sequence_path==''):

                        chromosome.sequence_path = seq_path # Set seq_path to be what is provided to wc_kb.io.Reader()
                        if idx !=0:
                            warnings.warn('Same sequence file is associated with mulitple chromosomes, make sure seq file is formatted accordingly!')

                    else:

                        if chromosome.sequence_path != seq_path:
                            warnings.warn('Sequence path ({}) provided in KB file ({}) is different from \
                                           seq_path provided to wc_kb.io.Reader ({}).'.format(chromosome.sequence_path, core_path, seq_path))


                # check that file has 1 knowledge base
                if len(objects[core.KnowledgeBase]) != 1:
                    raise ValueError('"{}" should define one knowledge base'.format(core_path))
                kb = objects[core.KnowledgeBase][0]

                # check that file has 0 or 1 cells
                if not objects[core.Cell]:
                    cell = None
                elif len(objects[core.Cell]) == 1:
                    cell = objects[core.Cell][0]
                else:
                    raise ValueError('"{}" should define zero or one cells'.format(core_path))

                # add implict relationships to `KnowledgeBase` and `Cell`
                kb.cell = cell

                for model, model_objects in objects.items():
                    for attr in model.Meta.attributes.values():
                        if isinstance(attr, obj_tables.RelatedAttribute) and attr.related_class == core.Cell:
                            for model_obj in model_objects:
                                setattr(model_obj, attr.name, cell)

                # link path to genome sequence to the DNA species types if rewrite_seq_path is True
                if rewrite_seq_path:
                    for dna in Bio.SeqIO.parse(seq_path, "fasta"):
                        species_type = kb.cell.species_types.get_one(id=dna.id)
                        species_type.sequence_path = seq_path

                # validate
                config = wc_kb.config.core.get_config()['wc_kb']['io']
                if (validate is not None and validate) or (validate is None and config['validate']):
                    objs = []
                    for cls_objs in objects.values():
                        objs.extend(cls_objs)

                    errors = obj_tables.Validator().validate(objs)
                    if errors:
                        raise ValueError(
                            indent_forest(['The knowledge base cannot be loaded because it fails to validate:', [errors]]))

                return objects


def convert(source_core, source_seq, dest_core, dest_seq, rewrite_seq_path=True):
    """ Convert among Excel (.xlsx), comma separated (.csv), and tab separated (.tsv) file formats

    Read a knowledge base from the `source` files(s) and write it to the `destination` files(s). A path to a
    delimiter separated set of knowledge base files must be represented by a Unix glob pattern (with a \\*) that
    matches all delimiter separated files.

    Args:
        source_core (:obj:`str`): path to the core of the source knowledge base
        source_seq (:obj:`str`): path to the genome sequence of the source knowledge base
        dest_core (:obj:`str`): path to save the converted core of the knowledge base
        dest_seq (:obj:`str`): path to save the converted genome sequence of the knowledge base
        rewrite_seq_path (:obj:`bool`, optional): if :obj:`True`, the path to genome sequence in the converted
            core of the knowledge base will be updated to the path of the converted genome sequence
    """
    kb = Reader().run(source_core, seq_path=source_seq)[core.KnowledgeBase][0]
    Writer().run(dest_core, kb, seq_path=dest_seq, rewrite_seq_path=rewrite_seq_path, data_repo_metadata=False)


def create_template(core_path, seq_path, extra_entries=10, data_repo_metadata=True):
    """ Create file with knowledge base template, including row and column headings

    Args:
        core_path (:obj:`str`): path to save template of core knowledge base
        seq_path (:obj:`str`): path to save genome sequence
        extra_entries (:obj:`int`, optional): additional entries to display
        data_repo_metadata (:obj:`bool`, optional): if :obj:`True`, try to write metadata information
            about the file's Git repo
    """
    kb = core.KnowledgeBase(
        id='template', name='Template', version=wc_kb.__version__)
    Writer().run(core_path, kb, seq_path=seq_path,
                 extra_entries=extra_entries,
                 data_repo_metadata=data_repo_metadata)
