""" Command line programs for managing knowledge bases for whole-cell models

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-04-20
:Copyright: 2018, Karr Lab
:License: MIT
"""

from cement.core.foundation import CementApp
from cement.core.controller import CementBaseController, expose
from wc_kb import io
import wc_kb
import wc_utils.workbook.io


class BaseController(CementBaseController):
    """ Base controller for command line application """

    class Meta:
        label = 'base'
        description = "Command line programs for managing knowledge bases for whole-cell models"

    @expose(help='Get version')
    def get_version(self):
        """ Get version """
        print(wc_kb.__version__)


class ValidateController(CementBaseController):
    """ Validate knowledge base and display errors """

    class Meta:
        label = 'validate'
        description = 'Validate knowledge base and display errors'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['core_path'], dict(type=str, help='Path to knowledge base core')),
            (['seq_path'], dict(type=str, help='Path to FASTA-formatted genome sequence')),
            (['--sloppy'], dict(dest='strict', default=True, action='store_false',
                                help='If set, do not validate the format of the knowledge base core file(s)')),
        ]

    @expose(hide=True)
    def default(self):
        args = self.app.pargs
        try:
            io.Reader().run(args.core_path, args.seq_path, strict=args.strict)
            print('Knowledge base is valid')
        except ValueError as exception:
            raise ValueError('Knowledge base is invalid: ' + str(exception))


class DifferenceController(CementBaseController):
    """ Display difference between two knowledge bases """

    class Meta:
        label = 'difference'
        description = 'Get difference between two knowledge bases'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['core_path_1'], dict(type=str, help='Path to core for first knowledge base')),
            (['seq_path_1'], dict(type=str, help='Path to FASTA-formatted genome sequence for first knowledge base')),
            (['core_path_2'], dict(type=str, help='Path to core for second knowledge base')),
            (['seq_path_2'], dict(type=str, help='Path to FASTA-formatted genome sequence for second knowledge base')),
            (['--compare-files'], dict(dest='compare_files', default=False, action='store_true',
                                       help='If true, compare knowledge bases; otherwise compare files directly')),
            (['--sloppy'], dict(dest='strict', default=True, action='store_false',
                                help='If set, do not validate the format of the knowledge base file(s)')),
        ]

    @expose(hide=True)
    def default(self):
        args = self.app.pargs

        if args.compare_files:
            kb1 = wc_utils.workbook.io.read(args.core_path_1)
            kb2 = wc_utils.workbook.io.read(args.core_path_2)
            diff = kb1.difference(kb2)

        else:
            kb1 = io.Reader().run(args.core_path_1, args.seq_path_1, strict=args.strict)
            kb2 = io.Reader().run(args.core_path_2, args.seq_path_2, strict=args.strict)
            diff = kb1.difference(kb2)

        if diff:
            print(diff)
        else:
            print('Knowledge bases are identical')


class NormalizeController(CementBaseController):
    """ Normalize knowledge base """

    class Meta:
        label = 'normalize'
        description = 'Normalize knowledge base'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['source_core'], dict(
                type=str,
                help='Path to core of the knowledge base')),
            (['source_seq'], dict(
                type=str,
                help='Path to FASTA-formatted genome sequence for the knowledge base')),
            (['--dest-core'], dict(
                default='', type=str,
                help='Path to save normalized core of the knowledge base')),
            (['--dest-seq'], dict(
                default='', type=str,
                help='Path to save normalized FASTA-formatted genome sequence for the knowledge base')),
            (['--sloppy'], dict(
                dest='strict', default=True, action='store_false',
                help='If set, do not validate the format of the knowledge base file(s)')),
        ]

    @expose(hide=True)
    def default(self):
        args = self.app.pargs
        kb = io.Reader().run(args.source_core, args.source_seq, strict=args.strict)
        if args.dest_core or args.dest_seq:
            io.Writer().run(kb, args.dest_core, args.dest_seq)
        else:
            io.Writer().run(kb, args.source_core, args.source_seq)


class ConvertController(CementBaseController):
    """ Convert knowledge base among Excel (.xlsx), comma separated (.csv), JavaScript Object Notation (.json),
    tab separated (.tsv), and Yet Another Markup Language (.yaml, .yml) formats """

    class Meta:
        label = 'convert'
        description = 'Convert knowledge base among .csv, .json, .tsv, .xlsx, .yaml, and .yml formats'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['source_core'], dict(type=str, help='Path to core of the knowledge base')),
            (['source_seq'], dict(type=str, help='Path to FASTA-formatted genome sqeuence of the knowledge base')),
            (['dest_core'], dict(type=str, help='Path to save the converted core of the knowledge base')),
            (['dest_seq'], dict(type=str, help='Path to save the converted FASTA-formatted genome sequence of the knowledge base')),
            (['--sloppy'], dict(dest='strict', default=True, action='store_false',
                                help='If set, do not validate the format of the knowledge base file(s)')),
        ]

    @expose(hide=True)
    def default(self):
        args = self.app.pargs
        io.convert(args.source_core, args.source_seq, args.dest_core, args.dest_seq, strict=args.strict)


class CreateTemplateController(CementBaseController):
    """ Create file with knowledge base template (i.e. create file with row and column labels) """

    class Meta:
        label = 'create-template'
        description = 'Create file with knowledge base template: blank file(s) with row and column labels'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['path_core'], dict(type=str, help='Path to save a template of the core of a knowledge base')),
            (['path_seq'], dict(type=str, help='Path to save a template of the genome sequence of a knowledge base')),
        ]

    @expose(hide=True)
    def default(self):
        args = self.app.pargs
        io.create_template(args.path_core, args.path_seq)


class UpdateWcKbVersionController(CementBaseController):
    """ Update wc_kb_version of a knowledge base """

    class Meta:
        label = 'update-wc-kb-version'
        description = 'Update wc_kb_version of a knowledge base'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['path_core'], dict(type=str, help='Path to the core of the knowledge base')),
            (['path_seq'], dict(type=str, help='Path to the FASTA-formatted genome sequence of a knowledge base')),
            (['--sloppy'], dict(dest='strict', default=True, action='store_false',
                                help='If set, do not validate the format of the knowledge base file(s)')),
        ]

    @expose(hide=True)
    def default(self):
        args = self.app.pargs
        kb = io.Reader().run(args.path_core, args.path_seq, strict=args.strict)
        kb.wc_kb_version = wc_kb.__version__
        io.Writer().run(kb, args.path_core, args.path_seq)


class App(CementApp):
    """ Command line application """
    class Meta:
        label = 'wc_kb'
        base_controller = 'base'
        handlers = [
            BaseController,
            ValidateController,
            DifferenceController,
            NormalizeController,
            ConvertController,
            CreateTemplateController,
            UpdateWcKbVersionController,
        ]


def main():
    with App() as app:
        app.run()
