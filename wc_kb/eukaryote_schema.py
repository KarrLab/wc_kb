""" Schema to represent a knowledge base to build models of eukaryotes

:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2018-09-10
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_kb import core
from wc_onto import kb_onto as kbOnt
from wc_utils.util import chem
from wc_utils.util.ontology import are_terms_equivalent
import Bio.Alphabet
import Bio.Seq
import enum
import obj_model
import re


#####################
#####################
# Attributes

class LocusAttribute(obj_model.ManyToManyAttribute):
    """ Start and end coordinates attribute """

    def __init__(self, related_name='', verbose_name='', verbose_related_name='', help=''):
        """
        Args:
            related_name (:obj:`str`, optional): name of related attribute on `related_class`
            verbose_name (:obj:`str`, optional): verbose name
            verbose_related_name (:obj:`str`, optional): verbose related name
            help (:obj:`str`, optional): help message
        """
        super(LocusAttribute, self).__init__(GenericLocus,
                                            related_name=related_name, min_related=0, min_related_rev=0,
                                            verbose_name=verbose_name, verbose_related_name=verbose_related_name, help=help)

    def serialize(self, coordinates, encoded=None):
        """ Serialize related object
        Args:
            coordinates (:obj:`list` of :obj:`Model`): a list of instances of GenericLocus Python representation
            encoded (:obj:`dict`, optional): dictionary of objects that have already been encoded
        Returns:
            :obj:`str`: simple Python representation
        """
        if not coordinates:
            return ''

        return ', '.join(obj_model.serialize() for obj_model in coordinates)

    def deserialize(self, value, objects, decoded=None):
        """ Deserialize value
        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model
            decoded (:obj:`dict`, optional): dictionary of objects that have already been decoded
        Returns:
            :obj:`tuple` of :obj:`list` of :obj:`GenericLocus`, :obj:`InvalidAttribute` or :obj:`None`: :obj:`tuple` of cleaned value
                and cleaning error
        """
        if not value:
            return ([], None)

        pattern = r'([0-9]+)\:([0-9]+)'
        if not re.match(pattern, value, flags=re.I):
            return (None, obj_model.InvalidAttribute(self, ['Incorrectly formatted list of coordinates: {}'.format(value)]))

        objs = []
        for pat_match in re.findall(pattern, value, flags=re.I):
            start = int(pat_match[0])
            end = int(pat_match[1])
            if self.related_class not in objects:
                objects[self.related_class] = {}
            serialized_value = self.related_class()._serialize(start, end)
            if serialized_value in objects[self.related_class]:
                obj = objects[self.related_class][serialized_value]
            else:
                obj = self.related_class(start=start, end=end)
                objects[self.related_class][serialized_value] = obj
            objs.append(obj)
        return (objs, None)


#####################
#####################
# Locus types


class GeneLocus(core.PolymerLocus):
    """ Knowledge of a gene

    Attributes:
        symbol (:obj:`str`): symbol
        type (:obj:`GeneType`): type of gene

    Related attributes:
        transcripts (:obj:`list` of :obj:`TranscriptSpeciesType`): transcripts
        regulatory_modules (:obj:`list` of `RegulatoryModule`): regulatory_modules
    """
    symbol = obj_model.StringAttribute()
    homologs = obj_model.LongStringAttribute()

    class Meta(obj_model.Model.Meta):
        verbose_name = 'Gene'
        attribute_order = ('id', 'name', 'synonyms', 'symbol', 'homologs', 'polymer', 'strand', 'start',
                           'end', 'identifiers', 'references', 'comments')


class RegulatoryElementLocus(core.PolymerLocus):
    """ Knowledge of a regulatory element of a gene

    Attributes:
        type (:obj:`RegulatoryElementType`): type of regulatory element
        activity (:obj:`ActivityLevel`): cell-type specific activity level
        bound_start (:obj:`int`): start coordinate of binding
        bound_end (:obj:`int`): end coordinate of binding
        motif_features (:obj:`list` of :obj:`ProteinSpeciesType`): proteins that bind to the site

    Related attributes:
        regulatory_modules (:obj:`list` of :obj:`RegulatoryModule`): regulatory_modules
    """
    bound_start = obj_model.PositiveIntegerAttribute()
    bound_end = obj_model.PositiveIntegerAttribute()
    motif_features = obj_model.ManyToManyAttribute('ProteinSpeciesType', related_name='regulatory_elements')
    type = obj_model.ontology.OntologyAttribute(kbOnt,
                                  terms = kbOnt['RegulatoryElementType'].rchildren(),
                                  none=True)
    activity = obj_model.ontology.OntologyAttribute(kbOnt,
                                  terms = kbOnt['ActivityLevelType'].rchildren(),
                                  none=True)

    class Meta(obj_model.Model.Meta):
        verbose_name = 'Regulatory element'
        attribute_order = ('id', 'name', 'synonyms', 'type', 'activity', 'polymer', 'strand', 'start', 'end',
                           'bound_start', 'bound_end', 'motif_features', 'identifiers', 'references', 'comments')


class RegulatoryModule(obj_model.Model):
    """ Knowledge about regulatory modules

    Attributes:
        id (:obj:`str`): identifier
        name (:obj:`str`): name
        gene (:obj:`GeneLocus`): gene
        regulatory_element (:obj:`RegulatoryElementLocus`): regulatory element
        binding_factor (:obj:`ProteinSpeciesType`): binding factor
        type (:obj:`RegulationType`): type of regulation (proximal or distal)
        direction (:obj:`RegulatoryDirection`): direction of regulation
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references
        identifiers (:obj:`list` of :obj:`Identifier`): identifiers
    """
    id = obj_model.SlugAttribute(primary=True, unique=True)
    name = obj_model.StringAttribute()
    gene = obj_model.ManyToOneAttribute(GeneLocus, related_name='regulatory_modules')

    binding_factor = obj_model.ManyToManyAttribute('ProteinSpeciesType', related_name='regulatory_modules')
    comments = obj_model.LongStringAttribute()
    references = obj_model.ManyToManyAttribute(core.Reference, related_name='regulatory_modules')
    identifiers = core.IdentifierAttribute(related_name='regulatory_modules')
    regulatory_element = obj_model.ManyToOneAttribute(RegulatoryElementLocus, related_name='regulatory_modules')
    type = obj_model.ontology.OntologyAttribute(kbOnt,
                                  terms = kbOnt['RegulationType'].rchildren(),
                                  none=True)
    direction = obj_model.ontology.OntologyAttribute(kbOnt,
                                  terms = kbOnt['DirectionType'].rchildren(),
                                  default = kbOnt['forward'],
                                  none=False)

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'gene', 'regulatory_element', 'binding_factor', 'type',
                            'direction', 'identifiers', 'references', 'comments')


class PtmSite(core.PolymerLocus):
    """ Knowledge of protein modification sites

    Attributes:
        id (:obj:`str`): identifier
        name (:obj:`str`): name
        modified_protein (:obj:`ProteinSpeciesType`): modified protein
        type (:obj:`str`): type of modification (phosphorylation, methylation, etc...)
        modified_residue (:obj:`str`): residue name and position in protein sequence
        abundance_ratio (:obj:`int`): ratio of modified protein abundance
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references
        identifiers (:obj:`list` of :obj:`Identifier`): identifiers
    """
    type = obj_model.StringAttribute()
    modified_protein = obj_model.ManyToOneAttribute('ProteinSpeciesType', related_name='ptm_sites')
    modified_residue = obj_model.StringAttribute()
    abundance_ratio = obj_model.FloatAttribute()

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'modified_protein', 'type', 'modified_residue',
                           'abundance_ratio', 'identifiers', 'references', 'comments')


class GenericLocus(obj_model.Model):
    """ Start and end coordinates of exons and CDSs

    Attributes:
        start (:obj:`int`): start coordinate
        end (:obj:`int`): end coordinate

    Related attributes:
        transcripts (:obj:`list` of :obj:`TranscriptSpeciesType`): transcripts
        proteins (:obj:`list` of :obj:`ProteinSpeciesType`): proteins
    """
    start = obj_model.PositiveIntegerAttribute()
    end = obj_model.PositiveIntegerAttribute()

    class Meta(obj_model.Model.Meta):
        attribute_order = ('start', 'end')
        tabular_orientation = obj_model.TabularOrientation.cell
        ordering = ('start', 'end')

    @staticmethod
    def _serialize(start, end):
        """ Generate string representation

        Args:
            start (:obj:`int`): start coordinate
            end (:obj:`int`): end coordinate

        Returns:
            :obj:`str`: value of primary attribute
        """
        return '{}:{}'.format(start, end)

    def serialize(self):
        """ Generate string representation

        Returns:
            :obj:`str`: value of primary attribute
        """
        return self._serialize(self.start, self.end)


#####################
#####################
# Species types


class TranscriptSpeciesType(core.PolymerSpeciesType):
    """ Knowledge of a transcript (spliced RNA) species

    Attributes:
        gene (:obj:`GeneLocus`): gene
        exons (:obj:`list` of :obj:`LocusAttribute`): exon coordinates

    Related attributes:
        protein (:obj:`ProteinSpeciesType`): protein
    """
    gene = obj_model.ManyToOneAttribute(GeneLocus, related_name='transcripts')
    exons = LocusAttribute(related_name='transcripts')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'gene', 'exons', 'identifiers', 'references', 'comments')

    def get_seq(self):
        """ Get the 5' to 3' sequence

        Returns:
            :obj:`Bio.Seq.Seq`: sequence
        """
        ordered_exons = sorted(self.exons, key=lambda x: x.start)

        dna_seq = self.gene.polymer.get_subseq(
                    start=ordered_exons[0].start, end=ordered_exons[-1].end)

        adjusted_exons = [(i.start - ordered_exons[0].start, i.end - ordered_exons[0].start + 1) \
            for i in ordered_exons]

        spliced_dna_seq = Bio.Seq.Seq('', alphabet=Bio.Alphabet.DNAAlphabet())
        for exon in adjusted_exons:
            spliced_dna_seq += dna_seq[exon[0]:exon[1]]

        if are_terms_equivalent(self.gene.strand, kbOnt['negative']):
            spliced_dna_seq = spliced_dna_seq.reverse_complement()

        return spliced_dna_seq.transcribe()

    def get_empirical_formula(self):
        """ Get the empirical formula for a transcript (spliced RNA) species with

        * 5' monophosphate
        * Deprotonated phosphate oxygens

        :math:`N_A * AMP + N_C * CMP + N_G * GMP + N_U * UMP - (L-1) * OH`

        Returns:
           :obj:`chem.EmpiricalFormula`: empirical formula
        """
        seq = self.get_seq()
        n_a = seq.count('A')
        n_c = seq.count('C')
        n_g = seq.count('G')
        n_u = seq.count('U')
        l = len(seq)

        formula = chem.EmpiricalFormula()
        formula.C = 10 * n_a + 9 * n_c + 10 * n_g + 9 * n_u
        formula.H = 12 * n_a + 12 * n_c + 12 * n_g + 11 * n_u - (l - 1)
        formula.N = 5 * n_a + 3 * n_c + 5 * n_g + 2 * n_u
        formula.O = 7 * n_a + 8 * n_c + 8 * n_g + 9 * n_u - (l - 1)
        formula.P = n_a + n_c + n_g + n_u

        return formula

    def get_charge(self):
        """ Get the charge for a transcript (spliced RNA) species with

        * 5' monophosphate
        * Deprotonated phosphate oxygens

        :math:`-L - 1`

        Returns:
           :obj:`int`: charge
        """
        return -self.get_len() - 1

    def get_mol_wt(self):
        """ Get the molecular weight for a transcript (spliced RNA) species with

        * 5' monophosphate
        * Deprotonated phosphate oxygens

        Returns:
            :obj:`float`: molecular weight (Da)
        """
        return self.get_empirical_formula().get_molecular_weight()


class ProteinSpeciesType(core.PolymerSpeciesType):
    """ Knowledge of a protein monomer

    Attributes:
        uniprot (:obj:`str`): uniprot id
        transcript (:obj:`TranscriptSpeciesType`): transcript
        coding_regions (:obj:`list` of :obj:`LocusAttribute`): CDS coordinates

    Related attributes:
        regulatory_elements (:obj:`list` of `RegulatoryElementLocus`): potential binding sites
        regulatory_modules (:obj:`list` of `RegulatoryModule`): regulatory DNA binding sites
        ptm_sites (:obj:list` of `PtmSite`): protein modification sites
    """

    uniprot = obj_model.StringAttribute()
    transcript = obj_model.OneToOneAttribute(TranscriptSpeciesType, related_name='protein')
    coding_regions = LocusAttribute(related_name='proteins')

    class Meta(obj_model.Model.Meta):
        verbose_name = 'Protein'
        attribute_order = ('id', 'name', 'uniprot', 'transcript', 'coding_regions', 
                           'identifiers', 'references', 'comments')

    def get_seq(self, table=1, cds=True):
        """ Get the 5' to 3' sequence

        Args:
            table (:obj:`int`, optional): NCBI identifier for translation table
                                        (default = standard table)
            cds (:obj:`bool`, optional): True indicates the sequence is a complete CDS

        Returns:
            :obj:`Bio.Seq.Seq`: sequence
        """
        ordered_cds = sorted(self.coding_regions, key=lambda x: x.start)

        dna_seq = self.transcript.gene.polymer.get_subseq(
                    start=ordered_cds[0].start, end=ordered_cds[-1].end)

        adjusted_cds = [(i.start - ordered_cds[0].start, i.end - ordered_cds[0].start + 1) \
            for i in ordered_cds]

        spliced_dna_seq = Bio.Seq.Seq('', alphabet=Bio.Alphabet.DNAAlphabet())
        for cds in adjusted_cds:
            spliced_dna_seq += dna_seq[cds[0]:cds[1]]

        if are_terms_equivalent(self.transcript.gene.strand, kbOnt['negative']):
            spliced_dna_seq = spliced_dna_seq.reverse_complement()

        return spliced_dna_seq.transcribe().translate(table=table, cds=cds)

    def get_empirical_formula(self, table=1, cds=True):
        """ Get the empirical formula

        Args:
            table (:obj:`int`, optional): NCBI identifier for translation table
                                        (default = standard table)
            cds (:obj:`bool`, optional): True indicates the sequence is a complete CDS

        Returns:
            :obj:`chem.EmpiricalFormula`: empirical formula
        """
        seq = self.get_seq(table=table, cds=cds)
        l = len(seq)

        n_a = seq.count('A')  # Ala: Alanine (C3 H7 N O2)
        n_r = seq.count('R')  # Arg: Arginine (C6 H14 N4 O2)
        n_n = seq.count('N')  # Asn: Asparagine (C4 H8 N2 O3)
        n_d = seq.count('D')  # Asp: Aspartic acid (C4 H7 N O4)
        n_c = seq.count('C')  # Cys: Cysteine (C3 H7 N O2 S)

        n_q = seq.count('Q')  # Gln: Glutamine (C5 H10 N2 O3)
        n_e = seq.count('E')  # Glu: Glutamic acid (C5 H9 N O4)
        n_g = seq.count('G')  # Gly: Glycine (C2 H5 N O2)
        n_h = seq.count('H')  # His: Histidine (C6 H9 N3 O2)
        n_i = seq.count('I')  # Ile: Isoleucine (C6 H13 N O2)

        n_l = seq.count('L')  # Leu: Leucine (C6 H13 N O2)
        n_k = seq.count('K')  # Lys: Lysine (C6 H14 N2 O2)
        n_m = seq.count('M')  # Met: Methionine (C5 H11 N O2 S)
        n_f = seq.count('F')  # Phe: Phenylalanine (C9 H11 N O2)
        n_p = seq.count('P')  # Pro: Proline (C5 H9 N O2)

        n_s = seq.count('S')  # Ser: Serine (C3 H7 N O3)
        n_t = seq.count('T')  # Thr: Threonine (C4 H9 N O3)
        n_w = seq.count('W')  # Trp: Tryptophan (C11 H12 N2 O2)
        n_y = seq.count('Y')  # Tyr: Tyrosine (C9 H11 N O3)
        n_v = seq.count('V')  # Val: Valine (C5 H11 N O2)

        formula = chem.EmpiricalFormula()

        formula.C = 3 * n_a + 6 * n_r + 4 * n_n + 4 * n_d + 3 * n_c + \
            5 * n_q + 5 * n_e + 2 * n_g + 6 * n_h + 6 * n_i + \
            6 * n_l + 6 * n_k + 5 * n_m + 9 * n_f + 5 * n_p + \
            3 * n_s + 4 * n_t + 11 * n_w + 9 * n_y + 5 * n_v

        formula.H = 7 * n_a + 14 * n_r + 8 * n_n + 7 * n_d + 7 * n_c + \
            10 * n_q + 9 * n_e + 5 * n_g + 9 * n_h + 13 * n_i + \
            13 * n_l + 14 * n_k + 11 * n_m + 11 * n_f + 9 * n_p + \
            7 * n_s + 9 * n_t + 12 * n_w + 11 * n_y + 11 * n_v - 2 * (l - 1)

        formula.N = 1 * n_a + 4 * n_r + 2 * n_n + 1 * n_d + 1 * n_c + \
            2 * n_q + 1 * n_e + 1 * n_g + 3 * n_h + 1 * n_i + \
            1 * n_l + 2 * n_k + 1 * n_m + 1 * n_f + 1 * n_p + \
            1 * n_s + 1 * n_t + 2 * n_w + 1 * n_y + 1 * n_v

        formula.O = 2 * n_a + 2 * n_r + 3 * n_n + 4 * n_d + 2 * n_c + \
            3 * n_q + 4 * n_e + 2 * n_g + 2 * n_h + 2 * n_i + \
            2 * n_l + 2 * n_k + 2 * n_m + 2 * n_f + 2 * n_p + \
            3 * n_s + 3 * n_t + 2 * n_w + 3 * n_y + 2 * n_v - (l - 1)

        formula.S = n_c + n_m
        return formula

    def get_charge(self, table=1, cds=True):
        """ Get the charge at physiological pH

        Args:
            table (:obj:`int`, optional): NCBI identifier for translation table
                                        (default = standard table)
            cds (:obj:`bool`, optional): True indicates the sequence is a complete CDS

        Returns:
            :obj:`int`: charge
        """
        seq = self.get_seq(table=table, cds=cds)

        n_r = seq.count('R')
        n_h = seq.count('H')
        n_k = seq.count('K')
        n_d = seq.count('D')
        n_e = seq.count('E')

        return (n_r + n_h + n_k) - (n_d + n_e)

    def get_mol_wt(self, table=1, cds=True):
        """ Get the molecular weight

        Args:
            table (:obj:`int`, optional): NCBI identifier for translation table
                                        (default = standard table)
            cds (:obj:`bool`, optional): True indicates the sequence is a complete CDS

        Returns:
            :obj:`float`: molecular weight
        """
        return self.get_empirical_formula(table=table, cds=cds).get_molecular_weight()
