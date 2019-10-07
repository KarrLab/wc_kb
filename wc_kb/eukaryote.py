""" Schema to represent a knowledge base to build models of eukaryotes

:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2018-09-10
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_kb import core
from wc_utils.util import chem
import Bio.Alphabet
import Bio.Seq
import enum
import obj_tables
import re

#####################
#####################
# Enumeration classes


class ActivityLevel(enum.Enum):
    """ Activity level of regulatory element"""
    active = 1
    poised = 2
    repressed = 3
    inactive = 4
    na = 5

class RegulationType(enum.Enum):
    """ Type of regulation between a regulatory element and a gene """
    proximal = 1
    distal = 2

class RegulatoryDirection(int, enum.Enum):
    """ The direction of regulation """
    activation = 1
    repression = -1
    unknown = 0

class TranscriptType(enum.Enum):
    """ Type of transcript """
    mRna = 1
    rRna = 2
    tRna = 3
    itRna = 4    

#####################
#####################
# Attributes

class LocusAttribute(obj_tables.ManyToManyAttribute):
    """ Start and end coordinates attribute """

    def __init__(self, related_name='', verbose_name='', verbose_related_name='', description=''):
        """
        Args:
            related_name (:obj:`str`, optional): name of related attribute on `related_class`
            verbose_name (:obj:`str`, optional): verbose name
            verbose_related_name (:obj:`str`, optional): verbose related name
            description (:obj:`str`, optional): description
        """
        super(LocusAttribute, self).__init__(GenericLocus,
                                            related_name=related_name, min_related=0, min_related_rev=0,
                                            verbose_name=verbose_name, verbose_related_name=verbose_related_name, description=description)

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

        return ', '.join(obj_tables.serialize() for obj_tables in coordinates)

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

        obj_pattern = r'([0-9]+) *\: *([0-9]+)'
        lst_pattern = r'^{}( *, *{})*$'.format(obj_pattern, obj_pattern)
        if not re.match(lst_pattern, value, flags=re.I):
            return (None, obj_tables.InvalidAttribute(self, ['Incorrectly formatted list of coordinates: {}'.format(value)]))

        objs = []
        for pat_match in re.findall(obj_pattern, value, flags=re.I):
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


class RegDirectionAttribute(obj_tables.ManyToManyAttribute):
    """ Regulatory direction attribute """

    def __init__(self, related_name='', verbose_name='', verbose_related_name='', description=''):
        """
        Args:
            related_name (:obj:`str`, optional): name of related attribute on `related_class`
            verbose_name (:obj:`str`, optional): verbose name
            verbose_related_name (:obj:`str`, optional): verbose related name
            description (:obj:`str`, optional): description
        """
        super(RegDirectionAttribute, self).__init__('TranscriptionFactorRegulation',
                                            related_name=related_name, 
                                            verbose_name=verbose_name, 
                                            verbose_related_name=verbose_related_name, description=description)

    def serialize(self, directions, encoded=None):
        """ Serialize related object
        Args:
            directions (:obj:`list` of :obj:`Model`): a list of instances of TFdirection Python representation
            encoded (:obj:`dict`, optional): dictionary of objects that have already been encoded
        Returns:
            :obj:`str`: simple Python representation
        """
        if not directions:
            return ''

        return ', '.join(obj_tables.serialize() for obj_tables in directions)

    def deserialize(self, value, objects, decoded=None):
        """ Deserialize value
        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model
            decoded (:obj:`dict`, optional): dictionary of objects that have already been decoded
        Returns:
            :obj:`tuple` of :obj:`list` of :obj:`RegDirection`, :obj:`InvalidAttribute` or :obj:`None`: :obj:`tuple` of cleaned value
                and cleaning error
        """
        if not value:
            return ([], None)

        tf_id = ProteinSpeciesType.id.pattern[1:-1]
        direction = '|'.join(rd.name for rd in RegulatoryDirection)
        obj_pattern = r'({}) *\: *({})'.format(tf_id, direction)
        lst_pattern = r'^{}( *, *{})*$'.format(obj_pattern, obj_pattern)
        if not re.match(lst_pattern, value, flags=re.I):
            return (None, obj_tables.InvalidAttribute(self, ['Incorrectly formatted list of transcription factor regulation: {}'.format(value)]))

        objs = []
        errors = []
        for pat_match in re.findall(obj_pattern, value, flags=re.I):
            tf = pat_match[0]
            direction = pat_match[-1]
            
            if self.related_class not in objects:
                objects[self.related_class] = {}           
            
            serialized_value = self.related_class._serialize(tf, direction)           
            
            if serialized_value in objects[self.related_class]:
                obj = objects[self.related_class][serialized_value]
            else:
                obj, error = self.related_class.deserialize(serialized_value, objects)
                if error:
                    errors.append(error)
                else:    
                    objects[self.related_class][serialized_value] = obj
                    objs.append(obj)
        
        if errors:
            return (None, obj_tables.InvalidAttribute(self, errors))

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
    symbol = obj_tables.StringAttribute()
    homologs = obj_tables.LongStringAttribute()

    class Meta(obj_tables.Model.Meta):
        verbose_name = 'Gene'
        verbose_name_plural = 'Genes'
        attribute_order = ('id', 'name', 'synonyms', 'symbol', 'homologs', 'polymer', 'strand', 'start',
                           'end', 'identifiers', 'references', 'comments')


class TranscriptionFactorRegulation(obj_tables.Model):
    """ Transcription factor and the direction of transcriptional regulation

    Attributes:
        transcription_factor (:obj:`ProteinSpeciesType`): transcription factor
        direction (:obj:`RegulatoryDirection`): regulatory direction

    Related attributes:
        regulatory_modules (:obj:`list` of `RegulatoryModule`): regulatory modules
    """

    transcription_factor = obj_tables.ManyToOneAttribute('ProteinSpeciesType', related_name='transcription_factor_regulation')
    direction = obj_tables.EnumAttribute(RegulatoryDirection)    

    class Meta(obj_tables.Model.Meta):
        attribute_order = ('transcription_factor', 'direction')
        frozen_columns = 1
        table_format = obj_tables.TableFormat.cell
        ordering = ('transcription_factor', 'direction')

    @staticmethod
    def _serialize(transcription_factor_id, direction_name):
        """ Generate string representation

        Args:
            transcription_factor_id (:obj:`str`): transcription factor id
            direction_name (:obj:`str`): regulatory direction name

        Returns:
            :obj:`str`: value of primary attribute
        """
        return '{}:{}'.format(transcription_factor_id, direction_name)

    def serialize(self):
        """ Generate string representation

        Returns:
            :obj:`str`: value of primary attribute
        """
        return self._serialize(self.transcription_factor.id, self.direction.name)

    @classmethod
    def deserialize(cls, value, objects):
        """ Deserialize value
        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model
        Returns:
            :obj:`tuple` of `list` of `TranscriptionFactorRegulation`, `InvalidAttribute` or `None`: tuple of cleaned value
                and cleaning error
        """
        if cls in objects and value in objects[cls]:
            return (objects[cls][value], None)

        tf_id = ProteinSpeciesType.id.pattern[1:-1]
        direction = '|'.join(rd.name for rd in RegulatoryDirection)
        pattern = r'^({}) *\: *({})$'.format(tf_id, direction)        
        match = re.match(pattern, value, flags=re.I)       

        if match:
            errors = []

            tf_reg_str = match.group(1)       
            if ProteinSpeciesType in objects and tf_reg_str in objects[ProteinSpeciesType]:
                transcription_factor = objects[ProteinSpeciesType][tf_reg_str]
            else:
                errors.append('Undefined transcription factor "{}"'.format(tf_reg_str))

            direction_str = match.group(match.lastindex)
            if direction_str in [i.name for i in list(RegulatoryDirection)]:
                direction = [i for i in list(RegulatoryDirection) if i.name==direction_str][0]
            else:
                errors.append('Undefined regulatory direction "{}"'.format(direction_str))

            if errors:
                return (None, obj_tables.InvalidAttribute(cls, errors)) 
            else:
                obj = cls(transcription_factor=transcription_factor, direction=direction)
                if cls not in objects:
                    objects[cls] = {}
                objects[cls][obj.serialize()] = obj
                return (obj, None)

        return (None, obj_tables.InvalidAttribute(cls, ['Invalid transcription factor regulation']))
            

class RegulatoryModule(obj_tables.Model):
    """ Knowledge about regulatory modules

    Attributes:
        id (:obj:`str`): identifier
        name (:obj:`str`): name
        gene (:obj:`GeneLocus`): gene
        promoter (:obj:`str`): promoter ensembl ID
        activity (:obj:`ActivityLevel`): cell-type specific activity level
        type (:obj:`RegulationType`): type of regulation (proximal or distal)
        transcription_factor_regulation (:obj:`TranscriptionFactorRegulation`): 
            transcription factor and direction of regulation
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references
        identifiers (:obj:`list` of :obj:`Identifier`): identifiers
    """
    id = obj_tables.SlugAttribute(primary=True, unique=True)
    name = obj_tables.StringAttribute()
    gene = obj_tables.ManyToOneAttribute(GeneLocus, related_name='regulatory_modules')    
    promoter = obj_tables.StringAttribute()
    activity = obj_tables.EnumAttribute(ActivityLevel)
    type = obj_tables.EnumAttribute(RegulationType)
    transcription_factor_regulation = RegDirectionAttribute(related_name='regulatory_modules')
    comments = obj_tables.LongStringAttribute()
    references = obj_tables.ManyToManyAttribute(core.Reference, related_name='regulatory_modules')
    identifiers = core.IdentifierAttribute(related_name='regulatory_modules')

    class Meta(obj_tables.Model.Meta):
        attribute_order = ('id', 'name', 'gene', 'promoter', 'activity', 'type',
                            'transcription_factor_regulation', 'identifiers', 'references', 'comments')


class PtmSite(core.PolymerLocus):
    """ Knowledge of protein modification sites

    Attributes:
        modified_protein (:obj:`ProteinSpeciesType`): modified protein
        type (:obj:`str`): type of modification (phosphorylation, methylation, etc...)
        modified_residue (:obj:`str`): residue name and position in protein sequence
        fractional_abundance (:obj:`int`): ratio of modified protein abundance
    """
    type = obj_tables.StringAttribute()
    modified_protein = obj_tables.ManyToOneAttribute('ProteinSpeciesType', related_name='ptm_sites')
    modified_residue = obj_tables.StringAttribute()
    fractional_abundance = obj_tables.FloatAttribute()

    class Meta(obj_tables.Model.Meta):
        attribute_order = ('id', 'name', 'modified_protein', 'type', 'modified_residue',
                           'fractional_abundance', 'identifiers', 'references', 'comments')


class GenericLocus(obj_tables.Model):
    """ Start and end coordinates of exons and CDSs

    Attributes:
        start (:obj:`int`): start coordinate
        end (:obj:`int`): end coordinate

    Related attributes:
        transcripts (:obj:`list` of :obj:`TranscriptSpeciesType`): transcripts
        proteins (:obj:`list` of :obj:`ProteinSpeciesType`): proteins
    """
    start = obj_tables.PositiveIntegerAttribute()
    end = obj_tables.PositiveIntegerAttribute()

    class Meta(obj_tables.Model.Meta):
        attribute_order = ('start', 'end')
        table_format = obj_tables.TableFormat.cell
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
        type (:obj:`TranscriptType`): type

    Related attributes:
        protein (:obj:`ProteinSpeciesType`): protein
    """
    gene = obj_tables.ManyToOneAttribute(GeneLocus, related_name='transcripts')    
    exons = LocusAttribute(related_name='transcripts')
    type = obj_tables.EnumAttribute(TranscriptType)

    class Meta(obj_tables.Model.Meta):
        verbose_name = 'Transcript'
        verbose_name_plural = 'Transcripts'
        attribute_order = ('id', 'name', 'gene', 'exons', 'type', 'identifiers', 'references', 'comments')

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

        if self.gene.strand == core.PolymerStrand.negative:
            spliced_dna_seq = spliced_dna_seq.reverse_complement()

        return spliced_dna_seq.transcribe()

    def get_empirical_formula(self, seq_input=None):
        """ Get the empirical formula for a transcript (spliced RNA) species with

        * 5' monophosphate
        * Deprotonated phosphate oxygens

        :math:`N_A * AMP + N_C * CMP + N_G * GMP + N_U * UMP - (L-1) * OH`

        Args:
            seq_input (:obj:`Bio.Seq.Seq`, optional): if provided, the method will use it
                instead of reading from fasta file to reduce IO operation 

        Returns:
           :obj:`chem.EmpiricalFormula`: empirical formula
        """
        if seq_input:
            seq = seq_input
        else:    
            seq = self.get_seq()
        
        n_a = seq.upper().count('A')
        n_c = seq.upper().count('C')
        n_g = seq.upper().count('G')
        n_u = seq.upper().count('U')
        l = len(seq)

        formula = chem.EmpiricalFormula()
        formula.C = 10 * n_a + 9 * n_c + 10 * n_g + 9 * n_u
        formula.H = 12 * n_a + 12 * n_c + 12 * n_g + 11 * n_u - (l - 1)
        formula.N = 5 * n_a + 3 * n_c + 5 * n_g + 2 * n_u
        formula.O = 7 * n_a + 8 * n_c + 8 * n_g + 9 * n_u - (l - 1)
        formula.P = n_a + n_c + n_g + n_u

        return formula

    def get_charge(self, seq_input=None):
        """ Get the charge for a transcript (spliced RNA) species with

        * 5' monophosphate
        * Deprotonated phosphate oxygens

        :math:`-L - 1`

        Args:
            seq_input (:obj:`Bio.Seq.Seq`, optional): if provided, the method will use it
                instead of reading from fasta file to reduce IO operation

        Returns:
           :obj:`int`: charge
        """
        if seq_input:
            length = len(seq_input)
        else:
            length = len(self.get_seq())    
        return -length - 1

    def get_mol_wt(self, seq_input=None):
        """ Get the molecular weight for a transcript (spliced RNA) species with

        * 5' monophosphate
        * Deprotonated phosphate oxygens

        Args:
            seq_input (:obj:`Bio.Seq.Seq`, optional): if provided, the method will use it
                instead of reading from fasta file to reduce IO operation        

        Returns:
            :obj:`float`: molecular weight (Da)
        """
        if seq_input:
            return self.get_empirical_formula(seq_input=seq_input).get_molecular_weight()
        else:
            return self.get_empirical_formula().get_molecular_weight()


class ProteinSpeciesType(core.PolymerSpeciesType):
    """ Knowledge of a protein monomer

    Attributes:
        uniprot (:obj:`str`): uniprot id
        transcript (:obj:`TranscriptSpeciesType`): transcript
        coding_regions (:obj:`list` of :obj:`LocusAttribute`): CDS coordinates

    Related attributes:
        transcription_factor_regulation (:obj:`list` of `TranscriptionFactorRegulation`): transcription factor regulation
        ptm_sites (:obj:list` of `PtmSite`): protein modification sites
    """

    uniprot = obj_tables.StringAttribute()
    transcript = obj_tables.OneToOneAttribute(TranscriptSpeciesType, related_name='protein')
    coding_regions = LocusAttribute(related_name='proteins')

    class Meta(obj_tables.Model.Meta):
        verbose_name = 'Protein'
        verbose_name_plural = 'Proteins'
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
        for i in adjusted_cds:
            spliced_dna_seq += dna_seq[i[0]:i[1]]

        if self.transcript.gene.strand == core.PolymerStrand.negative:
            spliced_dna_seq = spliced_dna_seq.reverse_complement()

        return spliced_dna_seq.transcribe().translate(table=table, cds=cds)

    def get_empirical_formula(self, table=1, cds=True, seq_input=None):
        """ Get the empirical formula

        Args:
            table (:obj:`int`, optional): NCBI identifier for translation table
                                        (default = standard table)
            cds (:obj:`bool`, optional): True indicates the sequence is a complete CDS
            seq_input (:obj:`Bio.Seq.Seq`, optional): if provided, the method will use it
                instead of reading from fasta file to reduce IO operation

        Returns:
            :obj:`chem.EmpiricalFormula`: empirical formula
        """
        if seq_input:
            seq = seq_input
        else:    
            seq = self.get_seq(table=table, cds=cds)
        l = len(seq) - seq.count('*')    

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

    def get_charge(self, table=1, cds=True, seq_input=None):
        """ Get the charge at physiological pH

        Args:
            table (:obj:`int`, optional): NCBI identifier for translation table
                                        (default = standard table)
            cds (:obj:`bool`, optional): True indicates the sequence is a complete CDS
            seq_input (:obj:`Bio.Seq.Seq`, optional): if provided, the method will use it
                instead of reading from fasta file to reduce IO operation

        Returns:
            :obj:`int`: charge
        """
        if seq_input:
            seq = seq_input
        else:    
            seq = self.get_seq(table=table, cds=cds)

        n_r = seq.count('R')
        n_h = seq.count('H')
        n_k = seq.count('K')
        n_d = seq.count('D')
        n_e = seq.count('E')

        return (n_r + n_h + n_k) - (n_d + n_e)

    def get_mol_wt(self, table=1, cds=True, seq_input=None):
        """ Get the molecular weight

        Args:
            table (:obj:`int`, optional): NCBI identifier for translation table
                                        (default = standard table)
            cds (:obj:`bool`, optional): True indicates the sequence is a complete CDS
            seq_input (:obj:`Bio.Seq.Seq`, optional): if provided, the method will use it
                instead of reading from fasta file to reduce IO operation

        Returns:
            :obj:`float`: molecular weight
        """
        if seq_input:
            return self.get_empirical_formula(table=table, cds=cds, seq_input=seq_input).get_molecular_weight()
        else:    
            return self.get_empirical_formula(table=table, cds=cds).get_molecular_weight()
