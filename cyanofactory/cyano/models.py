'''
Cyanofactory data model. Contains the database layout.

Based on Wholecell KB.
'''

from __future__ import unicode_literals
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from django.contrib.auth.models import User
from django.contrib.auth.models import Group
from django.core import validators
from django.core.exceptions import ValidationError, ObjectDoesNotExist
from django.core.urlresolvers import reverse
from django.db import models
from django.db.models import F, Model, OneToOneField, CharField, IntegerField, URLField, PositiveIntegerField, FloatField, BooleanField, SlugField, TextField, DateTimeField, options, permalink, SET_NULL, Min, Max
from django.db.models.query import EmptyQuerySet
from django.utils.http import urlencode
from itertools import chain
from cyano.templatetags.templatetags import set_time_zone
import math
import re
import settings
import subprocess
from datetime import datetime
from cyano.exceptions import WidAlreadyExist
from cyano.cache import Cache
import StringIO
from django.template import loader, Context
from cyano.history import HistoryForeignKey as ForeignKey
from cyano.history import HistoryManyToManyField as ManyToManyField
from django.db.models.loading import get_model
from django.dispatch.dispatcher import receiver
from django.db.models.signals import m2m_changed

def enum(**enums):
    return type(str('Enum'), (), enums)

PermissionEnum = enum(
    FULL_ACCESS = "Full",
    READ_NORMAL = "Read Normal",
    READ_DELETE = "Read Delete",
    READ_PERMISSION = "Read Permission",
    READ_HISTORY = "Read History",
    WRITE_NORMAL = "Write Normal",
    WRITE_DELETE = "Write Delete",
    WRITE_PERMISSION = "Write Permission"
)

''' BEGIN: choices '''

CHOICES_DBOPERATION = (
    ('I', 'Insert'),
    ('U', 'Update'),
    ('D', 'Delete'),
    ('X', 'Unknown')
)

CHOICES_DIRECTION = (
    ('f', 'Forward'),
    ('r', 'Reverse'),
)

CHOICES_REFERENCE_TYPE = (
    ('article', 'Article'),
    ('book', 'Book'),
    ('thesis', 'Thesis'),
    ('misc', 'Miscellaneous'),
)

CHOICES_STRANDEDNESS = (
    ('dsDNA', 'dsDNA'),
    ('ssDNA', 'ssDNA'),
    ('xsDNA', 'xsDNA'),
)

CHOICES_VMAX_UNITS = (
    ('1/min', '1/min'),
    ('U/mg', 'U/mg'),
)

CHOICES_HOMOLOG_SPECIES = (
    ('B. subtilis', 'B. subtilis'),
    ('E. coli', 'E. coli'),
    ('M. hyopneumoniae', 'M. hyopneumoniae'),
    ('M. mobile', 'M. mobile'),
    ('M. pneumoniae', 'M. pneumoniae'),
    ('S. coelicolor', 'S. coelicolor'),
    ('S. oneidensis', 'S. oneidensis'),
)

HOMOLOG_SPECIES_URLS = {
    "B. subtilis": "http://www.genome.jp/dbget-bin/www_bget?bsu:%s",
    "E. coli": "http://www.genome.jp/dbget-bin/www_bget?eco:%s",
    "M. hyopneumoniae": "http://www.genome.jp/dbget-bin/www_bget?mhj:%s",
    "M. mobile": "http://www.genome.jp/dbget-bin/www_bget?mmo:%s",
    "M. pneumoniae": "http://www.genome.jp/dbget-bin/www_bget?mpn:%s",
    "S. coelicolor": "http://www.genome.jp/dbget-bin/www_bget?sco:%s",
    "S. oneidensis": "http://www.genome.jp/dbget-bin/www_bget?son:%s",
}

CHOICES_CROSS_REFERENCE_SOURCES = (
    ('ATCC', 'ATCC'),
    ('BiGG', 'BiGG'),
    ('BioCyc', 'BioCyc'),
    ('BioProject', 'BioProject'), #http://www.ncbi.nlm.nih.gov/bioproject/%s
    ('CAS', 'CAS'),
    ('ChEBI', 'ChEBI'),
    ('CMR', 'CMR'),
    ('EC', 'EC'),
    ('GenBank', 'GenBank'),
    ('ISBN', 'ISBN'),
    ('KEGG', 'KEGG'),
    ('KNApSAcK', 'KNApSAcK'),
    ('LipidBank', 'LipidBank'),
    ('LIPIDMAPS', 'LIPIDMAPS'),
    ('PDB', 'PDB'),
    ('PDBCCD', 'PDBCCD'),
    ('PubChem', 'PubChem'),
    ('PubMed', 'PubMed'),
    ('RefSeq', 'RefSeq'),
    ('SABIO-RK', 'SABIO-RK'),
    ('SwissProt', 'SwissProt'),
    ('Taxonomy', 'Taxonomy'),
    ('ThreeDMET', 'ThreeDMET'),    
    ('URL', 'URL'),
)

CROSS_REFERENCE_SOURCE_URLS = {
    "ATCC": "http://www.atcc.org/ATCCAdvancedCatalogSearch/ProductDetails/tabid/452/Default.aspx?Template=bacteria&ATCCNum=%s",
    "BiGG": "http://bigg.ucsd.edu/bigg/postMet.pl?organism=3307911&organism=1461534&organism=222668&organism=3277493&organism=2795088&organism=2423527&organism=1823466&compartment_list=any&pathway_list=any&name_text=%s",
    "BioCyc": "http://biocyc.org/MGEN243273/NEW-IMAGE?type=GENE&object=%s",
    "BioProject": "http://www.ncbi.nlm.nih.gov/bioproject/%s",
    "CAS": "http://www.ncbi.nlm.nih.gov/sites/entrez?db=pccompound&term=%s",
    "ChEBI": "http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:%s",
    "CMR": "http://cmr.jcvi.org/tigr-scripts/CMR/shared/GenePage.cgi?locus=%s",
    "EC": "http://www.expasy.ch/enzyme/%s",
    "GenBank": "http://www.ncbi.nlm.nih.gov/sites/gquery?term=%s",
    "ISBN": "http://isbndb.com/search-all.html?kw=%s",
    "KEGG": "http://www.genome.jp/dbget-bin/www_bget?cpd:%s",
    "KNApSAcK": "http://kanaya.naist.jp/knapsack_jsp/information.jsp?word=%s",
    "LipidBank": "http://lipidbank.jp/cgi-bin/detail.cgi?id=%s",
    "LIPIDMAPS": "http://www.lipidmaps.org/data/get_lm_lipids_dbgif.php?LM_ID=%s",
    "PDB": "http://www.pdb.org/pdb/explore/explore.do?structureId=%s",
    "PDBCCD": "http://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/%s",
    "PubChem": "http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?viewopt=PubChem&sid=%s",
    "PubMed": "http://www.ncbi.nlm.nih.gov/pubmed/%s",
    "RefSeq": "http://www.ncbi.nlm.nih.gov/nuccore/%s",
    "SABIO-RK": "http://sabio.villa-bosch.de/kineticLawEntry.jsp?kinlawid=%s&viewData=true",
    "SwissProt": "http://www.uniprot.org/uniprot/%s",
    "Taxonomy": "http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=%s",
    "ThreeDMET": "http://www.3dmet.dna.affrc.go.jp/bin2/show_data.e?acc=%s",    
    "URL": "%s",
}

CHOICES_GENETIC_CODE = (
    ('1', 'Standard'),
    ('2', 'Vertebrate'),
    ('3', 'Yeast'),
    ('4', 'Mold, protozoa, coelenterate mitochondria, mycoplasma, and spiroplasma'),
    ('5', 'Invertebrate mitochondria'),
    ('6', 'Ciliate, dasycladacean and hexamita'),
    ('9', 'Echinoderm and flatworm mitochondria'),
    ('10', 'Euplotid'),
    ('11', 'Bacteria, archaea and plant plastids'),
    ('12', 'Alternative yeast'),
    ('13', 'Ascidian mitochondria'),
    ('14', 'Alternative flatworm mitochondria'),
    ('15', 'Blepharisma'),
    ('16', 'Chlorophycean mitochondria'),
    ('21', 'Trematode mitochondria'),
    ('22', 'Scenedesmus obliquus mitochondria'),
    ('23', 'Thraustochytrium mitochondria'),
    ('24', 'Pterobranchia mitochondria'),
)

CHOICES_REACTION_DIRECTION = (
    ('f', 'Forward'),
    ('b', 'Backward'),
    ('r', 'Reversible'),
)

CHOICES_SIGNAL_SEQUENCE_LOCATION = (
    ('N', 'N-terminus'),
    ('C', 'C-terminus'),
)

CHOICES_SIGNAL_SEQUENCE_TYPE = (
    ('lipoprotein', 'Lipoprotein'),
    ('secretory', 'Secretory'),
)

''' END: CHOICES '''

# add model options
options.DEFAULT_NAMES = options.DEFAULT_NAMES + ('listing', 'concrete_entry_model', 'fieldsets', 'field_list', 'facet_fields', 'clean', 'validate_unique')

''' BEGIN: validators '''
def validate_dna_sequence(seq):
    validators.RegexValidator(regex=r'^[ACGT]+$', message='Enter a valid DNA sequence consisting of only the letters A, C, G, and T')(seq)
    
def validate_kinetics(reaction, direction):
    if direction == 'f':
        prop_name = 'kinetics_forward'        
    else:
        prop_name = 'kinetics_backward'
    kinetics = reaction[prop_name]
    
    wids = []
    for s in reaction['stoichiometry']:
        if (direction == 'f' and s['coefficient'] < 0) or (direction == 'r' and s['coefficient'] > 0):
            wids.append(s['molecule'])

    #rate law
    usedKmMax = 0
    usedVmax = 0
    for match in re.finditer(r'([a-z][a-z0-9_]*)', kinetics['rate_law'], flags=re.I):
        if match.group(1) == 'Vmax':
            usedVmax = 1
        elif match.group(1)[0:2] == 'Km':
            if len(match.group(1)) == 2:
                usedKmMax = max(usedKmMax, 1);
            elif match.group(1)[2:].isnumeric():
                usedKmMax = max(usedKmMax, int(float(match.group(1)[2:])))
            else:
                raise ValidationError({prop_name: 'Invalid rate law'})
        elif match.group(1) not in wids:
            raise ValidationError({prop_name: 'Invalid rate law, unknown molecule: "%s"' % (match.group(1), )})

    #Vmax
    if usedVmax and kinetics['vmax'] is None:
        raise ValidationError({prop_name: 'Invalid rate law'})
    if usedVmax and kinetics['vmax_unit'] not in ["1/min", "U/mg"]:
        raise ValidationError({prop_name: 'Invalid rate law'})

    #Km
    if not ((kinetics['rate_law'] is None or kinetics['rate_law'] == '') or (kinetics['km'] == '' and usedKmMax == 0) or (kinetics['km'] != '' and usedKmMax == len(kinetics['km'].split(', ')))):
        raise ValidationError({prop_name: 'Invalid rate law'})
    
'''
Test cases:
tests = [
    "!(ATP)",
    "  (  ATP ) ",
    "CTP  (  ATP ) ",
    "CTP | (  ATP ) ",
    "CTP | (  ATP & !ATP & GTP ) ",
    "CTP | (  ATP & !ATP & GTP)( ) ",
    "CTP | (  ATP & !ATP & GTP ) | ()",
    "CTP | (  ATP & !ATP & GTP ) | (H)",
    "CTP | (  ATP & !ATP & GTP ) | (H) K",
    "!!CTP",
    "CTP | !!(  ATP & !ATP & GTP ) | (H)  & K",
    "CTP | !(  ATP & !ATP & GTP ) | (H)  & K",
    "CTP | !(  ATP & (!ATP | UTP) & GTP ) | (H)  & K",
    ]
for test in tests:
    try:
        models.parse_regulatory_rule(test)
    except:
        print test
'''
def parse_regulatory_rule(equation, all_obj_data, species_wid):
    from cyano.helpers import getModel, getEntry
    import settings
    
    pre = ''
    blocks = []
    posts = []
    pattern = '%s'
    begin = 0
    end = 0
    sense = 0
            
    equation = equation or ''
    equation = equation.replace(" ", "")
    
    match = re.match(r'^([!\-]{0,1})\((.+)\)$', equation)
    if match:
        pattern = pattern.replace('%s', match.group(1) + '(%s)')
        equation = match.group(2)
    
    match = re.match(r'^([!\-]{0,1})(.+)$', equation)
    if match:
        pattern = pattern.replace('%s', match.group(1) + '%s')
        equation = match.group(2)
    
    if equation == '':
        raise ValidationError({'regulatory_rule': 'Invalid regulatory rule'})

    #split into blocks    
    for i in range(len(equation)):
        if equation[i] == '(':
            sense += 1
        if equation[i] == ')':
            sense -= 1
        if sense < 0:
            raise ValidationError({'regulatory_rule': 'Invalid regulatory rule'})

        if sense == 0:
            if len(blocks) == 0:
                pre = equation[0:begin]
            if equation[i] in ["|", "&", "-", "+"]:                
                blocks.append(equation[begin:i])
                posts.append(equation[i])
                begin = i + 1
            elif equation[i:i+2] in [">=", "<=", "=="]:                
                blocks.append(equation[begin:i])
                posts.append(equation[i:i+1])
                begin = i + 2
            elif equation[i] in [">", "<"]:                
                blocks.append(equation[begin:i])
                posts.append(equation[i])
                begin = i + 1
    blocks.append(equation[begin:])
    posts.append('')
    
    #check parenthesis match
    if sense != 0:
        raise ValidationError({'regulatory_rule': 'Invalid regulatory rule'})

    #if no further substructure
    if len(blocks) == 1:    
        if '(' in blocks[0]:
            raise ValidationError({'regulatory_rule': 'Invalid regulatory rule'})
        elif '!' in blocks[0]:
            raise ValidationError({'regulatory_rule': 'Invalid regulatory rule'})
        elif blocks[0].isnumeric() or blocks[0] in ["true", "false"]:
            return pattern % (pre + ' '.join([block + ' ' + post for block, post in zip(blocks, posts)]))
            
        wid = blocks[0]
        molecule = None
        if all_obj_data is None:
            molecule = getEntry(species_wid=species_wid, wid=wid)
        elif all_obj_data.has_key(wid):
            molecule = all_obj_data[wid]
        
        if molecule is not None:
            if isinstance(molecule, Entry):
                molecule_model = molecule.__class__
            else:
                molecule_model = getModel(molecule['model_type'])
            if issubclass(molecule_model, Molecule):
                return pattern % (pre + ' '.join([block + ' ' + post for block, post in zip(blocks, posts)]))
            else:
                raise ValidationError({'regulatory_rule': 'Invalid regulatory rule, referenced invalid object %s' % wid})            
        try:
            obj = Molecule.objects.get(species__wid=species_wid, wid=wid)
            blocks[0] = '<a href="%s">%s</a>' % (obj.get_absolute_url(species_wid), wid)
            return pattern % (pre + ' '.join([block + ' ' + post for block, post in zip(blocks, posts)]))
        except:
            pass
        raise ValidationError({'regulatory_rule': 'Invalid object "%s" in regulatory rule' % wid})
    
    #recurse on blocks
    for i in range(len(blocks)):
        blocks[i] = parse_regulatory_rule(blocks[i], all_obj_data, species_wid)
        
    return pattern % (pre + ' '.join([block + ' ' + post for block, post in zip(blocks, posts)]))
    
''' END: validators '''

class Permission(Model):
    name = CharField(max_length=255, blank=False, default='', verbose_name = "Permission name")
    description = CharField(max_length=255, blank=False, default='', verbose_name = "Permission description")
    
    permission_mapping = None
    permission_types = ["FULL_ACCESS",
                        "READ_NORMAL",
                        "READ_DELETE",
                        "READ_PERMISSION",
                        "READ_HISTORY",
                        "WRITE_NORMAL",
                        "WRITE_DELETE",
                        "WRITE_PERMISSION"]
    
    @staticmethod
    def get_by_name(name):
        if Permission.permission_mapping == None:
            Permission.permission_mapping = {}
            permissions = Permission.objects.all()
            for i, t in enumerate(Permission.permission_types, 1):
                Permission.permission_mapping[t] = permissions.get(pk = i)
        
        if not name in Permission.permission_mapping:
            return None
        return Permission.permission_mapping[name]
    
    @staticmethod
    def get_by_pk(pk):
        if pk < 1 or pk > len(Permission.permission_types):
            return None
        
        return Permission.get_by_name(Permission.permission_types[pk - 1])
    
    def __unicode__(self):
        return self.name

class UserPermission(Model):
    entry = ForeignKey("Entry", related_name = "user_permissions")
    user = ForeignKey("UserProfile", related_name = 'permissions')
    allow = ManyToManyField(Permission, verbose_name = 'Allowed Permissions', related_name = 'user_permission_allow')
    deny = ManyToManyField(Permission, verbose_name = 'Denied Permissions', related_name = 'user_permission_deny')
    
    def __unicode__(self):
        perm_allow = [str(1) if Permission.get_by_pk(x) in self.allow.all() else str(0) for x in range(1, 9)]
        perm_deny = [str(1) if Permission.get_by_pk(x) in self.deny.all() else str(0) for x in range(1, 9)]
        
        return "[{}, {}]".format(
            "".join(perm_allow), "".join(perm_deny)
        )

class GroupPermission(Model):
    entry = ForeignKey("Entry", related_name = "group_permissions")
    group = ForeignKey("GroupProfile", related_name = 'permissions')
    allow = ManyToManyField(Permission, verbose_name = 'Allowed Permissions', related_name = 'group_permission_allow')
    deny = ManyToManyField(Permission, verbose_name = 'Denied Permissions', related_name = 'group_permission_deny')

    def __unicode__(self):
        perm_allow = [str(1) if Permission.get_by_pk(x) in self.allow.all() else str(0) for x in range(1, 9)]
        perm_deny = [str(1) if Permission.get_by_pk(x) in self.deny.all() else str(0) for x in range(1, 9)]
        
        return "[{}, {}]".format(
            "".join(perm_allow), "".join(perm_deny)
        )

class ProfileBase(Model):    
    def has_full_access(self, entry):
        return self.has_permission(entry, PermissionEnum.FULL_ACCESS)
    
    def can_read(self, entry):
        return self.has_permission(entry, PermissionEnum.READ_NORMAL)
    
    def can_read_delete(self, species, entry):
        return self.has_permission(entry, PermissionEnum.READ_DELETE)

    def can_read_permission(self, species, entry):
        return self.has_permission(entry, PermissionEnum.READ_PERMISSION)
    
    def can_read_history(self, species, entry):
        return self.has_permission(entry, PermissionEnum.READ_HISTORY)

    def can_write(self, species, entry):
        return self.has_permission(entry, PermissionEnum.WRITE_NORMAL)
    
    def can_delete(self, species, entry):
        return self.has_permission(entry, PermissionEnum.WRITE_DELETE)

    def can_write_permission(self, species, entry):
        return self.has_permission(entry, PermissionEnum.WRITE_PERMISSION)

    class Meta:
        abstract = True

class GroupProfile(ProfileBase):
    """
    Additional information for groups
    """
    group = OneToOneField(Group, related_name = "profile")
    description = CharField(max_length=255, blank=True, default='', verbose_name = "Group description")
    
    def get_permissions(self, entry):
        """
        Returns the allow/deny list for an entry
        
        :param entry: Entry to get permissions from
        :type entry: Entry
        
        :returns: Tuple (allow, deny). (None, None) if no permission object
         was found.
        """
        try:
            group_perm = entry.group_permissions.get(group = self)
            return group_perm.allow.all(), group_perm.deny.all()
        except ObjectDoesNotExist:
            return None, None

    def has_permission(self, entry, permissions):
        """
        Checks if the group has allow or deny permissions set for a list
        specified.
        
        :param entry: Entry to check permissions against
        :type entry: Entry
        
        :param permissions: List containing the permissions to check
        :type permissions: PermissionEnum
        
        :return: Tuple (allow, deny). allow is True if the group has allow
         permission for the list specified (same for deny). (None, None)
         if no permission object was found
        """
        allow, deny = self.get_permissions(entry)
        if allow == None:
            return None, None
        allow = all(perm in permissions for perm in allow)
        deny = all(perm in permissions for perm in deny)
        return allow, deny
    
    def __unicode__(self):
        return self.group.name

class UserProfile(ProfileBase):
    """
    Additional information for users
    """
    user = OneToOneField(User, related_name = "profile")
    affiliation = CharField(max_length=255, blank=True, default='', verbose_name='Affiliation')
    website = URLField(max_length=255, blank=True, default='', verbose_name='Website')
    phone = CharField(max_length=255, blank=True, default='', verbose_name='Phone')
    address = CharField(max_length=255, blank=True, default='', verbose_name='Address')
    city = CharField(max_length=255, blank=True, default='', verbose_name='City')
    state = CharField(max_length=255, blank=True, default='', verbose_name='State')
    zip = CharField(max_length=255, blank=True, default='', verbose_name='Zip')
    country = CharField(max_length=255, blank=True, default='', verbose_name='Country')
    
    def get_name(self):
        return self.user.first_name + " " + self.user.last_name
    
    def get_website(self):
        if not self.website.startswith("http://") and not self.website.startswith("https://"):
            return "http://" + self.website
        return self.website
    
    def get_groups(self):
        """
        Returns all groups where the user is in, including "Everybody" and
        (if the user isn't a guest) "Registred"
        
        :returns: GroupProfile List
        """
        groups = self.user.groups.all()
        profiles = map(lambda g: GroupProfile.objects.get(group = g), groups)
        if self.user.username != "guest":
            profiles += [GroupProfile.objects.get(group__name = "Registred")]
        
        profiles += [GroupProfile.objects.get(group__name = "Everybody")]
        return profiles
    
    def get_permissions(self, entry):
        """
        Calculates the complete allow and deny list for an entry.
        This includes the permissions for groups the user is in including
        group "Everybody" and (except user guest) "Registred"
        
        :param entry: Entry to get permissions from
        :type entry: Entry
        
        :returns: Tuple (allow, deny). (None, None) if no permission object
         was found.
        """
        # Use None here so it's possible to distinguish between 0 permission
        # and no permission at all
        allow = None
        deny = None
        try:
            user_perm = entry.user_permissions.get(user = self)
            allow = list(user_perm.allow.all())
            deny = list(user_perm.deny.all())
        except ObjectDoesNotExist:
            pass
        
        groups = self.get_groups()
        for g in groups:
            try:
                group_perm = entry.group_permissions.get(group = g)
                if allow == None:
                    allow = []
                    deny = []
                allow += list(group_perm.allow.all())
                deny += list(group_perm.deny.all())
            except ObjectDoesNotExist:
                pass
        
        return allow, deny
    
    def has_permission(self, entry, permissions):
        """
        Checks if the user has allow or deny permissions set for a list
        specified. This includes the permissions for groups the user is in
        including group "Everybody" and (except user guest) "Registred".
        Deny usually has a higher priority then Allow. If Deny is set 
        permission should be denied by the implementation.
        
        :param entry: Entry to check permissions against
        :type entry: Entry
        
        :param permissions: List containing the permissions to check
        :type permissions: PermissionEnum
        
        :return: Tuple (allow, deny). allow is True if the user has allow
         permission for the bitmask specified (same for deny). (None, None)
         if no permission object was found.
        """
        
        
        if isinstance(entry, Entry): 
            allow, deny = self.get_permissions(entry)
            if allow == None:
                return None, None

            allow = permissions in allow
            deny = permissions in deny
            
            return allow, deny
        elif isinstance(entry, list):
            # Mass permission check
            return 0
        
        return 0, 0
    
    def is_admin(self):
        """
        Checks if the user is in the administrator group
        """
        admin = Group.objects.get(name = "Administrator")
        return self.user.groups.filter(pk = admin.pk).exists()
    
    def __unicode__(self):
        return self.user.username

    class Meta:
        verbose_name='User profile'
        verbose_name_plural = 'User profiles'
        ordering = ['user__last_name', 'user__first_name']
        get_latest_by = 'user__date_joined'

''' BEGIN: helper models '''

class TableMeta(Model):
    table_name = CharField(max_length=255, verbose_name = "Name of the table", unique=True)
    model_name = CharField(max_length=255, verbose_name = "Name of the model associated with the table", unique=True)
    
    @staticmethod
    def get_by_table_name(name):
        model_type_key = "model/model_type/tbl/" + name
        cache_model_type = Cache.try_get(model_type_key, lambda: TableMeta.objects.get(table_name = name))
        return cache_model_type
        
    @staticmethod
    def get_by_model_name(name):
        model_type_key = "model/model_type/" + name
        cache_model_type = Cache.try_get(model_type_key, lambda: TableMeta.objects.get(model_name = name))
        return cache_model_type
    
    @staticmethod
    def get_by_model(model):
        return TableMeta.get_by_model_name(model._meta.object_name)
    
    def __unicode__(self):
        return self.table_name + " - " + self.model_name

class TableMetaColumn(Model):
    table = ForeignKey(TableMeta, related_name = 'columns')
    column_name = CharField(max_length=255, verbose_name = "Name of the column")
    column_id = IntegerField(verbose_name = "Index of the column")
    
    @staticmethod
    def get_by_field(field):
        model_type_key = "column/%s/%s" % (field.model._meta.object_name, field.column)
        
        return Cache.try_get(model_type_key,
                lambda: TableMetaColumn.objects.prefetch_related("table").get(table__model_name = field.model._meta.object_name, column_name = field.column))

    def __unicode__(self):
        return "{}: {} ({})".format(self.table.model_name, self.column_name, self.column_id)

    class Meta:
        unique_together = ('table', 'column_name')

class TableMetaManyToMany(Model):
    m2m_table = ForeignKey(TableMeta, related_name = '+', verbose_name = "Data about the M2M-Table", unique = True)
    source_table = ForeignKey(TableMeta, related_name = 'm2ms_source', verbose_name = "The table the m2m references from")
    target_table = ForeignKey(TableMeta, related_name = 'm2ms_target', verbose_name = "The table the m2m references to")
    
    @staticmethod
    def get_by_m2m_table_name(name):
        model_type_key = "model/model_type_m2m/tbl/" + name
        cache_model_type = Cache.try_get(model_type_key, lambda: TableMetaManyToMany.objects.prefetch_related("m2m_table", "source_table", "target_table").get(m2m_table__table_name = name))
        return cache_model_type
        
    @staticmethod
    def get_by_m2m_model_name(name):
        model_type_key = "model/model_type_m2m/" + name
        cache_model_type = Cache.try_get(model_type_key, lambda: TableMetaManyToMany.objects.prefetch_related("m2m_table", "source_table", "target_table").get(m2m_table__model_name = name))
        return cache_model_type

    def __unicode__(self):
        return self.m2m_table.table_name + " - " + self.m2m_table.model_name
    
    class Meta:
        unique_together = ('m2m_table', 'source_table', 'target_table')

class RevisionDetail(Model):
    user = ForeignKey(UserProfile, verbose_name = "Modified by", related_name = '+', editable = False)
    date = DateTimeField(default=datetime.now, verbose_name = "Modificiation date")
    reason = CharField(max_length=255, blank=True, default='', verbose_name='Reason for edit')

class Revision(Model):
    """To allow reverting of edits all edit operations are stored in this table.
    
    Always only the new value of a single cell is stored to save memory.
    
    :Columns:
        * ``current``: Reference to the latest entry
        * ``detail``: Contains additional information about the edit 
        * ``action``: Type of the operation (Update, Insert, Delete)
        * ``column``: Table and column where the modification occured
        * ``new_value``: New value in the cell of the column
    """
    current = ForeignKey("Entry", verbose_name = "Current version of entry", related_name = 'revisions', editable = False)
    detail = ForeignKey(RevisionDetail, verbose_name = 'Details about operation', related_name = 'revisions', editable = False, null = True)
    action = CharField(max_length=1, choices=CHOICES_DBOPERATION)
    column = ForeignKey(TableMetaColumn, verbose_name = 'Table and column where the modification occured', related_name = '+')
    new_value = TextField(blank = True, null = True)

    def __unicode__(self):
        return str(self.current)
    
class RevisionManyToMany(Model):
    """To allow reverting of edits all edit operations are stored in this table.
    
    Always only the new value of a single cell is stored to save memory.
    
    :Columns:
        * ``current``: Reference to the latest entry
        * ``detail``: Contains additional information about the edit 
        * ``action``: Type of the operation (Only insert and delete)
        * ``table``: Table where the modification occured
        * ``new_value``: New value in the cell of the column
    """
    current = ForeignKey("Entry", verbose_name = "Current version of entry points to source of m2m", related_name = 'revisions_m2m', editable = False)
    detail = ForeignKey(RevisionDetail, verbose_name = 'Details about operation', related_name = 'revisions_m2m', editable = False, null = True)
    action = CharField(max_length=1, choices=CHOICES_DBOPERATION)
    table = ForeignKey(TableMetaManyToMany, verbose_name = 'M2M Table where the modification occured', related_name = '+')
    new_value = IntegerField(blank = True, null = True)

    def __unicode__(self):
        return str(self.current)

class Evidence(Model):
    value = TextField(blank=True, default='', verbose_name='Value')    
    units = CharField(max_length=255, blank=True, null=True, verbose_name='Units')
    is_experimentally_constrained = BooleanField(verbose_name='Is experimentally <br/>constrained')
    species = CharField(max_length=255, blank=True, null=True, verbose_name='Species')
    media = CharField(max_length=255, blank=True, null=True, verbose_name='Media')
    pH = FloatField(blank=True, null=True, verbose_name='pH')
    temperature = FloatField(blank=True, null=True, verbose_name='Temperature (C)')
    comments = TextField(blank=True, default='', verbose_name='Comments')
    references = ForeignKey("PublicationReference", blank=True, null=True, related_name='evidence', verbose_name='References')
        
    species_component = ManyToManyField('SpeciesComponent', blank=True, verbose_name='Species component', related_name='+')
    
    def __unicode__(self):
        arr = []
        for field in self.__class__._meta.fields:
            if not field.auto_created and getattr(self, field.name) is not None and not (isinstance(getattr(self, field.name), (str, unicode, )) and getattr(self, field.name) == ''):
                arr.append('%s: %s' % (field.name, unicode(getattr(self, field.name))))
        for field in self.__class__._meta.many_to_many:
            arr.append('%s: %s' % (field.name, unicode(getattr(self, field.name).all())))
        return ', '.join(arr)
    
    class Meta:
        ordering = ['value', 'units']
        verbose_name = 'Evidence'
        verbose_name_plural = 'Evidence'
        
class EntryData(Model):
    def __unicode__(self):
        arr = []
        txt = unicode('')
        nFields = 0        
        for field in self.__class__._meta.fields:
            if not field.auto_created:
                nFields += 1
                try:
                    txt = unicode(getattr(self, field.name))
                    arr.append('%s: %s' % (field.name, txt))
                except ObjectDoesNotExist:
                    pass
        for field in self.__class__._meta.many_to_many:
            nFields += 1
            txt = unicode(getattr(self, field.name).all())
            arr.append('%s: %s' % (field.name, txt))
                    
        if nFields == 1:
            return txt
        else:
            return ', '.join(arr)
        
    class Meta:
        abstract = True
        
class EvidencedEntryData(EntryData):
    evidence = ManyToManyField(Evidence, blank=True, null=True, verbose_name='Evidence')
    detail = ForeignKey(RevisionDetail, verbose_name = 'Last edit', related_name = '+', editable = False)
    
    class Meta:
        abstract = True
        
#Basic entry data
class EntryBooleanData(EvidencedEntryData):
    value = BooleanField(verbose_name='Value')
    
    class Meta:
        ordering = ['value']
        verbose_name = 'Entry Boolean data'
        verbose_name_plural = 'Entry Boolean data'    

class EntryCharData(EvidencedEntryData):
    value = CharField(max_length = 255, blank=True, default='', verbose_name='Value')
    units = CharField(max_length = 255, blank=True, default='', verbose_name='Units')
    
    class Meta:
        ordering = ['value', 'units']
        verbose_name = 'Entry char data'
        verbose_name_plural = 'Entry char data'    
                
class EntryFloatData(EvidencedEntryData):
    value = FloatField(verbose_name='Value')
    units = CharField(max_length = 255, blank=True, default='', verbose_name='Units')
    
    class Meta:
        ordering = ['value', 'units']
        verbose_name = 'Entry float data'
        verbose_name_plural = 'Entry float data'
        
class EntryPositiveFloatData(EvidencedEntryData):
    value = FloatField(verbose_name='Value', validators=[validators.MinValueValidator(0)])
    units = CharField(max_length = 255, blank=True, default='', verbose_name='Units')
    
    class Meta:
        ordering = ['value', 'units']
        verbose_name = 'Entry positive float data'
        verbose_name_plural = 'Entry positive float data'

class EntryTextData(EvidencedEntryData):
    value = TextField(blank=True, default='', verbose_name='Value')
    units = CharField(max_length = 255, blank=True, default='', verbose_name='Units')
    
    class Meta:
        ordering = ['value', 'units']
        verbose_name = 'Entry text data'
        verbose_name_plural = 'Entry text data'            
        
#Complex entry data
class BindingSite(EvidencedEntryData):
    coordinate = PositiveIntegerField(verbose_name='Coordinate (nt)')
    length = PositiveIntegerField(verbose_name='Length (nt)')
    direction = CharField(max_length=10, choices=CHOICES_DIRECTION, verbose_name='Direction')    
            
    #getter
    def get_chromosome(self):
        return self.transcriptional_regulations.all()[0].transcription_unit.get_chromosome()
    
    class Meta:
        ordering = ['coordinate', 'length']
        verbose_name = 'Binding site'
        verbose_name_plural = 'Binding sites'    
        
class BiomassComposition(EvidencedEntryData):
    concentration = FloatField(verbose_name='Concentration (mmol gDCW<sup>-1</sup>)', validators=[validators.MinValueValidator(0)])
    compartment = ForeignKey('Compartment', related_name='biomass_compositions', verbose_name='Compartment')
    
    class Meta:
        ordering = ['-concentration']
        verbose_name = 'Biomass composition'
        verbose_name_plural = 'Biomass composition'        

class Codon(EvidencedEntryData):
    sequence = CharField(max_length=3, verbose_name='Sequence', validators=[
        validate_dna_sequence,
        validators.MinLengthValidator(3),
        ])
    
    class Meta:
        ordering = ['sequence']
        verbose_name='Codon'
        verbose_name_plural = 'Codons'
        
class CoenzymeParticipant(EvidencedEntryData):
    metabolite = ForeignKey('Metabolite', related_name='coenzyme_participants', verbose_name='Metabolite')
    compartment = ForeignKey('Compartment', related_name='+', verbose_name='Compartment')
    coefficient = FloatField(blank=True, null=True, verbose_name='Coefficient', validators=[validators.MinValueValidator(0)])
    
    class Meta:
        ordering = []
        verbose_name='Coenzyme participant'
        verbose_name_plural = 'Coenzyme participants'

class DisulfideBond(EvidencedEntryData):
    protein_monomer = ForeignKey('ProteinMonomer', related_name = 'disulfide_bonds', verbose_name='Protein monomer')
    residue_1 = PositiveIntegerField(verbose_name='Residue-1')
    residue_2 = PositiveIntegerField(verbose_name='Residue-2')
    
    class Meta:
        ordering = []
        verbose_name = 'Disulfide bond'
        verbose_name_plural = 'Disulfide bonds'        
        
class DNAFootprint(EvidencedEntryData):
    length = PositiveIntegerField(null=True, blank=True, verbose_name='Length (nt)')
    binding = CharField(max_length=10, blank=True, default='', verbose_name='Binding', choices=CHOICES_STRANDEDNESS)
    region = CharField(max_length=10, blank=True, default='', verbose_name='Region', choices=CHOICES_STRANDEDNESS)    
    
    class Meta:
        ordering = []
        verbose_name = 'DNA footprint'
        verbose_name_plural = 'DNA footprints'
        
class EnzymeParticipant(EvidencedEntryData):
    protein = ForeignKey('Protein', related_name='enzyme_participants', verbose_name='Protein')
    compartment = ForeignKey('Compartment', related_name='+', verbose_name='Compartment')
    
    class Meta:
        ordering = []
        verbose_name='Enzyme participant'
        verbose_name_plural = 'Enzyme participants'
        
class Homolog(EvidencedEntryData):
    xid = CharField(max_length=255, verbose_name='External ID')
    species = CharField(max_length=20, choices=CHOICES_HOMOLOG_SPECIES, verbose_name='Species')
    
    class Meta:
        ordering = ['xid']
        verbose_name='Homolog'
        verbose_name_plural = 'Homologs'
        
class Kinetics(EvidencedEntryData):
    rate_law = CharField(blank=True, default='', max_length=255, verbose_name='Rate law')
    km = CharField(max_length=255, blank=True, verbose_name='K<sub>m</sub> (&mu;M)', validators=[validators.RegexValidator(r'^([0-9\.]+)(, [0-9\.]+)*$')])
    vmax = FloatField(blank=True, null=True, verbose_name='V<sub>max</sub>', validators=[validators.MinValueValidator(0)])
    vmax_unit = CharField(blank=True, max_length=255, choices=CHOICES_VMAX_UNITS, verbose_name='V<sub>max</sub> Unit')
    
    def get_vmax_normalized(self):
        from cyano.helpers import getModel
        
        if vmax_unit.name  == 'U/mg':
            enz = self.reactions.all()[0].enzyme
            enz = getModel(enz.model_type).objects.get(id=enz.id)            
            return self.vmax * enz.get_molecular_weight() * 1e-3
        else:
            return self.vmax
    
    class Meta:
        ordering = []
        verbose_name='Kinetics'
        verbose_name_plural = 'Kinetics'
    
class MediaComposition(EvidencedEntryData):
    concentration = FloatField(verbose_name='Concentration (mM)', validators=[validators.MinValueValidator(0)])
    is_diffused = BooleanField(verbose_name='Is diffused')
    
    class Meta:
        ordering = ['-concentration']
        verbose_name='Media composition'
        verbose_name_plural = 'Media composition'
        
class MetaboliteMapCoordinate(EntryData):
    compartment = ForeignKey('Compartment', related_name='+', verbose_name='Compartment')
    x = FloatField(verbose_name='X')
    y = FloatField(verbose_name='Y')
    
    class Meta:
        ordering = ['x', 'y', 'compartment']
        verbose_name='Metabolite map coordinate'
        verbose_name_plural = 'Metabolite map coordinates'
        
class ModificationReaction(EvidencedEntryData):
    molecule = ForeignKey('Molecule', related_name='modification_reactions', verbose_name='Molecule')
    compartment = ForeignKey('Compartment', related_name='+', verbose_name='Compartment')
    position = PositiveIntegerField(blank=True, null=True, verbose_name='Position')
    
    class Meta:
        ordering = []
        verbose_name='Protein monomer'
        verbose_name_plural = 'Protein monomers'
        
class ProstheticGroupParticipant(EvidencedEntryData):
    metabolite = ForeignKey('Metabolite', related_name='prosthetic_group_participants', verbose_name='Metabolite')
    compartment = ForeignKey('Compartment', related_name='+', verbose_name='Compartment')
    coefficient = PositiveIntegerField(blank=True, null=True, verbose_name='Coefficient')
    
    class Meta:
        ordering = []
        verbose_name='Prosthetic group participant'
        verbose_name_plural = 'Prosthetic group participants'
        
class ProteinComplexBiosythesisParticipant(EvidencedEntryData):
    molecule = ForeignKey('Molecule', related_name='protein_complex_biosythesis_participants', verbose_name='Molecule')
    residue = PositiveIntegerField(blank=True, null=True, verbose_name='Residue')
    coefficient = FloatField(verbose_name='Coefficient')
    compartment = ForeignKey('Compartment', related_name='+', verbose_name='Compartment')    
    
    class Meta:
        ordering = []
        verbose_name='Protein complex biosythesis participant'
        verbose_name_plural = 'Protein complex biosythesis participants'
        
class ReactionMapCoordinate(EntryData):
    path = TextField(verbose_name='Path')
    value_x = FloatField(verbose_name='Label-X')
    value_y = FloatField(verbose_name='Label-Y')
    label_x = FloatField(verbose_name='Value-X')
    label_y = FloatField(verbose_name='Value-Y')
    
    class Meta:
        ordering = ['value_x', 'value_y', 'label_x', 'label_y']
        verbose_name='Reaction map coordinate'
        verbose_name_plural = 'Reaction map coordinates'
        
class ReactionStoichiometryParticipant(EvidencedEntryData):
    molecule = ForeignKey('Molecule', related_name='reaction_stoichiometry_participants', verbose_name='Molecule')
    coefficient = FloatField(verbose_name='Coefficient')
    compartment = ForeignKey('Compartment', related_name='+', verbose_name='Compartment')
    
    class Meta:
        ordering = []
        verbose_name='Molecule coefficient compartment'
        verbose_name_plural = 'Molecule coefficient compartments'
        
class SignalSequence(EvidencedEntryData):
    type = CharField(max_length=20, choices=CHOICES_SIGNAL_SEQUENCE_TYPE, verbose_name='Type')
    location = CharField(max_length=1, choices=CHOICES_SIGNAL_SEQUENCE_LOCATION, verbose_name='Location')
    length = PositiveIntegerField(verbose_name='Length (nt)')
    
    class Meta:
        ordering = ['type', 'location', 'length']
        verbose_name = 'Signal sequence'
        verbose_name_plural = 'Signal sequences'
        
class Synonym(EntryData):
    name = CharField(max_length=255, verbose_name='Name')
    
    class Meta:
        ordering = ['name']
        verbose_name = 'Synonym'
        verbose_name_plural = 'Synonyms'

''' END: helper models '''

''' BEGIN: Base classes for all knowledge base objects '''
@receiver(m2m_changed)
def m2m_changed_save(sender, instance, action, reverse, model, pk_set, **kwargs):
    if reverse:
        print "WARN: Reverse support not implemented:",sender,instance,pk_set
        return
    
    if not pk_set is None and len(pk_set) > 0 and isinstance(instance, Entry):
        if action == "post_add":
            print "m2m-add:",TableMetaManyToMany.get_by_m2m_model_name(sender._meta.object_name),instance,pk_set

            table_meta = TableMetaManyToMany.get_by_m2m_model_name(sender._meta.object_name)
    
            save_list = []
    
            for pk in pk_set:
                # FIXME: detail_id is the current one
                save_list.append(RevisionManyToMany(current_id = instance.pk, detail_id = instance.detail_id, action = "I", table = table_meta, new_value = pk))
    
            if len(save_list) > 0:
                print instance.wid + ": revisioning", len(save_list), "items"
                # Don't update the detail when actually nothing changed for that entry
                RevisionManyToMany.objects.bulk_create(save_list)
        elif action == "post_delete":
            print "m2m-del:",TableMetaManyToMany.get_by_m2m_model_name(sender._meta.object_name),instance,pk_set

            table_meta = TableMetaManyToMany.get_by_m2m_model_name(sender._meta.object_name)
    
            save_list = []
    
            for pk in pk_set:
                # FIXME: detail_id is the current one
                save_list.append(RevisionManyToMany(current_id = instance.pk, detail_id = instance.detail_id, action = "D", table = table_meta, new_value = pk))
    
            if len(save_list) > 0:
                print instance.wid + ": deleting", len(save_list), "items"
                # Don't update the detail when actually nothing changed for that entry
                RevisionManyToMany.objects.bulk_create(save_list)
            pass

class Entry(Model):
    """Base class for all knowledge base objects.
    
    :Columns:
        * ``model_type``: Specifies the child type of the item (used to find the child table for inheritance)
        * ``wid``: Warehouse Identifier. Unique identifier in the warehouse, used in the URLs e.g.
        * ``name``: Short, human readable name of that entry
        * ``synonyms``: Other names for that entry
        * ``comments``: Detailed notes about this entry
    """
    model_type = ForeignKey(TableMeta)
    wid = SlugField(max_length=150, unique = True, verbose_name='WID', validators=[validators.validate_slug])
    name = CharField(max_length=255, blank=True, default='', verbose_name='Name')
    synonyms = ManyToManyField(Synonym, blank=True, null=True, related_name='entry', verbose_name='Synonyms')
    comments = TextField(blank=True, default='', verbose_name='Comments')

    created_detail = ForeignKey(RevisionDetail, verbose_name='Entry created revision', related_name='entry_created_detail')
    detail = ForeignKey(RevisionDetail, verbose_name='Last Revision entry', related_name='entry_detail')
    
    def __unicode__(self):
        return self.wid
    
    def natural_key(self):
        return self.wid

    def save(self, *args, **kwargs):
        # Optimized to reduce number of database accesses to a minimum
        
        from itertools import ifilter
        from cyano.helpers import slugify

        if self.wid != slugify(self.wid):
            # Slugify the WID
            self.wid = slugify(self.wid)

        rev_detail = kwargs["revision_detail"]
        # Prevent error on save later
        del kwargs["revision_detail"]

        old_item = None
        if self.pk != None:
            # can this be optimized? Probably not
            old_item = self._meta.concrete_model.objects.get(pk = self.pk)
            super(Entry, self).save(*args, **kwargs)
        else:
            # New entry (no primary key)
            # The latest entry is not revisioned to save space (and time)
            model_type_key = "model/model_type/" + str(self._meta.object_name)
            cache_model_type = Cache.try_get(model_type_key, lambda: TableMeta.objects.get(model_name = self._meta.object_name))
            self.created_detail = rev_detail
            self.detail = rev_detail
            self.model_type = cache_model_type

            ##print "CREATE: No revision needed for", str(self.wid)
            super(Entry, self).save(*args, **kwargs)
            return

        save_list = []

        fields = self._meta.fields
        # Remove primary keys and some fields that don't need revisioning:
        fields = ifilter(lambda x: (not x.primary_key) and
                         (not x.name in ["detail", "created_detail"]), fields)

        for field in fields:
            if hasattr(self, field.name + "_id"):
                new_value = getattr(self, field.name + "_id")
            else:
                new_value = getattr(self, field.name)

            if old_item == None:
                old_value = None
            else:
                if hasattr(old_item, field.name + "_id"):
                    old_value = getattr(old_item, field.name + "_id")
                else:
                    old_value = getattr(old_item, field.name)

            if old_value != new_value:
                column = TableMetaColumn.get_by_field(field)

                if new_value == None:
                    ##print "delete"
                    action = "D"
                else:
                    if Revision.objects.filter(current_id = self.pk, column = column).exists():
                        ##print "update", self.wid, tm.model_name, column
                        action = "U"
                    else:
                        ##print "insert", self.wid, tm.model_name, column
                        action = "I"

                # Assumption: When nothing changed saving isn't needed at all
                save_list.append(Revision(current_id = self.pk, detail_id = old_item.created_detail_id if action == "I" else old_item.detail_id, action = action, column = column, new_value = old_value))

        if len(save_list) > 0:
            ##print self.wid + ": revisioning", len(save_list), "items"
            # Don't update the detail when actually nothing changed for that entry
            self.detail = rev_detail
            Revision.objects.bulk_create(save_list)
    
    #html formatting
    def get_as_html_synonyms(self, species, is_user_anonymous):
        return format_list_html([x.name for x in self.synonyms.all()], comma_separated=True)
    
    def get_as_html_cross_references(self, species, is_user_anonymous):
        results = []
        for cr in self.cross_references.all():
            results.append('%s: <a href="%s">%s</a>' % (cr.source, reverse("db_xref.views.dbxref", kwargs={"source" : cr.source, "xid" : cr.xid}), cr.xid ))
        return format_list_html(results, separator=', ')
    
    def get_as_html_created_user(self, species, is_user_anonymous):
        if is_user_anonymous:
            return '%s' % (self.created_date.strftime("%Y-%m-%d %H:%M:%S"))
        else:
            detail = self.created_detail
            user = detail.user.user
            return '<a href="%s">%s %s</a> on %s' % (user.get_absolute_url(), user.first_name, user.last_name, detail.date.strftime("%Y-%m-%d %H:%M:%S"))
    
    def get_as_html_last_updated_user(self, species, is_user_anonymous):
        if is_user_anonymous:
            return '%s' % (self.last_updated_date.strftime("%Y-%m-%d %H:%M:%S"))
        else:
            detail = self.detail
            user = detail.user.user
            return '<a href="%s">%s %s</a> on %s' % (user.get_absolute_url(), user.first_name, user.last_name, detail.date.strftime("%Y-%m-%d %H:%M:%S"))        
    
    #meta information
    class Meta:
        concrete_entry_model = False
        fieldsets = [
            ('Type', {'fields': ['model_type']}), 
            ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), 
            ('Comments', {'fields': ['comments']}), 
            ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]}),
            ]
        field_list = [
            'id', 'wid', 'name', 'synonyms', 'cross_references', 'comments', 
            ]
        listing = ['wid', 'name']
        facet_fields = []
        ordering = ['wid']
        get_latest_by = 'createdDate'
        verbose_name = 'Entry'
        verbose_name_plural = 'Entries'

class CrossReferenceMeta(Model):
    name = CharField(max_length=255, verbose_name = "Identifier for crossreference (DB name or URL)", editable = False)

class SpeciesComponent(Entry):
    '''
    Contains all components that belong to an organism.
    Improves the lookup speed and allows inheritance.
    
    * ``cross_references``: Databases referencing that entry
    * ``publication_references``: Publications referencing that entry
    '''
    species = ManyToManyField("Species", verbose_name = "Organism containing that entry", related_name = "species")
    type = ManyToManyField('Type', blank=True, null=True, related_name='members', verbose_name='Type')
    cross_references = ManyToManyField("CrossReference", blank=True, null=True, related_name='cross_referenced_entries', verbose_name='Cross references')
    publication_references = ManyToManyField("PublicationReference", blank=True, null=True, related_name='publication_referenced_entries', verbose_name='Publication references')

    #getters
    @permalink
    def get_absolute_url(self, species, history_id = None):
        model_type_key = "model/model_type/" + str(self.model_type_id)
        cache_model_type = Cache.try_get(model_type_key, lambda: self.model_type)
        
        dic = {'species_wid':species.wid, 'model_type':cache_model_type.model_name, 'wid': self.wid}
        
        if history_id is None:
            view = 'cyano.views.detail'
        else:
            view = 'cyano.views.history_detail'
            dic['detail_id'] = history_id

        return (view, (), dic)
        
    def get_all_references(self):
        return self.publication_references.all() | PublicationReference.objects.filter(evidence__species_component__id = self.id)
        
    #html formatting
    def get_as_html_parameters(self, species, is_user_anonymous):
        results = []
        for p in self.parameters.all():
            results.append('<a href="%s">%s</a>: <i>%s</i> = %s %s' % (p.get_absolute_url(species), p.wid, p.name, p.value.value, p.value.units))
        return format_list_html(results)
        
    def get_as_html_comments(self, species, is_user_anonymous):
        txt = self.comments
        
        #provide links to references
        return re.sub(r'\[(PUB_\d{4,4})(, PUB_\d{4,4})*\]', 
            lambda match: '[' + ', '.join(['<a href="%s">%s</a>' % (reverse('cyano.views.detail', kwargs={'species_wid':self.species.wid, 'wid': x}), x, ) for x in match.group(0)[1:-1].split(', ')]) + ']',
            txt)
    
    def get_as_html_references(self, species, is_user_anonymous):
        results = {}
        for r in self.get_all_references():
            key = r.authors + ' ' + r.editors
            results[key] = r.get_citation(True)
            
        keys = results.keys()
        keys.sort()
        ordered_results = []
        for key in keys:
            ordered_results.append(results[key])
        return format_list_html(ordered_results, numbered=True, force_list=True)
    
    #meta information
    class Meta:
        concrete_entry_model = False
        fieldsets = [
            ('Type', {'fields': ['model_type']}),
            ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), 
            ('Classification', {'fields': ['type']}), 
            ('Comments', {'fields': ['comments', 'publication_references']}),
            ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]}),
            ]
        field_list = [
            'id', 'wid', 'name', 'synonyms', 'cross_references', 'type',  'comments', 'publication_references',  'created_user', 'created_date', 'last_updated_user', 'last_updated_date', 
            ]
        facet_fields = ['type']
        verbose_name = 'Species component'
        verbose_name_plural = 'Species components'

class Molecule(SpeciesComponent):
    #parent pointer
    parent_ptr_species_component = OneToOneField(SpeciesComponent, related_name='child_ptr_molecule', parent_link=True, verbose_name='Species component')
    
    #additional fields
    
    #getters
    def get_empirical_formula(self):
        from cyano.helpers import EmpiricalFormula
        return EmpiricalFormula()
        
    def get_molecular_weight(self):
        return 0 # FIXME
        #return self.get_empirical_formula().get_molecular_weight()
        
    def get_atoms(self):
        return sum(self.get_empirical_formula().values())
        
    def get_absorbance_factor(self):
        return 1 / self.get_extinction_coefficient()
    
    #html formatting        
    def get_as_html_empirical_formula(self, species, is_user_anonymous):
        return ""# FIXME
        #return self.get_empirical_formula().get_as_html()
        
    def get_as_html_molecular_weight(self, species, is_user_anonymous):
        return ""# FIXME
        return self.get_molecular_weight()
        
    def get_as_html_extinction_coefficient(self, species, is_user_anonymous):
        return ""# FIXME
        return self.get_extinction_coefficient()
        
    def get_as_html_absorbance_factor(self, species, is_user_anonymous):
        return ""# FIXME
        return self.get_absorbance_factor()
        
    def get_as_html_pi(self, species, is_user_anonymous):
        return ""# FIXME
        if hasattr(self, 'pi'):
            return self.pi
        else:
            return self.get_pi()
            
    def get_as_html_modification_reactions(self, species, is_user_anonymous):
        return ""# FIXME
        results = []
        for obj in self.modification_reactions.all():
            if len(obj.reactions.all()) == 0:
                continue
            rxn = obj.reactions.all()[0]
            results.append('<a href="%s">%s</a><br/>%s' % (rxn.get_absolute_url(species), rxn.name, rxn.get_as_html_stoichiometry(is_user_anonymous)))
        return format_list_html(results, vertical_spacing=True)
            
    def get_as_html_reaction_stoichiometry_participants(self, species, is_user_anonymous):
        return ""# FIXME
        results = []
        for obj in self.reaction_stoichiometry_participants.all():
            if len(obj.reactions.all()) == 0:
                continue
            rxn = obj.reactions.all()[0]
            results.append('<a href="%s">%s</a><br/>%s' % (rxn.get_absolute_url(species), rxn.name, rxn.get_as_html_stoichiometry(is_user_anonymous)))
        return format_list_html(list(set(results)), vertical_spacing=True)
        
    def get_as_html_protein_complex_biosythesis_participants(self, species, is_user_anonymous):
        return ""# FIXME
        results = []
        for obj in self.protein_complex_biosythesis_participants.all():            
            if len(obj.protein_complexes.all()) == 0:
                continue
            pc = obj.protein_complexes.all()[0]
            results.append('<a href="%s">%s</a><br/>%s' % (pc.get_absolute_url(species), pc.name, pc.get_as_html_biosynthesis(is_user_anonymous)))
        return format_list_html(results)
    
    #meta information
    class Meta:
        concrete_entry_model = False
        fieldsets = [
            ('Type', {'fields': ['model_type']}),
            ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), 
            ('Classification', {'fields': ['type']}), 
            ('Function', {'fields': [
                {'verbose_name': 'Reaction participant', 'name':'reaction_stoichiometry_participants'},
                {'verbose_name': 'Complex subunit', 'name':'protein_complex_biosythesis_participants'},
                ]}),
            ('Parameters', {'fields': ['parameters']}),
            ('Comments', {'fields': ['comments', 'publication_references']}),
            ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]}),
            ]
        field_list = [
            'id', 'wid', 'name', 'synonyms', 'cross_references', 'type',  'comments', 'publication_references',  'created_user', 'created_date', 'last_updated_user', 'last_updated_date', 
            ]
        facet_fields = ['type']
        verbose_name = 'Molecule'
        verbose_name_plural = 'Molecules'
        
class Protein(Molecule):
    #parent pointer
    parent_ptr_molecule = OneToOneField(Molecule, related_name='child_ptr_protein', parent_link=True, verbose_name='Molecule')
    
    #additional fields
    prosthetic_groups = ManyToManyField(ProstheticGroupParticipant, blank=True, null=True, related_name='proteins', verbose_name='Prosthetic groups')
    chaperones = ManyToManyField('self', symmetrical=False, blank=True, null=True, related_name='chaperone_substrates', verbose_name='Chaperones')
    dna_footprint = ForeignKey(DNAFootprint, null=True, blank=True, related_name='proteins', verbose_name='DNA footprint')
    regulatory_rule = ForeignKey(EntryCharData, null=True, blank=True, on_delete=SET_NULL, verbose_name='Regulatory rule', related_name='+')    
    
    #html formatting
    def get_as_html_prosthetic_groups(self, species, is_user_anonymous):
        results = []
        for p in self.prosthetic_groups.all():
            results.append(format_with_evidence(list_item = True, obj = p, txt = '(%s) <a href="%s">%s</a> [<a href="%s">%s</a>]' % (p.coefficient, p.metabolite.get_absolute_url(species), p.metabolite.wid, p.compartment.get_absolute_url(species), p.compartment.wid)))
        return format_list_html(results, force_list=True)
        
    def get_as_html_chaperones(self, species, is_user_anonymous):
        results = [];
        for c in self.chaperones.all():
            results.append('<a href="%s">%s</a>' % (c.get_absolute_url(species), c.wid))
        return format_list_html(results, comma_separated=True)
        
    def get_as_html_chaperone_substrates(self, species, is_user_anonymous):
        results = [];
        for p in self.chaperone_substrates.all():
            results.append('<a href="%s">%s</a>' % (p.get_absolute_url(species), p.wid))
        return format_list_html(results, comma_separated=True)
        
    def get_as_html_dna_footprint(self, species, is_user_anonymous):
        if self.dna_footprint is None:
            return None
        return format_with_evidence(obj = self.dna_footprint, txt = 'Length: %s (nt), Binding: %s, Region: %s' % (self.dna_footprint.length, self.dna_footprint.binding, self.dna_footprint.region))
        
    def get_as_html_regulatory_rule(self, species, is_user_anonymous):
        if self.regulatory_rule is not None and self.regulatory_rule.value is not None and self.regulatory_rule.value != '':
            return parse_regulatory_rule(self.regulatory_rule.value, {}, self.species.wid)
        return ''
        
    def get_as_html_transcriptional_regulations(self, species, is_user_anonymous):
        results = [];
        for r in self.transcriptional_regulations.all():
            results.append('<a href="%s">%s</a>: <a href="%s">%s</a>' % (r.get_absolute_url(species), r.wid, r.transcription_unit.get_absolute_url(species), r.transcription_unit.wid))
        return format_list_html(results)
        
    def get_as_html_enzyme_participants(self, species, is_user_anonymous):
        results = [];
        for obj in self.enzyme_participants.all():
            if len(obj.reactions.all()) == 0:
                continue
            rxn = obj.reactions.all()[0]
            results.append('<a href="%s">%s</a><br/>%s' % (rxn.get_absolute_url(species), rxn.name, rxn.get_as_html_stoichiometry(is_user_anonymous)))
        return format_list_html(results, vertical_spacing=True)
        
    def get_as_html_half_life(self, species, is_user_anonymous):
        return self.get_half_life()
                    
    #meta information
    class Meta:        
        concrete_entry_model = False
        fieldsets = [
            ('Type', {'fields': ['model_type']}),
            ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), 
            ('Classification', {'fields': ['type']}), 
            ('Structure', {'fields': ['prosthetic_groups', 'chaperones', 'dna_footprint']}), 
            ('Regulation', {'fields': ['regulatory_rule']}), 
            ('Function', {'fields': [                
                {'verbose_name': 'Enzyme', 'name': 'enzyme_participants'},
                {'verbose_name': 'Transcriptional regulation', 'name': 'transcriptional_regulations'},
                {'verbose_name': 'Protein folding substrates', 'name': 'chaperone_substrates'},    
                {'verbose_name': 'Reaction participant', 'name':'reaction_stoichiometry_participants'},
                {'verbose_name': 'Complex subunit', 'name':'protein_complex_biosythesis_participants'},
                ]}), 
            ('Parameters', {'fields': ['parameters']}),
            ('Comments', {'fields': ['comments', 'publication_references']}),
            ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]}),
            ]            
        field_list = [
            'id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'prosthetic_groups', 'chaperones', 'dna_footprint', 'regulatory_rule', 'comments', 'publication_references',  'created_user', 'created_date', 'last_updated_user', 'last_updated_date', 
            ]
        facet_fields = ['type', 'chaperones', 'dna_footprint__binding', 'dna_footprint__region']
        verbose_name='Protein'
        verbose_name_plural = 'Proteins'
        
        def clean(model, obj_data, all_obj_data=None, all_obj_data_by_model=None):
            #regulatory rule
            if obj_data['regulatory_rule'] is not None and obj_data['regulatory_rule']['value'] is not None and obj_data['regulatory_rule']['value'] != '':
                parse_regulatory_rule(obj_data['regulatory_rule']['value'], all_obj_data, obj_data['species'])
            
            #DNA footprint
            if obj_data['dna_footprint'] is not None:
                chr_lens = []
                
                if all_obj_data_by_model is None:
                    species_id = Species.objects.values('id').get(wid=obj_data['species'])['id']
                    for obj in Chromosome.objects.values('length').filter(species__id=species_id):
                        chr_lens.append(obj['length'])
                else:
                    for obj in all_obj_data_by_model['Chromosome']:
                        if isinstance(obj, Entry):
                            chr_lens.append(obj.length)
                        else:
                            chr_lens.append(obj['length'])
                
                if obj_data['dna_footprint']['length'] > max(chr_lens):
                    raise ValidationError({'dna_footprint': 'Length must be less than chromosome length'})
        
''' END: base classes '''
        
''' 
BEGIN: Specific data types 
'''

class Chromosome(Molecule):
    #parent pointer
    parent_ptr_molecule = OneToOneField(Molecule, related_name='child_ptr_chromosome', parent_link=True, verbose_name='Molecule')
    
    #additional fields
    sequence = TextField(blank=True, default='', verbose_name='Sequence', validators=[validate_dna_sequence])
    length = PositiveIntegerField(verbose_name='Length (nt)')
    
    #getters
    def get_length(self):
        return self.length
        
    def get_gc_content(self):
        seq = self.sequence
        return float(seq.count('G') + seq.count('C')) / float(len(seq))
        
    def get_transcription_units(self):
        return list(set([g.transcription_units.all()[0] for g in self.genes.all()]))    
        
    #http://www.owczarzy.net/extinct.htm
    def get_extinction_coefficient(self):    
        from cyano.helpers import ExtinctionCoefficient
        
        seq = self.sequence
        
        value = 0;
        for i in range(len(seq) - 1):
            value += ExtinctionCoefficient.pairwise_dna[seq[i]][seq[i+1]]
        value += ExtinctionCoefficient.pairwise_dna[seq[-1]][seq[0]]
        return value / 2
        
    def get_pi(self):
        return calculate_nucleic_acid_pi(self.sequence)
    
    #html formatting
    def get_as_html_sequence(self, species, is_user_anonymous):
        from cyano.helpers import format_sequence_as_html
        return format_sequence_as_html(self.sequence)

    def get_as_html_diff_sequence(self, species, new_obj, is_user_anonymouse):
            from Bio import pairwise2
            pairwise2.MAX_ALIGNMENTS = 1
            align = pairwise2.align.globalxx(self.sequence[:100], new_obj[:100])
            output = StringIO.StringIO()
            
            for o, n in zip(align[0][0], align[0][1]):
                output.write("|" if o == n else " ")
            
            align_graphic = output.getvalue()
            output.close()
            output = StringIO.StringIO()
            
            old_align = re.sub(r'(.{80})', r'\1\n', align[0][0]).split("\n")
            new_align = re.sub(r'(.{80})', r'\1\n', align[0][1]).split("\n")
            align_graphic = re.sub(r'(.{80})', r'\1\n', align_graphic).split("\n")
            output.write("""<div class="sequence">""")
            output.write("Alignment of first 500 bp.")
            if self.sequence == new_obj:
                output.write(" Sequences are completly identical.")
            else:
                output.write(" Sequences differ.")
            output.write("<div></div><div>")
            for ol, ag, ne in zip(old_align, align_graphic, new_align):
                output.write(ol + "<br/>")
                output.write(ag.replace(" ", "&nbsp;") + "<br/>")
                output.write(ne + "<br/>")
                output.write("<br/>")
            output.write("</div></div>")
            
            return None, output.getvalue()
        
    def get_as_html_structure(self, species, is_user_anonymous, start_coordinate = None, end_coordinate = None, highlight_wid = None, zoom = 0):
        if zoom == 0:
            return self.get_as_html_structure_global(species, is_user_anonymous)
        else:
            return self.get_as_html_structure_local(species, is_user_anonymous, start_coordinate = start_coordinate, end_coordinate = end_coordinate, highlight_wid = highlight_wid)
    
    def get_as_html_structure_global(self, species, is_user_anonymous):
        ntPerSegment = 1e4
        segmentHeight = 27
        geneHeight = 10
        featureHeight = 5
        nSegments = int(math.ceil(self.length / ntPerSegment))
        W = 636
        H = segmentHeight * nSegments
        segmentLeft = 35
        segmentW = W - 4 - segmentLeft
        chrTop = -12
        
        #style
        colors = ['3d80b3', '3db34a', 'cd0a0a', 'e78f08']
        style = ''
        for i in range(len(colors)):
            style += '.color-%s{fill:#%s; stroke: #%s;}' % (i, colors[i], colors[i], )

        #chromosome
        chrStyle = '\
            .chr text{fill:#222; text-anchor:end; alignment-baseline:middle; font-size:10px}\
            .chr line{stroke:#666; stroke-width:0.5px;}\
        '
        
        chr = StringIO.StringIO()
        for i in range(nSegments):
            x1 = segmentLeft
            x2 = segmentLeft + ((min(self.length, (i+1) * ntPerSegment) - 1) % ntPerSegment) / ntPerSegment * segmentW
            y = chrTop + (i + 1) * segmentHeight
            chr.write('<text x="%s" y="%s">%d</text>' % (segmentLeft - 2, y, i * ntPerSegment + 1))
            chr.write('<line x1="%s" x2="%s" y1="%s" y2="%s"/>' % (x1, x2, y, y))
        
        #genes
        geneStyle = '\
            .genes g polygon{stroke-width:1px; fill-opacity:0.5;}\
            .genes g text{text-anchor:middle; alignment-baseline:middle; font-size:8px; fill: #222}\
        '

        genes = StringIO.StringIO()


        genesList = list(self.genes.prefetch_related('transcription_units', 'transcription_units__transcriptional_regulations').all())

        nTus = 0
        iTUs = {}
        tus = []

        for i, gene in enumerate(genesList):
            if gene.transcription_units.count() > 0:
                tu = gene.transcription_units.all()[:1][0]
                if iTUs.has_key(tu.wid):
                    iTu = iTUs[tu.wid]
                else:
                    tus.append(tu)
                    iTu = nTus
                    iTUs[tu.wid] = iTu
                    nTus += 1
                    
                iSegment = math.floor((gene.coordinate - 1) / ntPerSegment)
                
                if gene.direction == 'f':
                    x1 = segmentLeft + ((gene.coordinate - 1) % ntPerSegment) / ntPerSegment * segmentW
                    x3 = min(segmentLeft + segmentW, x1 + gene.length / ntPerSegment * segmentW)
                    x2 = max(x1, x3 - 5)
                else:
                    x3 = segmentLeft + ((gene.coordinate - 1) % ntPerSegment) / ntPerSegment * segmentW
                    x1 = min(segmentLeft + segmentW, x3 + gene.length / ntPerSegment * segmentW)
                    x2 = min(x1, x3 + 5)
                    
                y2 = chrTop + (iSegment + 1) * segmentHeight - 2
                y1 = y2 - geneHeight
                
                if math.fabs(x3 - x1) > len(gene.wid) * 5:
                    label = gene.wid
                else:
                    label = ''
                    
                if gene.name:
                    tip_title = gene.name
                else:
                    tip_title = gene.wid
                tip_content = 'Transcription unit: %s' % tu.name
                tip_title = tip_title.replace("'", "\'")
                tip_content = tip_content.replace("'", "\'")
                    
                gene_abs_url = gene.get_absolute_url(species)
                
                genes.write('<g>\
                    <a xlink:href="%s">\
                        <polygon class="color-%s" points="%s,%s %s,%s %s,%s %s,%s %s,%s" onmousemove="javascript: showToolTip(evt, \'%s\', \'%s\')" onmouseout="javascript: hideToolTip(evt);"/>\
                    </a>\
                    <a xlink:href="%s">\
                        <text x="%s" y="%s" onmousemove="javascript: showToolTip(evt, \'%s\', \'%s\')" onmouseout="javascript: hideToolTip(evt);">%s</text>\
                    </a>\
                    </g>' % (                
                    gene_abs_url,
                    iTu % len(colors),
                    x1, y1,
                    x2, y1,
                    x3, (y1 + y2) / 2, 
                    x2, y2,
                    x1, y2,
                    tip_title, tip_content,
                    gene_abs_url,
                    (x1 + x3) / 2, (y1 + y2) / 2 + 1, 
                    tip_title, tip_content,
                    label,                
                    ))
        
        #promoters
        promoterStyle = '.promoters rect{fill:#%s; opacity:0.5}' % (colors[0], )
        tfSiteStyle = '.tfSites rect{fill:#%s; opacity:0.5}' % (colors[1], )
        
        promoters = StringIO.StringIO()
        tfSites = StringIO.StringIO()
        for tu in tus:
            if tu.promoter_35_coordinate is not None:
                tu_coordinate = tu.get_coordinate() + tu.promoter_35_coordinate
                tu_length = tu.promoter_35_length
                
                iSegment = math.floor((tu_coordinate - 1) / ntPerSegment)
            
                x = segmentLeft + ((tu_coordinate - 1) % ntPerSegment) / ntPerSegment * segmentW
                w = max(1, x1 - min(segmentLeft + segmentW, x1 + tu_length / ntPerSegment * segmentW))
                
                y = chrTop + (iSegment + 1) * segmentHeight + 2
                    
                if tu.name:
                    tip_title = tu.name
                else:
                    tip_title = tu.wid
                tip_title = tip_title.replace("'", "\'")
                    
                promoters.write('<a xlink:href="%s"><rect x="%s" y="%s" width="%s" height="%s" onmousemove="javascript: showToolTip(evt, \'%s\', \'%s\')" onmouseout="javascript: hideToolTip(evt);"/></a>' % (
                    tu.get_absolute_url(species), x, y, w, featureHeight, tip_title, 'Promoter -35 box'))
                    
            if tu.promoter_10_coordinate is not None:
                tu_coordinate = tu.get_coordinate() + tu.promoter_10_coordinate
                tu_length = tu.promoter_10_length
                
                iSegment = math.floor((tu_coordinate - 1) / ntPerSegment)
            
                x = segmentLeft + ((tu_coordinate - 1) % ntPerSegment) / ntPerSegment * segmentW
                w = max(1, x1 - min(segmentLeft + segmentW, x1 + tu_length / ntPerSegment * segmentW))
                
                y = chrTop + (iSegment + 1) * segmentHeight + 2
                
                if tu.name:
                    tip_title = tu.name
                else:
                    tip_title = tu.wid
                tip_title = tip_title.replace("'", "\'")
                    
                promoters.write('<a xlink:href="%s"><rect x="%s" y="%s" width="%s" height="%s" onmousemove="javascript: showToolTip(evt, \'%s\', \'%s\')" onmouseout="javascript: hideToolTip(evt);"/></a>' % (
                    tu.get_absolute_url(species), x, y, w, featureHeight, tip_title, 'Promoter -10 box'))

            for tr in tu.transcriptional_regulations.all():
                if tr.binding_site is not None:    
                    iSegment = math.floor((tr.binding_site.coordinate - 1) / ntPerSegment)
            
                    x = segmentLeft + ((tr.binding_site.coordinate - 1) % ntPerSegment) / ntPerSegment * segmentW
                    w = max(1, x1 - min(segmentLeft + segmentW, x1 + tr.binding_site.length / ntPerSegment * segmentW))
                    
                    y = chrTop + (iSegment + 1) * segmentHeight + 2    
                    
                    if tr.name:
                        tip_title = tr.name
                    else:
                        tip_title = tr.wid
                    tip_title = tip_title.replace("'", "\'")
                        
                    tfSites.write('<a xlink:href="%s"><rect x="%s" y="%s" width="%s" height="%s" onmousemove="javascript: showToolTip(evt, \'%s\', \'%s\')" onmouseout="javascript: hideToolTip(evt);"/></a>' % (
                        tr.get_absolute_url(species), x, y, w, featureHeight, tip_title, 'Transcription factor binding site'))
        
        #features
        featureStyle = '.features rect{fill:#%s;}' % (colors[2], )
        features = StringIO.StringIO()
        
        for feature in self.features.all():        
            iSegment = math.floor((feature.coordinate - 1) / ntPerSegment)
            
            x = segmentLeft + ((feature.coordinate - 1) % ntPerSegment) / ntPerSegment * segmentW
            w = max(1, x1 - min(segmentLeft + segmentW, x1 + feature.length / ntPerSegment * segmentW))
                
            y = chrTop + (iSegment + 1) * segmentHeight + 2
                        
            if feature.type.all().count() > 0:
                type_ = feature.type.all()[0].name
            else:
                type_ = ''
                
            if feature.name:
                tip_title = feature.name
            else:
                tip_title = feature.wid
            tip_title = tip_title.replace("'", "\'")
            
            features.write('<a xlink:href="%s"><rect x="%s" y="%s" width="%s" height="%s" onmousemove="javascript: showToolTip(evt, \'%s\', \'%s\')" onmouseout="javascript: hideToolTip(evt);"/></a>' % (
                feature.get_absolute_url(species), x, y, w, featureHeight, tip_title, type_, ))
    
        return '<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="%s" height="%s" viewport="0 0 %s %s"><style>%s%s%s%s%s%s</style><g class="chr">%s</g><g class="genes">%s</g><g class="promoters">%s</g><g class="tfSites">%s</g><g class="features">%s</g></svg>' % (
            W, H, W, H, style, chrStyle, geneStyle, promoterStyle, tfSiteStyle, featureStyle, chr.getvalue(), genes.getvalue(), promoters.getvalue(), tfSites.getvalue(), features.getvalue())
            
    def get_as_html_structure_local(self, species, is_user_anonymous, start_coordinate = None, end_coordinate = None, highlight_wid = None):
        length = end_coordinate - start_coordinate + 1
        
        W = 636        
        geneHeight = 20
        featureHeight = 10
        H = 2 + geneHeight + 2 + 4 + 1 * (2 + featureHeight) + 2
        
        geneY = 2
        chrY = geneY + geneHeight + 4
        promoterY = chrY + 1 + 2
        featureY = promoterY
        
        #style
        colors = ['3d80b3', '3db34a', 'cd0a0a', 'e78f08']
        style = ''
        for i in range(len(colors)):
            style += '.color-%s{fill:#%s; stroke: #%s;}' % (i, colors[i], colors[i], )
        
        #chromosome
        chrL = 4.5 * len('%s' % start_coordinate) + 4
        chrR = W - 4.5 * len('%s' % end_coordinate) - 2 - 6
        chrW = chrR - chrL        
        
        chrStyle = '\
        .chr text{fill:#222; alignment-baseline:middle; font-size:10px}\
        .chr line{stroke:#666; stroke-width:0.5px;}\
        '
        
        chr = '<text x="%s" y="%s" style="text-anchor:end;">%s</text><line x1="%s" x2="%s" y1="%s" y2="%s" /><text x="%s" y="%s" style="text-anchor:start;">%s</text>' % (
            chrL - 4, chrY, start_coordinate,
            chrL, chrR, chrY, chrY,
            chrR + 2, chrY, end_coordinate)
        
        #genes
        geneStyle = '\
            .genes g polygon{}\
            .genes g text{text-anchor:middle; alignment-baseline:middle; font-size:8px; fill: #222}\
        '

        genes = StringIO.StringIO()      
        genesList = self.genes.filter(coordinate__lte = end_coordinate, coordinate__gte = start_coordinate + 1 - F("length")).prefetch_related('transcription_units', 'transcription_units__transcriptional_regulations').all()

        nTus = 0
        iTUs = {}
        tus = []

        for i, gene in enumerate(genesList):
            if gene.transcription_units.count() > 0:
                tu = gene.transcription_units.all()[:1][0]
                if iTUs.has_key(tu.wid):
                    iTu = iTUs[tu.wid]
                else:
                    tus.append(tu)
                    iTu = nTus
                    iTUs[tu.wid] = iTu
                    nTus += 1
                
                if gene.direction == 'f':
                    x1 = chrL + float(gene.coordinate - start_coordinate) / length * chrW
                    x3 = chrL + float(gene.coordinate + gene.length - 1 - start_coordinate) / length * chrW
                    x2 = max(x1, x3 - 5)
                else:
                    x3 = chrL + float(gene.coordinate - start_coordinate) / length * chrW
                    x1 = chrL + float(gene.coordinate + gene.length - 1 - start_coordinate) / length * chrW
                    x2 = min(x1, x3 + 5)
                    
                x1 = max(chrL, min(chrR, x1))
                x2 = max(chrL, min(chrR, x2))
                x3 = max(chrL, min(chrR, x3))
                            
                y1 = geneY
                y2 = geneY + geneHeight
                
                if highlight_wid is None or gene.wid in highlight_wid:
                    fillOpacity = 0.75
                    strokeOpacity = 1
                    strokeWidth = 3
                else:
                    fillOpacity = 0.15
                    strokeOpacity = 0.35
                    strokeWidth = 1
                
                if math.fabs(x3 - x1) > len(gene.wid) * 5:
                    label = gene.wid
                else:
                    label = ''
                
                if gene.name:
                    tip_title = gene.name
                else:
                    tip_title = gene.wid
                tip_content = 'Transcription unit: %s' % tu.name
                tip_title = tip_title.replace("'", "\'")
                tip_content = tip_content.replace("'", "\'")
                
                gene_abs_url = gene.get_absolute_url(species)
                    
                genes.write('<g style="">\n\
                    <a xlink:href="%s">\n\
                        <polygon class="color-%s" points="%s,%s %s,%s %s,%s %s,%s %s,%s" onmousemove="javascript: showToolTip(evt, \'%s\', \'%s\')" onmouseout="javascript: hideToolTip(evt);" style="fill-opacity: %s; stroke-opacity: %s; stroke-width: %spx"/>\n\
                    </a>\n\
                    <a xlink:href="%s">\n\
                        <text x="%s" y="%s" onmousemove="javascript: showToolTip(evt, \'%s\', \'%s\')" onmouseout="javascript: hideToolTip(evt);">%s</text>\n\
                    </a>\n\
                    </g>\n' % (                
                    gene_abs_url,
                    iTu % len(colors),
                    x1, y1,
                    x2, y1,
                    x3, (y1 + y2) / 2, 
                    x2, y2,
                    x1, y2,
                    tip_title, tip_content,
                    fillOpacity, strokeOpacity, strokeWidth,
                    gene_abs_url,
                    (x1 + x3) / 2, (y1 + y2) / 2 + 1, 
                    tip_title, tip_content,
                    label,
                    ))
                
        #promoters
        promoterStyle = '.promoters rect{fill:#%s; opacity:0.5}' % (colors[0], )
        tfSiteStyle = '.tfSites rect{fill:#%s; opacity:0.5}' % (colors[1], )

        promoters = StringIO.StringIO()
        tfSites = StringIO.StringIO()
        for tu in tus:
            if tu.promoter_35_coordinate is not None:
                tu_coordinate = tu.get_coordinate() + tu.promoter_35_coordinate
                tu_length = tu.promoter_35_length
                
                if not (tu_coordinate > end_coordinate or tu_coordinate + tu_length - 1 < start_coordinate):            
                    x1 = chrL + float(tu_coordinate - start_coordinate) / length * chrW
                    x2 = chrL + float(tu_coordinate + tu_length - 1 - start_coordinate) / length * chrW
                    
                    x1 = max(chrL, min(chrR, x1))
                    x2 = max(chrL, min(chrR, x2))
                        
                    if highlight_wid is None or tu.wid in highlight_wid:
                        opacity = 1
                    else:
                        opacity = 0.25
                        
                    if tu.name:
                        tip_title = tu.name
                    else:
                        tip_title = tu.wid
                    tip_title = tip_title.replace("'", "\'")
                        
                    promoters.write('<a xlink:href="%s"><rect x="%s" y="%s" width="%s" height="%s" onmousemove="javascript: showToolTip(evt, \'%s\', \'%s\')" onmouseout="javascript: hideToolTip(evt);" style="opacity: %s;"/></a>' % (
                        tu.get_absolute_url(species), x1, promoterY, x2 - x1, featureHeight, tip_title, 'Promoter -35 box', opacity))
                    
            if tu.promoter_10_coordinate is not None:
                tu_coordinate = tu.get_coordinate() + tu.promoter_10_coordinate
                tu_length = tu.promoter_10_length
                
                if not (tu_coordinate > end_coordinate or tu_coordinate + tu_length - 1 < start_coordinate):            
                    x1 = chrL + float(tu_coordinate - start_coordinate) / length * chrW
                    x2 = chrL + float(tu_coordinate + tu_length - 1 - start_coordinate) / length * chrW
                    
                    x1 = max(chrL, min(chrR, x1))
                    x2 = max(chrL, min(chrR, x2))
                        
                    if highlight_wid is None or tu.wid in highlight_wid:
                        opacity = 1
                    else:
                        opacity = 0.25
                        
                    if tu.name:
                        tip_title = tu.name
                    else:
                        tip_title = tu.wid
                    tip_title = tip_title.replace("'", "\'")
                    
                    promoters.write('<a xlink:href="%s"><rect x="%s" y="%s" width="%s" height="%s" onmousemove="javascript: showToolTip(evt, \'%s\', \'%s\')" onmouseout="javascript: hideToolTip(evt);" style="opacity: %s;"/></a>' % (
                        tu.get_absolute_url(species), x1, promoterY, x2 - x1, featureHeight, tip_title, 'Promoter -10 box', opacity))

            for tr in tu.transcriptional_regulations.all():
                if tr.binding_site is not None and not (tr.binding_site.coordinate > end_coordinate or tr.binding_site.coordinate + tr.binding_site.length - 1 < start_coordinate):    
                    x1 = chrL + float(tr.binding_site.coordinate - start_coordinate) / length * chrW
                    x2 = chrL + float(tr.binding_site.coordinate + tr.binding_site.length - 1 - start_coordinate) / length * chrW
                    
                    x1 = max(chrL, min(chrR, x1))
                    x2 = max(chrL, min(chrR, x2))
                        
                    if highlight_wid is None or tr.wid in highlight_wid:
                        opacity = 1
                    else:
                        opacity = 0.25
                        
                    if tr.name:
                        tip_title = tr.name
                    else:
                        tip_title = tr.wid
                    tip_title = tip_title.replace("'", "\'")
                    
                    tfSites.write('<a xlink:href="%s"><rect x="%s" y="%s" width="%s" height="%s" onmousemove="javascript: showToolTip(evt, \'%s\', \'%s\')" onmouseout="javascript: hideToolTip(evt);" style="opacity: %s;"/></a>' % (
                        tr.get_absolute_url(species), x1, promoterY, x2 - x1, featureHeight, tip_title, 'Transcription factor binding site', opacity))
                
        #features
        featureStyle = '.features rect{fill:#%s;}' % (colors[2], )
        features = StringIO.StringIO()

        for feature in self.features.all():        
            if feature.coordinate > end_coordinate or feature.coordinate + feature.length - 1 < start_coordinate:
                continue
        
            x1 = chrL + float(feature.coordinate - start_coordinate) / length * chrW
            x2 = chrL + float(feature.coordinate + feature.length - 1 - start_coordinate) / length * chrW
            
            x1 = max(chrL, min(chrR, x1))
            x2 = max(chrL, min(chrR, x2))
                            
            if feature.type.all().count() > 0:
                type_ = feature.type.all()[0].name
            else:
                type_ = ''
                
            if highlight_wid is None or feature.wid in highlight_wid:
                opacity = 1
            else:
                opacity = 0.25
                
            if feature.name:
                tip_title = feature.name
            else:
                tip_title = feature.wid
            tip_title = tip_title.replace("'", "\'")
            
            features.write('<a xlink:href="%s"><rect x="%s" y="%s" width="%s" height="%s" onmousemove="javascript: showToolTip(evt, \'%s\', \'%s\')" onmouseout="javascript: hideToolTip(evt);" style="opacity: %s;"/></a>' % (
                feature.get_absolute_url(species), x1, featureY, x2 - x1, featureHeight, tip_title, type_, opacity))
        
        return '<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="%s" height="%s" viewport="0 0 %s %s"><style>%s%s%s%s%s%s</style><g class="chr">%s</g><g class="genes">%s</g><g class="promoters">%s</g><g class="tfSites">%s</g><g class="features">%s</g></svg>' % (
            W, H, W, H, style, chrStyle, geneStyle, promoterStyle, tfSiteStyle, featureStyle, chr, genes.getvalue(), promoters.getvalue(), tfSites.getvalue(), features.getvalue())
        
    def get_as_html_genes(self, species, is_user_anonymous):
        return ""
        results = []
        for g in self.genes.all():
            results.append('<a href="%s">%s</a>' % (g.get_absolute_url(species), g.wid))
        return format_list_html(results, comma_separated=True)
        
    def get_as_html_features(self, species, is_user_anonymous):
        return ''
        results = []
        for f in self.features.all():
            results.append('<a href="%s">%s</a>' % (f.get_absolute_url(species), f.wid))
        return format_list_html(results, comma_separated=True)
        
    #meta information
    class Meta:
        concrete_entry_model = True
        fieldsets = [
            ('Type', {'fields': ['model_type']}),
            ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), 
            ('Classification', {'fields': ['type']}),
            ('Sequence', {'fields': [
                {'verbose_name': 'Structure', 'name': 'structure'},
                {'verbose_name': 'Extinction coefficient <br/>(260 nm, 25C, pH 7.0)', 'name': 'extinction_coefficient'},
                {'verbose_name': 'pI', 'name': 'pi'},
                ]}), 
            ('Features', {'fields': ['genes', {'verbose_name': 'Other features', 'name': 'features'}]}), 
            ('Function', {'fields': [
                {'verbose_name': 'Reaction participant', 'name':'reaction_stoichiometry_participants'},
                {'verbose_name': 'Complex subunit', 'name':'protein_complex_biosythesis_participants'},
                ]}),
            ('Parameters', {'fields': ['parameters']}),
            ('Comments', {'fields': ['comments', 'publication_references']}),
            ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]}),
            ]
        field_list = [
            'id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'sequence', 'length', 'comments', 'publication_references',  'created_user', 'created_date', 'last_updated_user', 'last_updated_date', 
            ]
        facet_fields = ['type']
        verbose_name='Chromosome'
        verbose_name_plural = 'Chromosomes'
        
        def clean(model, obj_data, all_obj_data=None, all_obj_data_by_model=None):
            if obj_data['sequence'] is not None and obj_data['sequence'] != '' and len(obj_data['sequence']) != obj_data['length']:
                raise ValidationError({'length': 'Length of sequence property must match length property'})

class ChromosomeFeature(SpeciesComponent):
    #parent pointer
    parent_ptr_species_component = OneToOneField(SpeciesComponent, related_name='child_ptr_chromosome_feature', parent_link=True, verbose_name='Species component')
    
    #additional fields
    chromosome = ForeignKey(Chromosome, related_name='features', verbose_name='Chromosome')
    coordinate = PositiveIntegerField(verbose_name='Coordinate (nt)')
    length = PositiveIntegerField(verbose_name='Length (nt)')
    direction = CharField(max_length=10, choices=CHOICES_DIRECTION, verbose_name='Direction')
    
    #getters
    def get_sequence(self):
        seq = self.chromosome.sequence[self.coordinate - 1:self.coordinate - 1 + self.length]
        if self.direction == 'r':
            seq = unicode(Seq(seq, IUPAC.unambiguous_dna).reverse_complement())
        return seq
        
    def get_genes(self):
        genes = []
        for g in self.chromosome.genes.all():
            if (
                g.coordinate <= self.coordinate                   and g.coordinate + g.length - 1 >= self.coordinate
                ) or (
                g.coordinate <= self.coordinate + self.length - 1 and g.coordinate + g.length - 1 >= self.coordinate + self.length - 1
                ) or (
                g.coordinate >= self.coordinate                   and g.coordinate + g.length - 1 <= self.coordinate + self.length - 1
                ):
                genes.append(g)
        return genes
        
    def get_transcription_units(self):
        tus = []
        for tu in self.chromosome.get_transcription_units():
            coordinate = tu.get_coordinate()
            length = tu.get_length()
            if (
                coordinate <= self.coordinate                   and coordinate + length - 1 >= self.coordinate
                ) or (
                coordinate <= self.coordinate + self.length - 1 and coordinate + length - 1 >= self.coordinate + self.length - 1
                ) or (
                coordinate >= self.coordinate                   and coordinate + length - 1 <= self.coordinate + self.length - 1
                ):
                tus.append(tu)
        return tus
        
    #html formatting
    def get_as_html_structure(self, species, is_user_anonymous):
        return self.chromosome.get_as_html_structure(is_user_anonymous, 
            zoom = 1, 
            start_coordinate = self.coordinate - 500, 
            end_coordinate = self.coordinate + self.length + 500, 
            highlight_wid = [self.wid])
            
    def get_as_html_sequence(self, species, is_user_anonymous):
        from cyano.helpers import format_sequence_as_html
        
        direction = CHOICES_DIRECTION[[x[0] for x in CHOICES_DIRECTION].index(self.direction)][1]        
        
        return 'Chromosome: <a href="%s">%s</a>, Coordinate: %s (nt), Length: %s (nt), Direction: %s, Sequence: %s' % (
            self.chromosome.get_absolute_url(species), self.chromosome.wid, 
            self.coordinate, self.length, direction, 
            format_sequence_as_html(self.get_sequence()))
            
    def get_as_html_genes(self, species, is_user_anonymous):
        results = []
        for g in self.get_genes():
            results.append('<a href="%s">%s</a>' % (g.get_absolute_url(species), g.wid))
        return format_list_html(results, comma_separated=True)
    
    def get_as_html_transcription_units(self, species, is_user_anonymous):
        results = []
        for tu in self.get_transcription_units():
            results.append('<a href="%s">%s</a>' % (tu.get_absolute_url(species), tu.wid))
        return format_list_html(results, comma_separated=True)
        
    #meta information
    class Meta:
        concrete_entry_model = True
        fieldsets = [
            ('Type', {'fields': ['model_type']}),
            ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), 
            ('Classification', {'fields': ['type']}),
            ('Structure', {'fields': [
                {'verbose_name': 'Structure', 'name': 'structure'}, 
                {'verbose_name': 'Sequence', 'name': 'sequence'}, 
                {'verbose_name': 'Genes', 'name': 'genes'}, 
                {'verbose_name': 'Transcription units', 'name': 'transcription_units'}, 
                ]}), 
            ('Comments', {'fields': ['comments', 'publication_references']}),
            ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]}),
            ]
        field_list = [
            'id', 'wid', 'name', 'synonyms', 'cross_references', 'type',  'chromosome', 'coordinate', 'length', 'direction', 'comments', 'publication_references',  'created_user', 'created_date', 'last_updated_user', 'last_updated_date', 
            ]
        facet_fields = ['type', 'chromosome', 'direction']
        verbose_name='Chromosome feature'
        verbose_name_plural = 'Chromosome features'
        
        def clean(model, obj_data, all_obj_data=None, all_obj_data_by_model=None):
            if all_obj_data is None:
                chr = Chromosome.objects.get(species__wid=obj_data['species'], wid=obj_data['chromosome'])
            else:
                chr = all_obj_data[obj_data['chromosome']]
                
            if isinstance(chr, Entry):
                chr_len = chr.length
            else:
                chr_len = chr['length']
            
            if obj_data['coordinate'] > chr_len:
                raise ValidationError({'coordinate': 'Coordinate must be less then chromosome length.'})
            if obj_data['length'] > chr_len:
                raise ValidationError({'length': 'Length must be less then chromosome length.'})
                
class Compartment(SpeciesComponent):
    #parent pointer
    parent_ptr_species_component = OneToOneField(SpeciesComponent, related_name='child_ptr_compartment', parent_link=True, verbose_name='Species component')
    
    #additional fields
    def get_protein_complexes(self):
        arr = []
        for obj in ProteinComplex.objects.all():
            if self.pk == obj.get_localization().pk:
                arr.append(obj)
        return arr
        
    #getters
    def get_as_html_biomass_compositions(self, species, is_user_anonymous):
        results = []
        for bm in self.biomass_compositions.all():
            if len(bm.metabolites.all()) == 0:
                continue
            m = bm.metabolites.all()[0]
            results.append('<a href="%s">%s</a>: %.4f' % (m.get_absolute_url(species), m.name, bm.concentration))
        return format_list_html(results, comma_separated=False)
        
    def get_as_html_protein_monomers(self, species, is_user_anonymous):
        results = []
        for p in self.protein_monomers.all():
            results.append('<a href="%s">%s</a>' % (p.get_absolute_url(species), p.wid))
        return format_list_html(results, comma_separated=True)
        
    def get_as_html_protein_complexes(self, species, is_user_anonymous):
        results = []
        for p in self.get_protein_complexes():
            results.append('<a href="%s">%s</a>' % (p.get_absolute_url(species), p.wid))
        return format_list_html(results, comma_separated=True)
    
    #meta information
    class Meta:
        concrete_entry_model = True
        fieldsets = [
            ('Type', {'fields': ['model_type']}),
            ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), 
            ('Classification', {'fields': ['type']}), 
            ('Content', {'fields': [
                {'verbose_name': 'Metabolites (mM)', 'name': 'biomass_compositions'},
                {'verbose_name': 'Protein monomers', 'name': 'protein_monomers'},
                {'verbose_name': 'Protein complexes', 'name': 'protein_complexes'},
                ]}),
            ('Comments', {'fields': ['comments', 'publication_references']}),
            ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]}),
            ]
        field_list = [
            'id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'comments', 'publication_references',  'created_user', 'created_date', 'last_updated_user', 'last_updated_date', 
            ]
        facet_fields = ['type']
        verbose_name='Compartment'
        verbose_name_plural = 'Compartments'
        
class Gene(Molecule):
    #parent pointer
    parent_ptr_molecule = OneToOneField(Molecule, related_name='child_ptr_gene', parent_link=True, verbose_name='Molecule')
    
    #additional fields
    symbol = CharField(max_length=255, blank=True, default='', verbose_name='Symbol')
    chromosome = ForeignKey(Chromosome, related_name='genes', verbose_name='Chromosome')
    coordinate = PositiveIntegerField(verbose_name='Coordinate (nt)')
    length = PositiveIntegerField(verbose_name='Length (nt)')
    direction = CharField(max_length=10, choices=CHOICES_DIRECTION, verbose_name='Direction')
    is_essential = ForeignKey(EntryBooleanData, blank=True, null=True, verbose_name='Is essential', related_name='+')
    expression = ForeignKey(EntryPositiveFloatData, blank=True, null=True, verbose_name='Relative expression', related_name='+')
    half_life = ForeignKey(EntryPositiveFloatData, blank=True, null=True, verbose_name='Half life', related_name='+')
    codons = ManyToManyField(Codon, blank=True, null=True, related_name='genes', verbose_name='Codons')
    amino_acid = ForeignKey('Metabolite', blank=True, null=True, on_delete=SET_NULL, related_name='genes', verbose_name='Amino acid')
    homologs = ManyToManyField(Homolog, blank=True, null=True, related_name='genes', verbose_name='Homologs')    
    
    #getters    
    def get_sequence(self):
        seq = self.chromosome.sequence[self.coordinate - 1:self.coordinate - 1 + self.length]
        if self.direction == 'r':
            seq = unicode(Seq(seq, IUPAC.unambiguous_dna).reverse_complement())
        return seq
        
    def get_length(self):
        return self.length
        
    def get_gc_content(self):
        seq = self.get_sequence()
        if len(seq) == 0:
            return 0
        return float(seq.count('G') + seq.count('C')) / float(len(seq))
        
    def get_empirical_formula(self):
        from cyano.helpers import EmpiricalFormula
        
        seq = self.get_sequence()
        return \
            + Metabolite.objects.get(species__id=self.species.id, wid='AMP').get_empirical_formula() * seq.count('A') \
            + Metabolite.objects.get(species__id=self.species.id, wid='GMP').get_empirical_formula() * seq.count('C') \
            + Metabolite.objects.get(species__id=self.species.id, wid='GMP').get_empirical_formula() * seq.count('G') \
            + Metabolite.objects.get(species__id=self.species.id, wid='UMP').get_empirical_formula() * seq.count('T') \
            - EmpiricalFormula(H=1, O=1) * (len(seq)-1)        

    #http://www.owczarzy.net/extinct.htm
    def get_extinction_coefficient(self):    
        from cyano.helpers import ExtinctionCoefficient
        
        seq = Seq(self.get_sequence(), IUPAC.unambiguous_dna).transcribe()
        
        value = 0;
        for i in range(len(seq) - 1):
            value += ExtinctionCoefficient.pairwise_rna[seq[i]][seq[i+1]]
        for i in range(len(seq)):
            value -= ExtinctionCoefficient.single_rna[seq[i]]
        return value
                
    def get_pi(self):
        return calculate_nucleic_acid_pi(self.get_sequence())

    def get_synthesis_rate(self):
        return math.log(2) / self.half_life.value * self.expression.value
        
    def get_decay_rate(self):        
        return math.log(2) / self.half_life.value
        
    #html formatting    
    def get_as_html_codons(self, species, is_user_anonymous):
        result = []
        for codon in self.codons.all():
            result.append(format_with_evidence(list_item = True, obj = codon, txt = '<tt>%s</tt>' % codon.sequence))        
        return format_list_html(result, force_list=True)
        
    def get_as_html_homologs(self, species, is_user_anonymous):
        results = []
        for h in self.homologs.all():
            results.append(format_with_evidence(list_item = True, obj = h, txt = '%s: <a href="%s">%s</a>' % (h.species, HOMOLOG_SPECIES_URLS[h.species] % h.xid, h.xid)))
        return format_list_html(results, force_list=True)
        
    def get_as_html_structure(self, species, is_user_anonymous):
        return self.chromosome.get_as_html_structure(species, is_user_anonymous, 
            zoom = 1, 
            start_coordinate = self.coordinate - 2500, 
            end_coordinate = self.coordinate + self.length + 2500, 
            highlight_wid = [self.wid])
        
    def get_as_html_sequence(self, species, is_user_anonymous):
        from cyano.helpers import format_sequence_as_html
        
        direction = CHOICES_DIRECTION[[x[0] for x in CHOICES_DIRECTION].index(self.direction)][1]        
        
        return 'Chromosome: <a href="%s">%s</a>, Coordinate: %s (nt), Length: %s (nt), Direction: %s, G/C content: %.1f%%, Sequence: %s' % (
            self.chromosome.get_absolute_url(species), self.chromosome.wid, 
            self.coordinate, self.length, direction, 
            self.get_gc_content() * 100,
            format_sequence_as_html(self.get_sequence()))
        
    #meta information
    class Meta:
        concrete_entry_model = True
        fieldsets = [
            ('Type', {'fields': ['model_type']}),
            ('Name', {'fields': ['wid', 'name', 'symbol', 'synonyms', {'verbose_name': 'Protein product', 'name': 'protein_monomers'}, 'cross_references', 'homologs']}),
            ('Classification', {'fields': ['type']}), 
            ('Structure', {'fields': [
                {'verbose_name': 'Structure', 'name': 'structure'}, 
                {'verbose_name': 'Sequence', 'name': 'sequence'}, 
                {'verbose_name': 'Transcription unit', 'name': 'transcription_units'},
                {'verbose_name': 'Empirical formula (pH 7.5)', 'name': 'empirical_formula'},
                {'verbose_name': 'Molecular weight (pH 7.5; Da)', 'name': 'molecular_weight'},
                ]}), 
            ('Functional genomics', {'fields': [
                'is_essential', 
                'expression', 
                'half_life', 
                'codons', 
                'amino_acid',
                {'verbose_name': 'Extinction coefficient <br/>(260 nm, 25C, pH 7.0)', 'name': 'extinction_coefficient'},
                {'verbose_name': 'pI', 'name': 'pi'},
                ]}), 
            ('Function', {'fields': [
                {'verbose_name': 'Reaction participant', 'name':'reaction_stoichiometry_participants'},
                {'verbose_name': 'Complex subunit', 'name':'protein_complex_biosythesis_participants'},
                ]}),
            ('Parameters', {'fields': ['parameters']}),
            ('Comments', {'fields': ['comments', 'publication_references']}),
            ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]}),
            ]
        field_list = [
            'id', 'wid', 'name', 'symbol', 'synonyms', 'cross_references', 'homologs', 'type', 'chromosome', 'coordinate', 'length', 'direction',  'is_essential', 'expression', 'half_life', 'codons', 'amino_acid', 'comments', 'publication_references', 
            ]
        listing = ['wid', 'symbol']
        facet_fields = ['type', 'chromosome', 'direction', 'is_essential', 'amino_acid']
        verbose_name='Gene'
        verbose_name_plural = 'Genes'
        
        #chromosome coordinate, length
        def clean(model, obj_data, all_obj_data=None, all_obj_data_by_model=None):
            if all_obj_data is None:
                chr = Chromosome.objects.get(species__wid=obj_data['species'], wid=obj_data['chromosome'])
            else:
                chr = all_obj_data[obj_data['chromosome']]

            if isinstance(chr, Entry):
                chr_len = chr.length
            else:
                chr_len = chr['length']
            
            if obj_data['coordinate'] > chr_len:
                raise ValidationError({'coordinate': 'Coordinate must be less then chromosome length.'})
            if obj_data['length'] > chr_len:
                raise ValidationError({'length': 'Length must be less then chromosome length.'})
        
class Metabolite(Molecule):
    #parent pointer
    parent_ptr_molecule = OneToOneField(Molecule, related_name='child_ptr_metabolite', parent_link=True, verbose_name='Molecule')
    
    #additional fields
    traditional_name = TextField(blank=True, default='', verbose_name='Traditional name')
    iupac_name = TextField(blank=True, default='', verbose_name='IUPAC name')
    empirical_formula = TextField(verbose_name='Empirical formula (pH 7.5)', validators=[
        validators.RegexValidator(regex=r'^([A-Z][a-z]*[0-9]*)+$', message='Invalid empirical formula')
        ])
    smiles = TextField(blank=True, default='', verbose_name='SMILES (pH 7.5)')
    charge = IntegerField(verbose_name='Charge (pH 7.5)')
    is_hydrophobic = BooleanField(verbose_name='Is hydrophobic')
    volume = FloatField(null=True, blank=True, verbose_name='van der Waals volume <br/>(pH 7.5; &#8491;<sup>3</sup> molecule<sup>-1</sup>)', validators=[validators.MinValueValidator(0)])
    deltag_formation = FloatField(null=True, blank=True, verbose_name="&Delta;<sub>f</sub>G<sup>'o</sup> (pH 7.5, 25C, I = 0; kJ mol<sup>-1</sup>)")
    pka = FloatField(null=True, blank=True, verbose_name='pK<sub>a</sub>', validators=[validators.MinValueValidator(0)])
    pi = FloatField(null=True, blank=True, verbose_name='pI', validators=[validators.MinValueValidator(0)])
    log_p = FloatField(null=True, blank=True, verbose_name='logP')
    log_d = FloatField(null=True, blank=True, verbose_name='logD (pH 7.5)')
    biomass_composition = ManyToManyField(BiomassComposition, blank=True, null=True, related_name='metabolites', verbose_name='Biomass composition (SP4 media, <br/>5% CO<sub>2</sub>, 37C; mmol gDCW<sup>-1</sup>)')
    media_composition = ForeignKey(MediaComposition, blank=True, null=True, on_delete=SET_NULL, related_name='metabolites', verbose_name='Media composition (SP4; mM)')
    map_coordinates = ManyToManyField(MetaboliteMapCoordinate, blank=True, null=True, related_name='metabolites', verbose_name='Map coordinates')
    
    #getters
    def get_empirical_formula(self):
        from cyano.helpers import EmpiricalFormula
        return EmpiricalFormula(self.empirical_formula)
        
    #calculations
    def calculate_properties(self):        
        subprocess.call('cxcalc name -t preferred name -t traditional logP logD -H 7.5 formalcharge -H 7.5 isoelectricpoint volume "%s"' % (settings.ROOT_DIR, self.smiles, ))
        junk, iupac_name, traditional_name, log_p, log_d, charge, pi, volume = sys.stdout[1]('\t')
        self.iupac_name = iupac_name
        self.traditional_name = traditional_name
        self.log_p = log_p
        self.log_d = log_d
        self.charge = charge
        self.pi = pi
        self.volume = volume
        
        self.deltag_formation = subprocess.call('gcm')
        
    #html formatting
    def get_as_html_structure(self, species, is_user_anonymous):
        from cyano.helpers import draw_molecule
        return draw_molecule(self.smiles, 'svg', 636, 150)
    
    def get_as_html_empirical_formula(self, species, is_user_anonymous):
        from cyano.helpers import EmpiricalFormula
        return EmpiricalFormula(self.empirical_formula).get_as_html()        
        
    def get_as_html_biomass_composition(self, species, is_user_anonymous):
        results = []
        for b in self.biomass_composition.all():
            results.append(format_with_evidence(list_item = True, obj = b, txt = '%s [<a href="%s">%s</a>]' % (b.concentration, b.compartment.get_absolute_url(species), b.compartment.wid)))
        return format_list_html(results, force_list=True)
        
    def get_as_html_media_composition(self, species, is_user_anonymous):
        m = self.media_composition
        if m is None:
            return
        if m.is_diffused:
            txt = '%s (diffused)' % (m.concentration, )
        else:
            txt = m.concentration
            
        return format_with_evidence(obj = m, txt = txt)
            
    def get_as_html_coenzyme_participants(self, species, is_user_anonymous):
        results = []
        for obj in self.coenzyme_participants.all():
            if len(obj.reactions.all()) == 0:
                continue
            rxn = obj.reactions.all()[0]
            results.append('<a href="%s">%s</a><br/>%s' % (rxn.get_absolute_url(species), rxn.name, rxn.get_as_html_stoichiometry(is_user_anonymous)))
        return format_list_html(results, vertical_spacing=True)
        
    def get_as_html_prosthetic_group_participants(self, species, is_user_anonymous):
        results = []
        for obj in self.prosthetic_group_participants.all():
            if len(obj.proteins.all()) == 0:
                continue
            p = obj.proteins.all()[0]
            results.append('<a href="%s">%s</a>' % (p.get_absolute_url(species), p.name))
        return format_list_html(results)
        
    #meta information
    class Meta:
        concrete_entry_model = True
        fieldsets = [
            ('Type', {'fields': ['model_type']}),
            ('Name', {'fields': ['wid', 'name', 'traditional_name', 'iupac_name', 'synonyms', 'cross_references']}),
            ('Classification', {'fields': ['type']}), 
            ('Structure', {'fields': [
                {'verbose_name': 'Structure', 'name': 'structure'}, 
                'empirical_formula', 
                'smiles', 
                'charge', 
                'is_hydrophobic', 
                {'verbose_name': 'Molecular weight (Da)', 'name': 'molecular_weight'}, 
                'volume',
                'deltag_formation', 
                'pka', 
                'pi', 
                'log_p', 
                'log_d']}),
            ('Concentrations', {'fields': ['biomass_composition', 'media_composition']}), 
            ('Function', {'fields': [
                {'verbose_name': 'Coenzyme', 'name':'coenzyme_participants'},                
                {'verbose_name': 'Complex subunit', 'name':'protein_complex_biosythesis_participants'},
                {'verbose_name': 'Prosthetic group', 'name':'prosthetic_group_participants'},
                {'verbose_name': 'Reaction participant', 'name':'reaction_stoichiometry_participants'},
                ]}),
            ('Parameters', {'fields': ['parameters']}),
            ('Comments', {'fields': ['comments', 'publication_references']}),
            ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]}),
            ]
        field_list = [
            'id', 'wid', 'name', 'traditional_name', 'iupac_name', 'synonyms', 'cross_references', 'type',
            'empirical_formula', 'smiles', 'charge', 'is_hydrophobic', 'volume', 'deltag_formation', 'pka', 'pi', 'log_p', 'log_d',
            'biomass_composition', 'media_composition', 
            'map_coordinates',
            'comments',
            'publication_references', 
            'created_user', 'created_date', 'last_updated_user', 'last_updated_date', 
            ]
        facet_fields = ['type', 'charge', 'is_hydrophobic']
        verbose_name='Metabolite'
        verbose_name_plural = 'Metabolites'
        
class Note(SpeciesComponent):
    #parent pointer
    parent_ptr_species_component = OneToOneField(SpeciesComponent, related_name='child_ptr_note', parent_link=True, verbose_name='Species component')
    
    #additional fields
    
    #getters
    
    #meta information
    class Meta:
        concrete_entry_model = True
        fieldsets = [
            ('Type', {'fields': ['model_type']}),
            ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), 
            ('Classification', {'fields': ['type']}), 
            ('Comments', {'fields': ['comments', 'publication_references']}),
            ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]}),
            ]
        field_list = [
            'id', 'wid', 'name', 'synonyms', 'cross_references',
            'type', 
            'comments',
            'publication_references', 
            'created_user', 'created_date', 'last_updated_user', 'last_updated_date', 
            ]
        facet_fields = ['type']
        verbose_name='Note'
        verbose_name_plural = 'Notes'
    
class Parameter(SpeciesComponent):
    #parent pointer
    parent_ptr_species_component = OneToOneField(SpeciesComponent, related_name='child_ptr_parameter', parent_link=True, verbose_name='Species component')
    
    #additional fields
    value = ForeignKey(EntryCharData, verbose_name='Value')
    reactions = ManyToManyField('Reaction', blank=True, null=True, related_name='parameters', verbose_name='Reactions')
    molecules = ManyToManyField('Molecule', blank=True, null=True, related_name='parameters', verbose_name='Molecules')
    state = ForeignKey('State', blank=True, null=True, on_delete=SET_NULL, related_name='parameters', verbose_name='State')
    process = ForeignKey('Process', blank=True, null=True, on_delete=SET_NULL, related_name='parameters', verbose_name='Process')
    
    #getters
    
    #html formatting    
    
    #meta information
    class Meta:
        concrete_entry_model = True
        fieldsets = [
            ('Type', {'fields': ['model_type']}),
            ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), 
            ('Classification', {'fields': ['type']}),
            ('Value', {'fields': ['value']}),
            ('Associations', {'fields': ['reactions', 'molecules', 'state', 'process']}), 
            ('Comments', {'fields': ['comments', 'publication_references']}),
            ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]}),
            ]
        field_list = [
            'id', 'wid', 'name', 'synonyms', 'cross_references',
            'type', 
            'value',
            'reactions', 'molecules', 'state', 'process', 
            'comments',
            'publication_references', 
            'created_user', 'created_date', 'last_updated_user', 'last_updated_date', 
            ]
        facet_fields = ['type', 'reactions', 'molecules', 'state', 'process']
        verbose_name='Misc. parameter'
        verbose_name_plural = 'Misc. parameters'
        
class Pathway(SpeciesComponent):
    #parent pointer
    parent_ptr_species_component = OneToOneField(SpeciesComponent, related_name='child_ptr_pathway', parent_link=True, verbose_name='Species component')
    
    #helpers
    def extract_ecs(self, text):
        """Extracts EC numbers out of a string and returns a list with all numbers"""
        return re.findall(r"[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+", text)
    
    def uniqify(self, seq, idfun=None):
        """Order preserving list uniqifier.
        Source: http://www.peterbe.com/plog/uniqifiers-benchmark
        """
        if idfun is None:
            def idfun(x): return x
        seen = {}
        result = []
        for item in seq:
            marker = idfun(item)
            if marker in seen: continue
            seen[marker] = 1
            result.append(item)
        return result
    
    def get_boehringer_hits(self, species):
        import boehringer.models as bmodels
        species_ecs = CrossReference.objects.filter(species = species, source = "EC").values_list('xid', flat=True)
        enzymes = bmodels.Enzyme.objects.filter(ec__in = species_ecs).order_by('title')

        species_metabolites = Metabolite.objects.filter(species = species).values_list('name', flat=True)
        metabolites = bmodels.Metabolite.objects.filter(title__in = species_metabolites).order_by('title')
        
        return enzymes, metabolites
    
    # TODO: The current logic hardcodes on Boehringer and uses KEGG otherwise
    # Not really nice solution...
    
    def get_as_html_navigator(self, species, is_user_anonymous):
        if self.wid != "Boehringer":
            return ""

        enzymes, metabolites = self.get_boehringer_hits(species)

        template = loader.get_template("cyano/pathway/navigator.html")

        c = Context({'enzymes': enzymes, 'metabolites': metabolites})

        rendered = template.render(c)
        return rendered

    #html formatting
    def get_as_html_reaction_map(self, species, is_user_anonymous):
        W = 731
        H = 600
        
        if self.wid == "Boehringer":
            enzymes, metabolites = self.get_boehringer_hits(species)
            template = loader.get_template("cyano/pathway/reaction_map_boehringer.html")
            c = Context({'enzymes': enzymes, 'metabolites': metabolites,
                         'width': W, 'height': H})
            return template.render(c)
        
        from PIL import Image
        from xml.etree.ElementTree import ElementTree, Element
        
        # All Pathways of the organism
        species_ecs = CrossReference.objects.filter(species = species, source = "EC").values_list('xid', flat=True)

        with open("{}/kegg/img/{}.png".format(settings.STATICFILES_DIRS[0], self.wid), "rb") as image:

            im = Image.open(image)
            iwidth = im.size[0]
            iheight = im.size[1]
    
        with open("{}/kegg/maps/{}.html".format(settings.ROOT_DIR, self.wid)) as infile:
            tree = ElementTree()
            tree.parse(infile)
            root = tree.getroot()
    
            root.tag = "svg"
            root.set("xmlns", "http://www.w3.org/2000/svg")
            root.set("xmlns:xlink", "http://www.w3.org/1999/xlink")
            root.set("width", str(W))
            root.set("height", str(min(iheight,H)))
            root.set("viewport", "0 0 {} {}".format(str(W), str(min(iheight,H))))

            script = Element("script")
            script.set("xlink:href", "{}kegg/js/SVGPan.js".format(settings.STATIC_URL))
            root.append(script)
            
            graphics = Element("g")
            graphics.set("id", "viewport")
            root.append(graphics)
            
            image = Element("image")
            image.set("x", "0")
            image.set("y", "0")
            image.set("width", str(iwidth))
            image.set("height", str(iheight))
            image.set("xlink:href", "{}kegg/img/{}.png".format(settings.STATIC_URL, self.wid))
            graphics.append(image)
    
            areas = tree.findall("area")
            
            # Draw polys before the rest because they are line segments and sometimes
            # overlap other elements
            pending = []
            
            for area in areas:
                shape = area.get("shape")
                coords = area.get("coords")
                url = area.get("href")
                title = area.get("title")
    
                if shape is None:
                    continue
                
                if coords is None:
                    continue
                
                if url is None:
                    continue
                
                if shape == "poly":
                    shape = "polygon"
                
                url = "http://www.kegg.jp" + url
            
                coords = coords.split(",")
            
                area.attrib.pop("shape")
                area.attrib.pop("coords")
                area.attrib.pop("href")
                area.attrib.pop("title")
            
                area.tag = shape
                
                elem = Element("a")
                
                # Pathways are blue
                # Something with EC numbers green
                # Everything else red 
                color_component = "rgb(255,0,0)"
                fill_opacity = "0.0"
                fill_color = "green" # not displayed with 0.0
                
                # Found objects in the organism show a background color
                
                # Repoint pathway links to our versions
                if "show_pathway?" in url:
                    pathway_name = url[url.index("?") + 1:]
                    color_component = "rgb(0,0,255)"
                    
                    try:
                        pw_obj = Pathway.objects.get(wid = pathway_name, species = species)
                        elem.set("xlink:href", pw_obj.get_absolute_url(species))
                        fill_opacity = "0.2"
                        fill_color = "blue"
                    except ObjectDoesNotExist:
                        elem.set("xlink:href", url)
                        elem.set("target", "_blank")
                        pass
                else:
                    elem.set("xlink:href", url)
                    elem.set("target", "_blank")
    
                elem.set("xlink:title", title)
                root.remove(area)
                elem.append(area)
                
                ecs = self.uniqify(self.extract_ecs(title))
                
                if len(ecs) > 0:
                    color_component = "rgb(0,255,0)"
                    for ec in ecs:
                        if ec in species_ecs:
                            fill_opacity = "0.2"

                if shape == "circle":
                    pending.append(elem)
    
                    area.set("cx", coords[0])
                    area.set("cy", coords[1])
                    area.set("r", str(int(coords[2]) + 1))
                elif shape == "rect":
                    pending.append(elem)
    
                    area.set("x", coords[0])
                    area.set("y", coords[1])
                    area.set("width", str(int(coords[2]) - int(coords[0])))
                    area.set("height", str(int(coords[3]) - int(coords[1])))
                elif shape == "polygon":
                    graphics.append(elem)
                       
                    points = zip(*2*[iter(coords)])
    
                    area.set("points", " ".join([",".join(x) for x in points]))

                area.set("style", "stroke-width:1;stroke:{};fill-opacity:{};fill:{}".format(color_component, fill_opacity, fill_color))
            
            for elem in pending:
                graphics.append(elem)
            
            out = StringIO.StringIO()
            out.write('<script type="text/javascript" src="' + settings.STATIC_URL + 'kegg/js/jquery-svgpan.js"></script>')
            tree.write(out)
            out.write('<script type="text/javascript">$("svg").svgPan("viewport");</script>')

            return out.getvalue()

    #meta information
    class Meta:
        concrete_entry_model = True
        fieldsets = [
            ('Type', {'fields': ['model_type']}),
            ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), 
            ('Classification', {'fields': ['type']}),
            ('Navigator', {'fields' : [
                {'verbose_name': 'Navigate to', 'name': 'navigator'}
                ]}), 
            ('Reactions', {'fields': [
                {'verbose_name': 'Reactions', 'name': 'reaction_map'},
                'reactions',
                ]}),
            ('Comments', {'fields': ['comments', 'publication_references']}),
            ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]}),
            ]            
        field_list = [
            'id', 'wid', 'name', 'synonyms', 'cross_references',
            'type', 
            'comments',
            'publication_references', 
            'created_user', 'created_date', 'last_updated_user', 'last_updated_date', 
            ]
        facet_fields = ['type']
        verbose_name='Pathway'
        verbose_name_plural = 'Pathways'
        
class Process(SpeciesComponent):
    #parent pointer
    parent_ptr_species_component = OneToOneField(SpeciesComponent, related_name='child_ptr_process', parent_link=True, verbose_name='Species component')
    
    #additional fields
    initialization_order = PositiveIntegerField(verbose_name='Initialization order')
    
    #getters
    
    #html formatting    
    def get_as_html_reactions(self, species, is_user_anonymous):
        results = []
        for reaction in self.reactions.all():
            results.append('<a href="%s">%s</a><br/>%s' % (reaction.get_absolute_url(species), reaction.name, reaction.get_as_html_stoichiometry(is_user_anonymous)))
        return format_list_html(results, vertical_spacing=True)
        
    def get_as_html_formed_complexes(self, species, is_user_anonymous):
        results = []
        for complex in self.formed_complexes.all():
            results.append('<a href="%s">%s</a><br/>%s' % (complex.get_absolute_url(species), complex.name, complex.get_as_html_biosynthesis(is_user_anonymous)))
        return format_list_html(results, vertical_spacing=True)
    
    #meta information
    class Meta:        
        concrete_entry_model = True
        fieldsets = [
            ('Type', {'fields': ['model_type']}),
            ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), 
            ('Classification', {'fields': ['type']}),
            ('Implementation', {'fields': ['initialization_order']}), 
            ('Reactions', {'fields': [{'verbose_name': 'Chemical reactions', 'name': 'reactions'}, {'verbose_name': 'Complex formation reactions', 'name': 'formed_complexes'}]}),
            ('Parameters', {'fields': ['parameters']}),
            ('Comments', {'fields': ['comments', 'publication_references']}),
            ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]}),
            ]
        field_list = [
            'id', 'wid', 'name', 'synonyms', 'cross_references',
            'type', 
            'initialization_order', 
            'comments',
            'publication_references', 
            'created_user', 'created_date', 'last_updated_user', 'last_updated_date', 
            ]
        facet_fields = ['type']
        verbose_name='Process'
        verbose_name_plural = 'Processes'
        
        def validate_unique(model, model_objects_data, all_obj_data=None, all_obj_data_by_model=None):
            order = []
            for obj in model_objects_data:
                if isinstance(obj, Entry):
                    order.append(obj.initialization_order)
                else:
                    order.append(obj['initialization_order'])        
                
            if len(set(order)) < len(order):
                raise ValidationError({'initialization_order': 'Initialization order must unique within a species'})
        
class ProteinComplex(Protein):
    #parent pointer
    parent_ptr_protein = OneToOneField(Protein, related_name='child_ptr_protein_complex', parent_link=True, verbose_name='Protein')
    
    #additional fields
    biosynthesis = ManyToManyField(ProteinComplexBiosythesisParticipant, related_name='protein_complexes', verbose_name='Biosynthesis')
    disulfide_bonds = ManyToManyField(DisulfideBond, blank=True, null=True, related_name='protein_complexes', verbose_name='Disulfide bonds (pH 7.5)')    
    formation_process = ForeignKey('Process', blank=True, null=True, on_delete=SET_NULL, related_name='formed_complexes', verbose_name='Formation process')
    
    #getters
    def get_num_subunits(self):
        n = 0
        for subunit in self.biosynthesis.all():
            if subunit.coefficient < 0:
                if subunit.molecule.model_type in ['TranscriptionUnit', 'Gene', 'ProteinMonomer']:
                    n -= subunit.coefficient
                elif subunit.molecule.model_type is 'ProteinComplex':
                    n -= subunit.coefficient * getEntry(species_wid = self.species.wid, wid = subunit.molecule.wid).get_num_subunits()
        return n

    def get_as_html_num_subunits(self, species, is_user_anonymous):
        return '%0.0f' % self.get_num_subunits()
        
    def get_empirical_formula(self):
        from cyano.helpers import EmpiricalFormula, getModel
        
        formula = EmpiricalFormula()
        for participant in self.biosynthesis.all():
            if not (participant.molecule.wid == self.wid and participant.coefficient == 1):
                molecule = participant.molecule
                molecule = getModel(molecule.model_type).objects.get(id=molecule.id)
                formula += molecule.get_empirical_formula() * -participant.coefficient
            
        return formula
    
    def get_localization(self):
        from cyano.helpers import getModel
        
        localizations = []
        for participant in self.biosynthesis.all():
            if participant.coefficient > 0:
                continue
                
            molecule = participant.molecule
            molecule = getModel(molecule.model_type).objects.get(id=molecule.id)
            if isinstance(molecule, ProteinMonomer):
                localizations.append(molecule.localization)
            elif isinstance(molecule, ProteinComplex):
                localizations.append(molecule.get_localization())
                
        localizations = list(set(localizations))
        if len(localizations) == 1:
            return localizations[0]
        
        if len(localizations) == 2 and len(set(['c', 'm']) & set([x.wid for x in localizations])) == 2:
            return Compartment.objects.get(species__id=self.species.id, wid='m')
        
        if len(localizations) == 3 and len(set(['c', 'm', 'e']) & set([x.wid for x in localizations])) == 3:
            return Compartment.objects.get(species__id=self.species.id, wid='m')
    
        raise TypeError(str(localizations))
        
    def get_half_life(self):
        from cyano.helpers import getModel
        
        val = 0
        for participant in self.biosynthesis.all():
            if participant.coefficient == 1 and participant.molecule.wid == self.wid:
                continue
                
            molecule = participant.molecule
            molecule = getModel(molecule.model_type).objects.get(id=molecule.id)
            
            if isinstance(molecule, Protein):
                val += -participant.coefficient * molecule.get_molecular_weight() * molecule.get_half_life()
            elif isinstance(molecule, Gene):
                val += -participant.coefficient * molecule.get_molecular_weight() * molecule.half_life.value
                
        return val / self.get_molecular_weight()
        
    def get_neg_aa(self):
        from cyano.helpers import getModel
        
        val = 0
        for participant in self.biosynthesis.all():
            if participant.coefficient == 1 and participant.molecule.wid == self.wid:
                continue
                
            molecule = participant.molecule
            molecule = getModel(molecule.model_type).objects.get(id=molecule.id)
            if isinstance(molecule, Protein):
                val += participant.coefficient * molecule.get_neg_aa()
                
        return val
        
    def get_pos_aa(self):
        from cyano.helpers import getModel
        
        val = 0
        for participant in self.biosynthesis.all():
            if participant.coefficient == 1 and participant.molecule.wid == self.wid:
                continue
                
            molecule = participant.molecule
            molecule = getModel(molecule.model_type).objects.get(id=molecule.id)
            if isinstance(molecule, Protein):
                val += participant.coefficient * molecule.get_neg_aa()
                
        return val
        
    def get_extinction_coefficient(self):
        from cyano.helpers import getModel
        
        val = 0
        for participant in self.biosynthesis.all():
            if participant.coefficient == 1 and participant.molecule.wid == self.wid:
                continue
                
            molecule = participant.molecule
            molecule = getModel(molecule.model_type).objects.get(id=molecule.id)
            if isinstance(molecule, (Protein, Gene, )):
                val += participant.coefficient * molecule.get_extinction_coefficient()
                
        return val
        
    #html formatting    
    def get_as_html_biosynthesis(self, species, is_user_anonymous):
        compartments = []
        for s in self.biosynthesis.all():
            compartments.append(s.compartment)
        compartments = list(set(compartments))
            
        pos = []
        neg = []
        for s in self.biosynthesis.all():
            if s.coefficient < 0:
                tmp = ''
                if s.coefficient != -1:
                    tmp += '(%d) ' % -s.coefficient
                tmp += '<a href="%s">%s</a>' % (s.molecule.get_absolute_url(species), s.molecule.wid)
                if len(compartments) > 1:
                    tmp += '[<a href="%s">%s</a>]' % (s.compartment.get_absolute_url(species), s.compartment.wid)
                pos.append(tmp)
            else:
                tmp = ''
                if s.coefficient != 1:
                    tmp += '(%d) ' % s.coefficient
                tmp += '<a href="%s">%s</a>' % (s.molecule.get_absolute_url(species), s.molecule.wid)
                if len(compartments) > 1:
                    tmp += '[<a href="%s">%s</a>]' % (s.compartment.get_absolute_url(species), s.compartment.wid)
                neg.append(tmp)
                    
        result = ''
        if len(compartments) == 1:
            result += '[<a href="%s">%s</a>]: ' % (compartments[0].get_absolute_url(species), compartments[0].wid)
        result += ' + '.join(pos)
        result += ' &rArr; '
        result += ' + '.join(neg)
        return format_with_evidence(obj = self.biosynthesis.all(), txt = result)
        
    def get_as_html_disulfide_bonds(self, species, is_user_anonymous):
        results = [];
        for b in self.disulfide_bonds.all():
            results.append(format_with_evidence(list_item = True, obj = b, txt = '<a href="%s">%s</a>: %s-%s' % (b.protein_monomer.get_absolute_url(species), b.protein_monomer.wid, b.residue_1, b.residue_2)))
        return format_list_html(results, force_list=True)
        
    def get_as_html_localization(self, species, is_user_anonymous):
        localization = self.get_localization()
        return '<a href="%s">%s</a>' % (localization.get_absolute_url(species), localization.wid, )  
                    
    #meta information
    class Meta:        
        concrete_entry_model = True
        fieldsets = [
            ('Type', {'fields': ['model_type']}),
            ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), 
            ('Classification', {'fields': ['type']}),
            ('Structure', {'fields': [
                'biosynthesis', 
                {'verbose_name': 'No. subunits', 'name': 'num_subunits'},
                'disulfide_bonds', 
                'prosthetic_groups', 
                'dna_footprint',
                {'verbose_name': 'Empirical formula (pH 7.5)', 'name': 'empirical_formula'},
                {'verbose_name': 'Molecular weight (pH 7.5; Da)', 'name': 'molecular_weight'},
                {'verbose_name': 'Extinction coefficient <br/>(260 nm, 25C, pH 7.0)', 'name': 'extinction_coefficient'},
                {'verbose_name': 'Half life (OD (600 nm) = 0.3, <br/>M9 media, 36C; min)', 'name': 'half_life'},
                ]}),
            ('Synthesis', {'fields': [
                'formation_process', 
                'chaperones', 
                {'verbose_name': 'Localization', 'name': 'localization'},
                ]}), 
            ('Regulation', {'fields': ['regulatory_rule']}),
            ('Function', {'fields': [
                {'verbose_name': 'Enzyme', 'name': 'enzyme_participants'},
                {'verbose_name': 'Transcriptional regulation', 'name': 'transcriptional_regulations'},
                {'verbose_name': 'Protein folding substrates', 'name': 'chaperone_substrates'},
                {'verbose_name': 'Reaction participant', 'name':'reaction_stoichiometry_participants'},
                {'verbose_name': 'Complex subunit', 'name':'protein_complex_biosythesis_participants'},
                ]}), 
            ('Parameters', {'fields': ['parameters']}),
            ('Comments', {'fields': ['comments', 'publication_references']}),
            ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]}),
            ]            
        field_list = [
            'id', 'wid', 'name', 'synonyms', 'cross_references',
            'type', 
            'biosynthesis', 'disulfide_bonds', 'prosthetic_groups', 'dna_footprint',
            'formation_process', 'chaperones', 
            'regulatory_rule',
            'comments',
            'publication_references', 
            'created_user', 'created_date', 'last_updated_user', 'last_updated_date', 
            ]
        facet_fields = ['type', 'dna_footprint__binding', 'dna_footprint__region', 'formation_process', 'chaperones']
        verbose_name='Protein complex'
        verbose_name_plural = 'Protein complexes'
        
        def clean(model, obj_data, all_obj_data=None, all_obj_data_by_model=None):
            from cyano.helpers import getModel, getEntry

            #biosynthesis
            coeff = 0
            for b in obj_data['biosynthesis']:
                if b['molecule'] == obj_data['wid']:
                    coeff += b['coefficient']
            
            if coeff != 1:
                raise ValidationError({'biosynthesis': 'Protein complex must appear on the right side of the biosynthesis reaction'})
                
            #disulfide bonds
            for dsfb in obj_data['disulfide_bonds']:
                mon_wid = dsfb['protein_monomer']
        
                if all_obj_data is None:
                    mon = ProteinMonomer.objects.get(species__wid=obj_data['species'], wid=mon_wid)
                else:
                    mon = all_obj_data[mon_wid]
                if isinstance(mon, Entry):
                    gene_wid = mon.gene.wid
                else:
                    gene_wid = mon['gene']
                
                if all_obj_data is None:
                    gene = Gene.objects.get(species_wid=obj_data['species'], wid=gene_wid)
                else:
                    gene = all_obj_data[gene_wid]
                if isinstance(gene, Entry):
                    mon_len = gene.length / 3
                else:
                    mon_len = gene['length'] / 3
                
                if dsfb['residue_1'] > mon_len:
                    raise ValidationError({'disulfide_bond': 'Residue-1 must be less then protein length'})
                if dsfb['residue_2'] > mon_len:
                    raise ValidationError({'disulfide_bond': 'Residue-2 must be less then protein length'})                
            
            #biosynthesis residues
            for b in obj_data['biosynthesis']:
                if all_obj_data is None:
                    molecule = getEntry(species_wid=obj_data['species'], wid=b['molecule'])
                else:
                    molecule = all_obj_data[b['molecule']]
                if isinstance(molecule, Entry):
                    molecule_type = molecule.model_type
                else:
                    molecule_type = molecule['model_type']
                molecule_len = None
                if molecule_type == 'Gene':
                    if isinstance(molecule, Entry):
                        molecule_len = molecule.length
                    else:
                        molecule_len = molecule['length']
                elif molecule_type == 'ProteinMonomer':
                    if isinstance(molecule, Entry):
                        gene_wid = molecule.gene.wid
                    else:
                        gene_wid = molecule['gene']
                    if all_obj_data is None:
                        gene = Gene.objects.get(species__wid=obj_data['species'], wid=gene_wid)
                    else:
                        gene = all_obj_data[gene_wid]
                    if isinstance(gene, Entry):
                        molecule_len = gene.length
                    else:
                        molecule_len = gene['length']
                
                if b['residue'] is not None and not issubclass(getModel(molecule_type), (ProteinMonomer, Gene, )):
                    raise ValidationError({'biosynthesis': 'Residue must be null'})
                
                if b['residue'] is not None and b['residue'] > molecule_len:
                    raise ValidationError({'biosynthesis': 'Residue must be less than molecule length'})
        
class ProteinMonomer(Protein):
    #parent pointer
    parent_ptr_protein = OneToOneField(Protein, related_name='child_ptr_protein_monomer', parent_link=True, verbose_name='Protein')
    
    #additional fields
    gene = ForeignKey(Gene, related_name='protein_monomers', verbose_name='Gene')    
    is_n_terminal_methionine_cleaved = ForeignKey(EntryBooleanData, null = True, verbose_name='Is N-terminal methionine cleaved', related_name='+')
    localization = ForeignKey(Compartment, null = True, related_name='protein_monomers', verbose_name='Localization')
    signal_sequence = ForeignKey(SignalSequence, blank=True, null=True, related_name='protein_monomers', on_delete=SET_NULL, verbose_name='Sequence sequence')

    #getters
    def get_sequence(self, species):
        return unicode(Seq(self.gene.get_sequence(), IUPAC.unambiguous_dna).translate(table=species.genetic_code))
        
    def get_length(self):
        return len(self.get_sequence())
        
    def get_neg_aa(self):
        seq = self.get_sequence()
        return seq.count('E') + seq.count('D')
        
    def get_pos_aa(self):
        seq = self.get_sequence()
        return seq.count('R') + seq.count('H') + seq.count('K')
        
    def get_n_terminal_aa(self):
        return self.get_sequence()[0]
        
    def get_empirical_formula(self):
        from cyano.helpers import EmpiricalFormula
        
        seq = self.get_sequence()
        return \
            + Metabolite.objects.get(species__id=self.species.id, wid='ALA').get_empirical_formula() * seq.count('A') \
            + Metabolite.objects.get(species__id=self.species.id, wid='ARG').get_empirical_formula() * seq.count('R') \
            + Metabolite.objects.get(species__id=self.species.id, wid='ASN').get_empirical_formula() * seq.count('N') \
            + Metabolite.objects.get(species__id=self.species.id, wid='ASP').get_empirical_formula() * seq.count('D') \
            + Metabolite.objects.get(species__id=self.species.id, wid='CYS').get_empirical_formula() * seq.count('C') \
            + Metabolite.objects.get(species__id=self.species.id, wid='GLU').get_empirical_formula() * seq.count('E') \
            + Metabolite.objects.get(species__id=self.species.id, wid='GLN').get_empirical_formula() * seq.count('Q') \
            + Metabolite.objects.get(species__id=self.species.id, wid='GLY').get_empirical_formula() * seq.count('G') \
            + Metabolite.objects.get(species__id=self.species.id, wid='HIS').get_empirical_formula() * seq.count('H') \
            + Metabolite.objects.get(species__id=self.species.id, wid='ILE').get_empirical_formula() * seq.count('I') \
            + Metabolite.objects.get(species__id=self.species.id, wid='LEU').get_empirical_formula() * seq.count('L') \
            + Metabolite.objects.get(species__id=self.species.id, wid='LYS').get_empirical_formula() * seq.count('K') \
            + Metabolite.objects.get(species__id=self.species.id, wid='MET').get_empirical_formula() * (seq.count('M')-1) + Metabolite.objects.get(species__id=self.species.id, wid='FMET').get_empirical_formula() * (1) \
            + Metabolite.objects.get(species__id=self.species.id, wid='PHE').get_empirical_formula() * seq.count('F') \
            + Metabolite.objects.get(species__id=self.species.id, wid='PRO').get_empirical_formula() * seq.count('P') \
            + Metabolite.objects.get(species__id=self.species.id, wid='SER').get_empirical_formula() * seq.count('S') \
            + Metabolite.objects.get(species__id=self.species.id, wid='THR').get_empirical_formula() * seq.count('T') \
            + Metabolite.objects.get(species__id=self.species.id, wid='TRP').get_empirical_formula() * seq.count('W') \
            + Metabolite.objects.get(species__id=self.species.id, wid='TYR').get_empirical_formula() * seq.count('Y') \
            + Metabolite.objects.get(species__id=self.species.id, wid='VAL').get_empirical_formula() * seq.count('V') \
            - EmpiricalFormula(H=2, O=1) * (len(seq)-1)    
        
    def get_pi(self):
        seq = self.get_sequence()
        
        numAsp = float(seq.count('D'))
        numGlu = float(seq.count('E'))
        numCys = float(seq.count('C'))
        numTyr = float(seq.count('Y'))
        numHis = float(seq.count('H'))
        numLys = float(seq.count('K'))
        numArg = float(seq.count('R'))

        pH = 6.5             #starting point pI = 6.5 - theoretically it should be 7, but average protein pI is 6.5 so we increase the probability
        pHprev = 0.0         #of finding the solution
        pHnext = 14.0        #0-14 is possible pH range
        E = 0.01             #epsilon means precision [pI = pH \pm E]

        #the infinite loop
        while True:
            # we are using pK values form Wikipedia as they give quite good approximation
            # if you want you can change it
            QN1 = -    1. / (1. + 10.**( 3.65 - pH)) #C-terminal charge
            QN2 = -numAsp / (1. + 10.**( 3.90 - pH)) #D charge
            QN3 = -numGlu / (1. + 10.**( 4.07 - pH)) #E charge
            QN4 = -numCys / (1. + 10.**( 8.18 - pH)) #C charge
            QN5 = -numTyr / (1. + 10.**(10.46 - pH)) #Y charge
            QP1 =  numHis / (1. + 10.**(pH - 6.04))  #H charge
            QP2 =      1. / (1. + 10.**(pH - 8.20))  #NH2 charge
            QP3 =  numLys / (1. + 10.**(pH -10.54))  #K charge
            QP4 =  numArg / (1. + 10.**(pH -12.48))  #R charge

            NQ = QN1 + QN2 + QN3 + QN4 + QN5 + QP1 + QP2 + QP3 + QP4 #net charge in given pH

            if pH >= 14.0:
                raise

            #%%%%%%%%%%%%%%%%%%%%%%%%%   BISECTION   %%%%%%%%%%%%%%%%%%%%%%%%

            #we are out of range, thus the new pH value must be smaller
            if NQ < 0.:
                temp = pH
                pH = pH - ((pH - pHprev) / 2.)
                pHnext = temp

            #we used to small pH value, so we have to increase it
            else:
                temp = pH
                pH = pH + ((pHnext - pH ) / 2.)
                pHprev = temp

            #terminal condition, finding isoelectric point with given precision
            if (pH - pHprev < E) and (pHnext - pH < E):
                break

        value = pH            
        
    def get_half_life(self):    
        return 20 * 60
        
    #http://ca.expasy.org/tools/protparam-doc.html
    def get_instability(self, species):
        from cyano.helpers import DipeptideInstabilityWeight
        
        seq = self.get_sequence(species)
        value = 0.;
        for i in range(len(seq)-1):
            if seq[i] != '*' and seq[i+1] != '*':
                value += DipeptideInstabilityWeight.value[seq[i]][seq[i+1]]
        return 10. / float(len(seq)) * value;
        
    #http://ca.expasy.org/tools/protparam-doc.html
    def get_is_stable(self, species):
        return self.get_instability(species) < 40.
        
    #http://ca.expasy.org/tools/protparam-doc.html
    def get_aliphatic(self, species):
        seq = self.get_sequence(species)
        return 100. * ( \
            + 1.0 * float(seq.count('A')) \
            + 2.9 * float(seq.count('V')) \
            + 3.9 * float(seq.count('I')) \
            + 3.9 * float(seq.count('L')) \
            ) / float(len(seq))
        
    #http://ca.expasy.org/tools/protparam-doc.html
    def get_gravy(self, species):        
        seq = self.get_sequence(species)
        return \
            ( \
            + 1.8 * float(seq.count('A')) \
            - 4.5 * float(seq.count('R')) \
            - 3.5 * float(seq.count('N')) \
            - 3.5 * float(seq.count('D')) \
            + 2.5 * float(seq.count('C')) \
            - 3.5 * float(seq.count('Q')) \
            - 3.5 * float(seq.count('E')) \
            - 0.4 * float(seq.count('G')) \
            - 3.2 * float(seq.count('H')) \
            + 4.5 * float(seq.count('I')) \
            + 3.8 * float(seq.count('L')) \
            - 3.9 * float(seq.count('K')) \
            + 1.9 * float(seq.count('M')) \
            + 2.8 * float(seq.count('F')) \
            - 1.6 * float(seq.count('P')) \
            - 0.8 * float(seq.count('S')) \
            - 0.7 * float(seq.count('T')) \
            - 0.9 * float(seq.count('W')) \
            - 1.3 * float(seq.count('Y')) \
            + 4.2 * float(seq.count('V')) \
             ) / float(len(seq))
        
    #Source: http://ca.expasy.org/tools/protparam-doc.html
    def get_extinction_coefficient(self):
        seq = self.get_sequence()
        return \
            + seq.count('W') * 5500 \
            + seq.count('Y') * 1490 \
            + seq.count('C') * 125    
        
    #html formatting
    def get_as_html_sequence(self, species, is_user_anonymous):
        from cyano.helpers import format_sequence_as_html
        return format_sequence_as_html(self.get_sequence(species))
        
    def get_as_html_signal_sequence(self, species, is_user_anonymous):
        ss = self.signal_sequence
        if ss is None:
            return
        return format_with_evidence(obj = ss, txt = 'Type: %s, Location: %s, Length: %s (nt)' % (ss.type, ss.location, ss.length))
        
    def get_as_html_disulfide_bonds(self, species, is_user_anonymous):
        results = [];
        for b in self.disulfide_bonds.all():
            results.append('<a href="%s">%s</a>: %s-%s' % (b.protein_complexes.all()[0].get_absolute_url(species), b.protein_complexes.all()[0].wid, b.residue_1, b.residue_2))
        return format_list_html(results)
        
    def get_as_html_instability(self, species, is_user_anonymous):
        return self.get_instability(species)
        
    def get_as_html_is_stable(self, species, is_user_anonymous):
        return self.get_is_stable(species)
    
    def get_as_html_aliphatic(self, species, is_user_anonymous):    
        return self.get_aliphatic(species)
    
    def get_as_html_gravy(self, species, is_user_anonymous):
        return self.get_gravy(species)        

    #meta information
    class Meta:    
        concrete_entry_model = True
        fieldsets = [
            ('Type', {'fields': ['model_type']}),
            ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), 
            ('Classification', {'fields': ['type']}),
            ('Genetics', {'fields': ['gene']}), 
            ('Structure', {'fields': [
                {'verbose_name': 'Sequence', 'name': 'sequence'}, 
                'is_n_terminal_methionine_cleaved', 
                'signal_sequence', 
                'prosthetic_groups', 
                'disulfide_bonds', 
                'dna_footprint',
                {'verbose_name': 'Empirical formula (pH 7.5)', 'name': 'empirical_formula'},
                {'verbose_name': 'Molecular weight (pH 7.5; Da)', 'name': 'molecular_weight'},
                {'verbose_name': 'Extinction coefficient <br/>(260 nm, 25C, pH 7.0)', 'name': 'extinction_coefficient'},
                {'verbose_name': 'pI', 'name': 'pi'},
                {'verbose_name': 'Instability index', 'name': 'instability'},
                {'verbose_name': 'Is stable', 'name': 'is_stable'},
                {'verbose_name': 'Aliphatic index', 'name': 'aliphatic'},
                {'verbose_name': 'GRAVY (25C, pH 7.0)', 'name': 'gravy'},
                {'verbose_name': 'Half life (OD (600 nm) = 0.3, <br/>M9 media, 36C; min)', 'name': 'half_life'},
                ]}), 
            ('Synthesis', {'fields': ['localization', 'chaperones',]}), 
            ('Regulation', {'fields': ['regulatory_rule']}), 
            ('Function', {'fields': [
                {'verbose_name': 'Enzyme', 'name': 'enzyme_participants'},
                {'verbose_name': 'Transcriptional regulation', 'name': 'transcriptional_regulations'},
                {'verbose_name': 'Protein folding substrates', 'name': 'chaperone_substrates'},    
                {'verbose_name': 'Reaction participant', 'name':'reaction_stoichiometry_participants'},
                {'verbose_name': 'Complex subunit', 'name':'protein_complex_biosythesis_participants'},
                ]}), 
            ('Parameters', {'fields': ['parameters']}),
            ('Comments', {'fields': ['comments', 'publication_references']}),
            ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]}),
            ]        
        field_list = [
            'id', 'wid', 'name', 'synonyms', 'cross_references',
            'type', 
            'gene', 
            'is_n_terminal_methionine_cleaved', 'signal_sequence', 'prosthetic_groups', 'dna_footprint', 
            'localization', 'chaperones',
            'regulatory_rule', 
            'comments',
            'publication_references', 
            'created_user', 'created_date', 'last_updated_user', 'last_updated_date', 
            ]
        facet_fields = ['type', 'is_n_terminal_methionine_cleaved__value', 'signal_sequence__type', 'signal_sequence__location', 'dna_footprint__binding', 'dna_footprint__region', 'localization', 'chaperones']
        verbose_name='Protein monomer'
        verbose_name_plural = 'Protein monomers'
        
        def clean(model, obj_data, all_obj_data=None, all_obj_data_by_model=None):
            if obj_data['signal_sequence'] is not None:
                if all_obj_data is None:
                    gene = Gene.objects.get(species__wid=obj_data['species'], wid=obj_data['gene'])
                else:
                    gene = all_obj_data[obj_data['gene']]
                if isinstance(gene, Entry):
                    mon_len = gene.length / 3
                else:
                    mon_len = gene['length'] / 3
                
                if obj_data['signal_sequence']['length'] > mon_len:
                    raise ValidationError({'signal_sequence': 'Length must be less than protein length'})

class Reaction(SpeciesComponent):
    #parent pointer
    parent_ptr_species_component = OneToOneField(SpeciesComponent, related_name='child_ptr_reaction', parent_link=True, verbose_name='Species component')
    
    #additional fields
    stoichiometry = ManyToManyField(ReactionStoichiometryParticipant, related_name='reactions', verbose_name='Stoichiometry')
    direction = CharField(max_length=1, choices=CHOICES_REACTION_DIRECTION, verbose_name='Direction') 
    modification = ForeignKey(ModificationReaction, blank=True, null=True, on_delete=SET_NULL, related_name='reactions', verbose_name='Modification')
    enzyme = ForeignKey(EnzymeParticipant, blank=True, null=True, on_delete=SET_NULL, related_name='reactions', verbose_name='Enzyme')
    coenzymes = ManyToManyField(CoenzymeParticipant, blank=True, null=True, related_name='reactions', verbose_name='Coenzymes')
    is_spontaneous = BooleanField(verbose_name='Is spontaneous (pH 7.5, 25C, <i>I</i> = 0)')
    delta_g = FloatField(blank=True, null=True, verbose_name='&Delta;G (pH 7.5, 25C, <i>I</i> = 0; kJ mol<sup>-1</sup>)')
    keq = ForeignKey(EntryPositiveFloatData, blank=True, null=True, on_delete=SET_NULL, verbose_name='K<sub>eq</sub>', related_name='+')
    kinetics_forward = ForeignKey(Kinetics, blank=True, null=True, on_delete=SET_NULL, related_name='reactions_forward', verbose_name='Forward kinetics')
    kinetics_backward = ForeignKey(Kinetics, blank=True, null=True, on_delete=SET_NULL, related_name='reactions_backward', verbose_name='Backward kinetics')
    optimal_ph = ForeignKey(EntryPositiveFloatData, blank=True, null=True, on_delete=SET_NULL, verbose_name='Optimal pH', related_name='+')
    optimal_temperature = ForeignKey(EntryFloatData, blank=True, null=True, on_delete=SET_NULL, verbose_name='Optimal temperature', related_name='+')        
    pathways = ManyToManyField('Pathway', blank=True, null=True, related_name='reactions', verbose_name='Pathways')
    processes = ForeignKey('Process', blank=True, null=True, on_delete=SET_NULL, related_name='reactions', verbose_name='Process')
    states = ForeignKey('State', blank=True, null=True, on_delete=SET_NULL, related_name='reactions', verbose_name='State')
    map_coordinates = ManyToManyField(ReactionMapCoordinate, blank=True, null=True, related_name='reactions', verbose_name='Map coordinates')

    #getters
            
    #html formatting
    def get_as_html_stoichiometry(self, species, is_user_anonymous):
        compartments = []
        for s in self.stoichiometry.all():
            compartments.append(s.compartment)
        compartments = list(set(compartments))
            
        pos = []
        neg = []
        for s in self.stoichiometry.all():
            if s.coefficient < 0:
                tmp = ''
                if s.coefficient != -1:
                    tmp += '(%d) ' % -s.coefficient
                tmp += '<a href="%s">%s</a>' % (s.molecule.get_absolute_url(species), s.molecule.name)
                if len(compartments) > 1:
                    tmp += '[<a href="%s">%s</a>]' % (s.compartment.get_absolute_url(species), s.compartment.wid)
                pos.append(tmp)
            else:
                tmp = ''
                if s.coefficient != 1:
                    tmp += '(%d) ' % s.coefficient
                tmp += '<a href="%s">%s</a>' % (s.molecule.get_absolute_url(species), s.molecule.name)
                if len(compartments) > 1:
                    tmp += '[<a href="%s">%s</a>]' % (s.compartment.get_absolute_url(species), s.compartment.wid)
                neg.append(tmp)
                    
        result = ''
        if len(compartments) == 1:
            result += '[<a href="%s">%s</a>]: ' % (compartments[0].get_absolute_url(species), compartments[0].wid)
        result += ' + '.join(pos)
        if self.direction == 'f':
            result += ' &rArr; '
        elif self.direction == 'b':
            result += ' &lArr; '
        elif self.direction == 'r':
            result += ' &hArr; '
        result += ' + '.join(neg)
        return format_with_evidence(obj = self.stoichiometry.all(), txt = result)
    
    def get_as_html_modification(self, species, is_user_anonymous):
        m = self.modification
        if m is None:
            return
        if m.position is None:
            txt = '<a href="%s">%s</a> [<a href="%s">%s</a>]' % (m.molecule.get_absolute_url(species), m.molecule.wid, m.compartment.get_absolute_url(species), m.compartment.wid)
        else:
            txt = '(%d) <a href="%s">%s</a> [<a href="%s">%s</a>]' % (m.position, m.molecule.get_absolute_url(species), m.molecule.wid, m.compartment.get_absolute_url(species), m.compartment.wid)
            
        return format_with_evidence(obj = m, txt = txt)
    
    def get_as_html_enzyme(self, species, is_user_anonymous):
        e = self.enzyme
        if e is None:
            return
        return format_with_evidence(obj = e, txt = '<a href="%s">%s</a> [<a href="%s">%s</a>]' % (e.protein.get_absolute_url(species), e.protein.wid, e.compartment.get_absolute_url(species), e.compartment.wid))
        
    def get_as_html_coenzymes(self, species, is_user_anonymous):
        results = []
        for c in self.coenzymes.all():
            if c.coefficient is None:
                results.append(format_with_evidence(list_item = True, obj = c, txt = '<a href="%s">%s</a> [<a href="%s">%s</a>]' % (c.metabolite.get_absolute_url(species), c.metabolite.wid, c.compartment.get_absolute_url(species), c.compartment.wid)))
            else:
                results.append(format_with_evidence(list_item = True, obj = c, txt = '(%d) <a href="%s">%s</a> [<a href="%s">%s</a>]' % (c.coefficient, c.metabolite.get_absolute_url(species), c.metabolite.wid, c.compartment.get_absolute_url(species), c.compartment.wid)))
                
        return format_list_html(results, force_list=True)
        
    def get_as_html_kinetics_forward(self, species, is_user_anonymous):
        k = self.kinetics_forward
        if k is None:
            return
        law = ''
        kms_html = ''
        vmax = ''
        if k.rate_law is not None and k.rate_law != '':
            law = k.rate_law
            law = law \
                .replace(' ', '') \
                .replace('*', ' * ') \
                .replace('+', ' + ') \
                .replace('-', ' - ') \
                .replace('/', ' / ')
            law = '<i>v</i> = %s;' % re.sub(r'([a-z0-9_]+)', sub_rate_law(self.species.wid), law, flags=re.I)            
        if k.km != '':
            kms = k.km.split(', ')
            if len(kms) == 1:
                kms_html = '<i>K</i><sub>m</sub> = %s (nM), ' % (kms[0], )
            else:            
                kms_html = ''
                for i in range(len(kms)):
                    kms_html += '<i>K</i><sub>m%s</sub> = %s (nM), ' % (i+1, kms[i], )
        if k.vmax is not None:
            vmax = '<i>V</i><sub>max</sub> = %s %s' % (k.vmax, k.vmax_unit, )
        return format_with_evidence(obj = k, txt = '%s %s %s' % (law, kms_html, vmax))
        
    def get_as_html_kinetics_backward(self, species, is_user_anonymous):
        k = self.kinetics_backward
        if k is None:
            return
        law = ''
        km = ''
        vmax = ''
        if k.rate_law is not None and k.rate_law != '':
            law = k.rate_law
            law = law \
                .replace(' ', '') \
                .replace('*', ' * ') \
                .replace('+', ' + ') \
                .replace('-', ' - ') \
                .replace('/', ' / ')
            law = '<i>v</i> = %s;' % re.sub(r'([a-z0-9_]+)', sub_rate_law(self.species.wid), law, flags=re.I)
        if k.km != '':
            kms = k.km.split(', ')            
            if len(kms) == 1:
                kms_html = '<i>K</i><sub>m</sub> = %s (nM), ' % (kms[0], )
            else:            
                kms_html = ''
                for i in range(len(kms)):
                    kms_html += '<i>K</i><sub>m%s</sub> = %s (nM), ' % (i+1, kms[i], )
        if k.vmax is not None:
            vmax = '<i>V</i><sub>max</sub> = %s %s' % (k.vmax, k.vmax_unit, )
        return format_with_evidence(obj = k, txt = '%s %s %s' % (law, km, vmax))
        
    #meta information
    class Meta:
        concrete_entry_model = True
        fieldsets = [
            ('Type', {'fields': ['model_type']}),
            ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), 
            ('Classification', {'fields': ['type']}),            
            ('Reaction', {'fields': ['stoichiometry', 'modification']}), 
            ('Catalysis', {'fields': ['enzyme', 'coenzymes', 'optimal_ph', 'optimal_temperature']}), 
            ('Energetics', {'fields': ['is_spontaneous', 'delta_g', 'keq']}), 
            ('Kinetics', {'fields': ['kinetics_forward', 'kinetics_backward']}), 
            ('Parameters', {'fields': ['parameters']}),
            ('Associations', {'fields': ['pathways', 'processes', 'states']}),
            ('Comments', {'fields': ['comments', 'publication_references']}),
            ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]}),
            ]
        field_list = [
            'id', 'wid', 'name', 'synonyms', 'cross_references',
            'type',             
            'stoichiometry', 'direction', 'modification', 
            'enzyme', 'coenzymes', 'optimal_ph', 'optimal_temperature', 
            'is_spontaneous', 'delta_g', 'keq', 
            'kinetics_forward', 'kinetics_backward',
            'pathways', 'processes', 'states',
            'map_coordinates',
            'comments',
            'publication_references', 
            'created_user', 'created_date', 'last_updated_user', 'last_updated_date', 
            ]
        facet_fields = ['type', 'direction', 'enzyme__protein', 'coenzymes__metabolite', 'is_spontaneous', 'pathways', 'processes', 'states']
        verbose_name='Reaction'
        verbose_name_plural = 'Reactions'    

        def clean(model, obj_data, all_obj_data=None, all_obj_data_by_model=None):
            from cyano.helpers import getEntry, getModel, EmpiricalFormula
            
            #stoichiometry
            formula = EmpiricalFormula()
            includes_macromolecules = False
            for s in obj_data['stoichiometry']:
                if all_obj_data is None:
                    molecule = getEntry(species_wid=obj_data['species'], wid=s['molecule'])
                else:
                    molecule = all_obj_data[s['molecule']]
                if isinstance(molecule, Entry):
                    molecule_type = molecule.model_type
                else:
                    molecule_type = molecule['model_type']
                if molecule_type == 'Metabolite':
                    if isinstance(molecule, Entry):
                        molecule_formula = EmpiricalFormula(molecule.empirical_formula)
                    else:
                        molecule_formula = EmpiricalFormula(molecule['empirical_formula'])
                else:
                    molecule_formula = EmpiricalFormula() #todo: implement
                
                formula += molecule_formula * s['coefficient']
                    
                if molecule_type not in ['Metabolite', 'Stimulus']:
                    includes_macromolecules = True
                if obj_data['modification'] is not None:
                    includes_macromolecules = True
            
            if len(formula) > 0 and not includes_macromolecules: #todo: remove "includes_macromolecules"
                raise ValidationError({'stoichiometry': 'Reaction imbalanced by %s' % formula.get_as_html()})
                
            #kinetics
            if obj_data['kinetics_forward']    is not None:
                validate_kinetics(obj_data, 'f')
            if obj_data['kinetics_backward'] is not None:
                validate_kinetics(obj_data, 'r')
            
            #modication
            if obj_data['modification'] is not None:
                mod = obj_data['modification']
                
                if all_obj_data is None:
                    molecule = getEntry(species_wid=obj_data['species'], wid=mod['molecule'])
                else:
                    molecule = all_obj_data[mod['molecule']]
                if isinstance(molecule, Entry):
                    molecule_type = molecule.model_type
                else:
                    molecule_type = molecule['model_type']
                molecule_len = None
                if molecule_type == 'Gene':
                    if isinstance(molecule, Entry):
                        molecule_len = molecule.length
                    else:
                        molecule_len = molecule['length']                    
                elif molecule_type == 'ProteinMonomer':
                    if isinstance(molecule, Entry):
                        gene_wid = molecule.gene.wid
                    else:
                        gene_wid = molecule['gene']
                    if all_obj_data is None:
                        gene = Gene.objects.get(species__wid=obj_data['species'], wid=gene_wid)
                    else:
                        gene = all_obj_data[gene_wid]
                    if isinstance(gene, Entry):
                        molecule_len = gene.length / 3
                    else:
                        molecule_len = gene['length'] / 3
                
                if mod['position'] is not None and not issubclass(getModel(molecule_type), (ProteinMonomer, Gene, )):
                    raise ValidationError({'modification': 'Position must be null'})
                
                if mod['position'] is not None and mod['position'] > molecule_len:
                    raise ValidationError({'modification': 'Position must be less than molecule length'})

class Species(Entry):
    #parent pointer
    parent_ptr_entry = OneToOneField(Entry, related_name='child_ptr_species', parent_link=True, verbose_name='Entry')
    
    #additional fields
    genetic_code = CharField(max_length=50, verbose_name='Genetic code', choices = CHOICES_GENETIC_CODE)
    
    cross_references = ManyToManyField("CrossReference", blank=True, null=True, related_name='cross_referenced_species', verbose_name='Cross references')
    publication_references = ManyToManyField("PublicationReference", blank=True, null=True, related_name='publication_referenced_species', verbose_name='Publication references')
    
    #getters
    @permalink
    def get_absolute_url(self):
        return ('cyano.views.species', (), {'species_wid': self.wid})
    
    #html formatting
    def get_as_html_comments(self, species, is_user_anonymous):
        txt = self.comments
        
        #provide links to references
        return re.sub(r'\[(PUB_\d{4,4})(, PUB_\d{4,4})*\]', 
            lambda match: '[' + ', '.join(['<a href="%s">%s</a>' % (reverse('cyano.views.detail', kwargs={'species_wid':self.wid, 'wid': x}), x, ) for x in match.group(0)[1:-1].split(', ')]) + ']',
            txt)
            
    def get_as_html_genetic_code(self, species, is_user_anonymous):
        return '<a href="http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG%s">%s</a>' % (self.genetic_code, self.genetic_code, )
    
    #meta information
    class Meta:
        concrete_entry_model = True
        fieldsets = [
            ('Type', {'fields': ['model_type']}), 
            ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}),             
            ('Physiology', {'fields': ['genetic_code']}), 
            ('Comments', {'fields': ['comments']}), 
            ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]}),
            ]        
        field_list = [
            'id', 'wid', 'name', 'synonyms', 'cross_references', 'genetic_code',
            'comments'
            ]
        facet_fields = []
        verbose_name='Species'
        verbose_name_plural = 'Species'

class State(SpeciesComponent):
    #parent pointer
    parent_ptr_species_component = OneToOneField(SpeciesComponent, related_name='child_ptr_state', parent_link=True, verbose_name='Species component')
    
    #getters
    
    #html formatting
    
    #meta information
    class Meta:
        concrete_entry_model = True
        fieldsets = [
            ('Type', {'fields': ['model_type']}),
            ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), 
            ('Classification', {'fields': ['type']}),             
            ('Reactions', {'fields': ['reactions']}),
            ('Parameters', {'fields': ['parameters']}),
            ('Comments', {'fields': ['comments', 'publication_references']}),
            ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]}),
            ]        
        field_list = [
            'id', 'wid', 'name', 'synonyms', 'cross_references',
            'type', 
            'comments',
            'publication_references', 
            'created_user', 'created_date', 'last_updated_user', 'last_updated_date', 
            ]
        facet_fields = ['type']
        verbose_name='State'
        verbose_name_plural = 'States'
        
class Stimulus(Molecule):
    #parent pointer
    parent_ptr_molecule = OneToOneField(Molecule, related_name='child_ptr_stimulus', parent_link=True, verbose_name='Molecule')
    
    #additional fields    
    value = ForeignKey(EntryFloatData, verbose_name='Value', related_name='+')
    
    #getters
    
    #html formatting    
    
    #meta information
    class Meta:
        concrete_entry_model = True
        fieldsets = [
            ('Type', {'fields': ['model_type']}),
            ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), 
            ('Classification', {'fields': ['type']}),
            ('Value', {'fields': ['value']}), 
            ('Function', {'fields': [
                {'verbose_name': 'Reaction participant', 'name':'reaction_stoichiometry_participants'},
                {'verbose_name': 'Complex subunit', 'name':'protein_complex_biosythesis_participants'},
                ]}),
            ('Parameters', {'fields': ['parameters']}),
            ('Comments', {'fields': ['comments', 'publication_references']}),
            ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]}),
            ]
        field_list = [
            'id', 'wid', 'name', 'synonyms', 'cross_references',
            'type', 
            'value', 
            'comments',
            'publication_references', 
            'created_user', 'created_date', 'last_updated_user', 'last_updated_date', 
            ]
        facet_fields = ['type', 'value__units']
        verbose_name='Stimulus'
        verbose_name_plural = 'Stimuli'
    
class TranscriptionUnit(Molecule):
    #parent pointer
    parent_ptr_molecule = OneToOneField(Molecule, related_name='child_ptr_transcription_unit', parent_link=True, verbose_name='Molecule')
    
    #additional fields
    genes = ManyToManyField('Gene', related_name='transcription_units', verbose_name='Genes')
    promoter_35_coordinate = IntegerField(null=True, blank=True, verbose_name='Promoter -35 box coordinate (nt)')
    promoter_35_length = IntegerField(null=True, blank=True, verbose_name='Promoter -35 box length (nt)')
    promoter_10_coordinate = IntegerField(null=True, blank=True, verbose_name='Promoter -10 box coordinate (nt)')
    promoter_10_length = IntegerField(null=True, blank=True, verbose_name='Promoter -10 box length (nt)')
    tss_coordinate = IntegerField(null=True, blank=True, verbose_name='Transcription start site coordinate (nt)')
    
    #getters
    def get_chromosome(self):
        chr = list(set([g.chromosome for g in self.genes.all()]))
        if len(chr) != 1:
            raise        
        return chr[0]
        
    def get_coordinate(self):
        return self.genes.all().aggregate(Min('coordinate'))['coordinate__min']
        
    def get_length(self):
        return max(self.genes.extra(select={"end_coordinate": "(coordinate + length - 1)"}).values('end_coordinate'))['end_coordinate'] - self.get_coordinate() + 1
        
    def get_direction(self):
        dir = list(set([g[0] for g in self.genes.values_list('direction').all()]))
        if len(dir) != 1:
            raise        
        return dir[0]
        
    def get_sequence(self):
        seq = self.get_chromosome().sequence[self.get_coordinate() - 1:self.get_coordinate() - 1 + self.get_length()]
        if self.get_direction() == 'r':
            seq = unicode(Seq(seq, IUPAC.unambiguous_dna).reverse_complement())
        return seq
        
    def get_gc_content(self):
        seq = self.get_sequence()
        return float(seq.count('G') + seq.count('C')) / float(len(seq))
        
    #html formatting
    def get_as_html_structure(self, species, is_user_anonymous):
        return self.get_chromosome().get_as_html_structure(species, is_user_anonymous, 
            zoom = 1, 
            start_coordinate = self.get_coordinate() - 500, 
            end_coordinate = self.get_coordinate() + self.get_length() + 500, 
            highlight_wid = [self.wid] + [g.wid for g in self.genes.all()])
            
    def get_as_html_genes(self, species, is_user_anonymous):
        results = []
        for g in self.genes.all():
            results.append('<a href="%s">%s</a>' % (g.get_absolute_url(species), g.wid, ))
        return format_list_html(results, comma_separated=True)
        
    def get_as_html_sequence(self, species, is_user_anonymous):
        from cyano.helpers import format_sequence_as_html
        
        direction = CHOICES_DIRECTION[[x[0] for x in CHOICES_DIRECTION].index(self.get_direction())][1]        
        
        return 'Chromosome: <a href="%s">%s</a>, Coordinate: %s, Length: %s, Direction: %s, Sequence: %s' % (
            self.get_chromosome().get_absolute_url(species), self.get_chromosome().wid, 
            self.get_coordinate(), self.get_length(), direction, 
            format_sequence_as_html(self.get_sequence()))
        
    def get_as_html_transcriptional_regulations(self, species, is_user_anonymous):
        results = []
        for r in self.transcriptional_regulations.all():
            results.append('<a href="%s">%s</a>: <a href="%s">%s</a>' % (r.get_absolute_url(species), r.wid, r.transcription_factor.get_absolute_url(species), r.transcription_factor.wid))
        return format_list_html(results)

    #meta information
    class Meta:
        concrete_entry_model = True
        fieldsets = [
            ('Type', {'fields': ['model_type']}),
            ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), 
            ('Classification', {'fields': ['type']}),
            ('Structure (Hayflick media, 37C)', {'fields': [
                {'verbose_name': 'Structure', 'name': 'structure'},
                'genes', 
                'promoter_35_coordinate', 
                'promoter_35_length', 
                'promoter_10_coordinate', 
                'promoter_10_length', 
                'tss_coordinate', 
                {'verbose_name': 'Sequence', 'name': 'sequence'},
                ]}),    
            ('Regulation', {'fields': [{'verbose_name': 'Regulation', 'name': 'transcriptional_regulations'}]}),
            ('Function', {'fields': [
                {'verbose_name': 'Reaction participant', 'name':'reaction_stoichiometry_participants'},
                {'verbose_name': 'Complex subunit', 'name':'protein_complex_biosythesis_participants'},
                ]}),
            ('Parameters', {'fields': ['parameters']}),
            ('Comments', {'fields': ['comments', 'publication_references']}),
            ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]}),
            ]
        field_list = [
            'id', 'wid', 'name', 'synonyms', 'cross_references',
            'type', 
            'genes', 'promoter_35_coordinate', 'promoter_35_length', 'promoter_10_coordinate', 'promoter_10_length', 'tss_coordinate',    
            'comments',
            'publication_references', 
            'created_user', 'created_date', 'last_updated_user', 'last_updated_date', 
            ]
        facet_fields = ['type']
        verbose_name='Transcription unit'
        verbose_name_plural = 'Transcription units'
        
        def clean(model, obj_data, all_obj_data=None, all_obj_data_by_model=None):
            if len(obj_data['genes']) == 1:
                return
            if len(obj_data['genes']) == 0:
                raise ValidationError({'genes': 'Transcription units most contain at least 1 gene'})
                
            chr_wids = []
            if all_obj_data is None:
                for gene_wid in obj_data['genes']:
                    chr_wids.append(Gene.objects.get(species__wid=obj_data['species'], wid=gene_wid).chromosome.wid)
            else:
                for gene_wid in obj_data['genes']:
                    if isinstance(all_obj_data[gene_wid], Entry):
                        chr_wids.append(all_obj_data[gene_wid].chromosome.wid)
                    else:
                        chr_wids.append(all_obj_data[gene_wid]['chromosome'])
                    
            if len(set(chr_wids)) > 1:
                raise ValidationError({'genes': 'Genes must all belong to the same chromosome'})
        
class TranscriptionalRegulation(SpeciesComponent):
    #parent pointer
    parent_ptr_species_component = OneToOneField(SpeciesComponent, related_name='child_ptr_transcriptional_regulation', parent_link=True, verbose_name='Species component')
    
    #additional fields
    transcription_unit = ForeignKey('TranscriptionUnit', related_name='transcriptional_regulations', verbose_name='Transcription unit')
    transcription_factor = ForeignKey('Protein', related_name='transcriptional_regulations',  verbose_name='Transcripton factor')
    binding_site = ForeignKey(BindingSite, null=True, blank=True, on_delete=SET_NULL, related_name='transcriptional_regulations', verbose_name='Binding site')
    affinity = ForeignKey(EntryPositiveFloatData, null=True, blank=True, on_delete=SET_NULL, verbose_name='Affinity', related_name='+')
    activity = ForeignKey(EntryPositiveFloatData, verbose_name='Fold-change activity', related_name='+')
    
    #getters
    def get_binding_site_sequence(self):
        bs = self.binding_site
        if bs is None:
            return None
            
        seq = self.transcription_unit.get_chromosome().sequence[bs.coordinate - 1:bs.coordinate - 1 + bs.length]
        if bs.direction == 'r':
            seq = seq[::-1] \
                .replace('A', 't') \
                .replace('C', 'g') \
                .replace('G', 'C') \
                .replace('T', 'A') \
                .replace('t', 'T') \
                .replace('g', 'G')
        return seq
    
    #html formatting
    def get_as_html_binding_site(self, species, is_user_anonymous):
        from cyano.helpers import format_sequence_as_html
        
        bs = self.binding_site
        if bs is None:
            return None
        
        direction = CHOICES_DIRECTION[[x[0] for x in CHOICES_DIRECTION].index(bs.direction)][1]
        
        chr = self.transcription_unit.get_chromosome()
        
        map = chr.get_as_html_structure(is_user_anonymous,
                zoom = 1, 
                start_coordinate = bs.coordinate - 500, 
                end_coordinate = bs.coordinate + bs.length + 500, 
                highlight_wid = [self.wid])
            
        txt = '%s<br/>Chromosome: <a href="%s">%s</a>, Coordinate: %s (nt), Length: %s (nt), Direction: %s, Sequence: %s' % (
            map, chr.get_absolute_url(species), chr.wid, 
            bs.coordinate, bs.length, direction, 
            format_sequence_as_html(self.get_binding_site_sequence()))
            
        return format_with_evidence(obj = bs, txt = txt)

    #meta information
    class Meta:
        concrete_entry_model = True
        fieldsets = [
            ('Type', {'fields': ['model_type']}),
            ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), 
            ('Classification', {'fields': ['type']}),
            ('Regulation', {'fields': [
                'transcription_unit', 
                'transcription_factor', 
                'binding_site', 
                'affinity', 
                'activity'
                ]}), 
            ('Comments', {'fields': ['comments', 'publication_references']}),
            ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]}),
            ]
        field_list = [
            'id', 'wid', 'name', 'synonyms', 'cross_references',
            'type', 
            'transcription_unit', 'transcription_factor', 'binding_site', 'affinity', 'activity', 
            'comments',
            'publication_references', 
            'created_user', 'created_date', 'last_updated_user', 'last_updated_date', 
            ]
        facet_fields = ['type', 'transcription_unit', 'transcription_factor']
        verbose_name='Transcriptional regulation'
        verbose_name_plural = 'Transcriptional regulation'
        
        def clean(model, obj_data, all_obj_data=None, all_obj_data_by_model=None):
            #gene wid
            if all_obj_data is None:                
                tu = TranscriptionUnit.objects.get(species__wid=obj_data['species'], wid=obj_data['transcription_unit'])
            else:
                tu = all_obj_data[obj_data['transcription_unit']]
            if isinstance(tu, Entry):
                gene_wid = tu.genes.all()[0].wid
            else:
                gene_wid = tu['genes'][0]
                
            #chr wid
            if all_obj_data is None:
                gene = Gene.objects.get(species__wid=obj_data['species'], wid=gene_wid)
            else:
                gene = all_obj_data[gene_wid]
            if isinstance(gene, Entry):
                chr_wid = gene.chromosome.wid
            else:
                chr_wid = gene['chromosome']
                
            #chr length
            if all_obj_data is None:
                chr = Chromosome.objects.get(species__wid=obj_data['species'], wid=chr_wid)
            else:
                chr = all_obj_data[chr_wid]
            if isinstance(chr, Entry):
                chr_len = chr.length
            else:
                chr_len = chr['length']
            
            #error check binding site coordinate, length
            if obj_data['binding_site'] is not None:
                if obj_data['binding_site']['coordinate'] > chr_len:
                    raise ValidationError({'binding_site': 'Coordinate must be less then chromosome length.'})
                if obj_data['binding_site']['length'] > chr_len:
                    raise ValidationError({'binding_site': 'Length must be less then chromosome length.'})
        
class Type(SpeciesComponent):
    #parent pointer
    parent_ptr_species_component = OneToOneField(SpeciesComponent, related_name='child_ptr_type', parent_link=True, verbose_name='Species component')
    
    #additional fields
    parent = ForeignKey('self', blank=True, null=True, on_delete=SET_NULL, related_name='children', verbose_name='Parent')

    #getters
    def get_all_members(self, species):
        members = []
        for m in self.members.filter(species = species):
            members.append(m)
        for c in self.children.filter(species = species):
            members += c.get_all_members()
        return members

    #html formatting    
    def get_as_html_parent(self, species, is_user_anonymous):
        if self.parent is not None:
            result = '<a href="%s">%s</a>' % (self.parent.get_absolute_url(species), self.parent.wid, )
            if self.parent.parent is not None:
                result = self.parent.get_as_html_parent(is_user_anonymous) + ' &#8250; ' + result
            return result
        
    def get_as_html_children(self, species, is_user_anonymous):
        results = []
        for c in self.children.filter(species = species):
            results.append('<a href="%s">%s</a>' % (c.get_absolute_url(species), c.wid))        
        return format_list_html(results, comma_separated=True)
        
    def get_as_html_members(self, species, is_user_anonymous):
        results = []
        for m in self.get_all_members(species):
            results.append('<a href="%s">%s</a>' % (m.get_absolute_url(species), m.wid))        
        return format_list_html(results, comma_separated=True)
    
    #meta information
    class Meta:
        concrete_entry_model = True
        fieldsets = [
            ('Type', {'fields': ['model_type']}),
            ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), 
            ('Classification', {'fields': ['type', 'parent', 'children', 'members']}),
            ('Comments', {'fields': ['comments', 'publication_references']}),
            ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]}),
            ]
        field_list = [
            'id', 'wid', 'name', 'synonyms', 'cross_references',
            'type', 'parent',
            'comments',
            'publication_references', 
            'created_user', 'created_date', 'last_updated_user', 'last_updated_date', 
            ]
        facet_fields = ['type', 'parent']
        verbose_name='Type'
        verbose_name_plural = 'Types'
            
class CrossReference(SpeciesComponent):
    #parent pointer
    parent_ptr_species_component = OneToOneField(SpeciesComponent, related_name='child_ptr_crossreference', parent_link=True, verbose_name='Species component')
    
    #additional fields
    xid = CharField(max_length=255, verbose_name='External ID', blank=True)
    source = CharField(max_length=20, choices=CHOICES_CROSS_REFERENCE_SOURCES, verbose_name='Source', blank=True)

    #getters
    def get_all_members(self, species):
        return self.cross_referenced_entries.filter(species = species)
    
    def get_as_html_members(self, species, is_user_anonymous):
        results = []
        for m in self.get_all_members(species):
            results.append('<a href="%s">%s</a>' % (m.get_absolute_url(species), m.wid))        
        return format_list_html(results, comma_separated=True)

    class Meta:
        concrete_entry_model = True
        fieldsets = [
            ('Type', {'fields': ['model_type']}),
            ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}),
            ('Classification', {'fields': [{'verbose_name': "Referenced by", 'name': 'members'}]}),
            ('Comments', {'fields': ['comments', 'publication_references']}),
            ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]}),
            ]
        field_list = [
            'id', 'wid', 'name', 'synonyms', 'cross_references',
            'comments',
            'publication_references', 
            'created_user', 'created_date', 'last_updated_user', 'last_updated_date', 
            ]
        facet_fields = ['type', 'source']
        ordering = ['xid']
        verbose_name='Cross reference'
        verbose_name_plural = 'Cross references'

class PublicationReference(SpeciesComponent):
    #parent pointer
    parent_ptr_species_component = OneToOneField(SpeciesComponent, related_name='child_ptr_publicationreference', parent_link=True, verbose_name='Species component')
    
    #additional fields
    authors = TextField(blank=True, default='', verbose_name='Author(s)')
    editors = TextField(blank=True, default='', verbose_name='Editor(s)')
    year = PositiveIntegerField(blank=True, null=True, verbose_name='Year')
    title = TextField(blank=True, default='', verbose_name='Title')
    publication = CharField(max_length=255, blank=True, default='', verbose_name='Publication')
    publisher = CharField(max_length=255, blank=True, default='', verbose_name='Publisher')
    volume = CharField(max_length=255, blank=True, default='', verbose_name='Volume')
    issue = CharField(max_length=255, blank=True, default='', verbose_name='Issue')
    pages = CharField(max_length=255, blank=True, default='', verbose_name='Page(s)')
    
    #getters
    def get_citation(self, species, cross_references = False):
        if self.type.exists():
            if self.type.all()[0].wid == 'article':
                txt = '%s. %s. <i>%s</i> <b>%s</b>, %s (%s).' % (self.authors, self.title, self.publication, self.volume, self.pages, self.year, )
            elif self.type.all()[0].wid == 'book':
                authors = ''
                editors = ''
                if self.authors != '':
                    authors = '%s.' % self.authors 
                if self.editors != '':
                    editors = 'Eds %s.' % self.editors
                txt = '%s %s <i>%s</i>. %s %s (%s).' % (authors, editors, self.title, self.publisher, self.pages, self.year)
            elif self.type.all()[0].wid == 'thesis':
                txt = '%s. <i>%s</i>. %s (%s).' % (self.authors, self.title, self.publisher, self.year)
            else:
                txt = '%s. <i>%s</i>. (%s).' % (self.authors, self.title, self.year)
        else:
            txt = '%s. <i>%s</i>. (%s) %s.' % (self.authors, self.title, self.publication, self.year)

        cr = self.get_as_html_cross_references(species, True)
        cr_spacer = ''
        if cr != '':
            cr_spacer = ', '
        return '%s CyanoFactory: <a href="%s">%s</a>%s%s' % (txt, self.get_absolute_url(species), self.wid, cr_spacer, cr, )
            
    def get_all_referenced_entries(self, species):
        entries = []
        for entry in self.publication_referenced_entries.filter(species = species):
            entries.append(entry)
        for ev in Evidence.objects.filter(references__id=self.id):
            entries.append(ev.species_component)
        return entries
    
    #html formatting    
    def get_as_html_citation(self, species, is_user_anonymous):
        return self.get_citation(species)
        
    def get_as_html_referenced_entries(self, species, is_user_anonymous):
        results = []
        for o in self.get_all_referenced_entries(species):
            results.append('<a href="%s">%s</a>' % (o.get_absolute_url(species), o.wid))
        return format_list_html(results)
        
    def get_as_bibtex(self):
        type = None
        props = []        
        if self.type.all()[0].wid == 'article':
            type = 'ARTICLE'            
            if self.authors != '':
                props.append(('AUTHOR', self.format_authors_bibtex(self.authors)))
            if self.title != '':
                props.append(('TITLE', self.title))
            if self.publication != '':
                props.append(('JOURNAL', self.publication))
            if self.year is not None:
                props.append(('YEAR', self.year))
            if self.volume != '':
                props.append(('VOLUME', self.volume))
            if self.issue != '':
                props.append(('NUMBER', self.issue))
            if self.pages != '':
                props.append(('PAGES', self.pages))
            for cr in self.cross_references.all():
                if cr.source == 'PubMed':
                    props.append(('eprint', cr.xid))
                    props.append(('eprinttype', 'pubmed'))
                    break
            for cr in self.cross_references.all():
                if cr.source == 'URL':
                    props.append(('URL', cr.xid))
                    break
        elif self.type.all()[0].wid == 'book':
            type = 'BOOK'
            if self.authors != '':
                props.append(('AUTHOR', self.format_authors_bibtex(self.authors)))
            elif self.editors != '':
                props.append(('EDITOR', self.format_authors_bibtex(self.editors)))
            if self.title != '':
                props.append(('TITLE', self.title))
            if self.year is not None:
                props.append(('YEAR', self.year))
            if self.volume != '':
                props.append(('VOLUME', self.volume))
            if self.publisher != '':
                props.append(('PUBLISHER', self.publisher))
            for cr in self.cross_references.all():
                if cr.source == 'ISBN':
                    props.append(('ISBN', cr.xid))
                    break
            for cr in self.cross_references.all():
                if cr.source == 'URL':
                    props.append(('URL', cr.xid))
                    break
        elif self.type.all()[0].wid == 'thesis':
            type = 'THESIS'
            if self.authors != '':
                props.append(('AUTHOR', self.format_authors_bibtex(self.authors)))
            if self.editors != '':
                props.append(('TITLE', self.editors))
            if self.year is not None:
                props.append(('YEAR', self.year))
            if self.publisher is not None:
                props.append(('SCHOOL', self.publisher))
            for cr in self.cross_references.all():
                if cr.source == 'URL':
                    props.append(('URL', cr.xid))
                    break
        else:
            type = 'MISC'
            if self.authors != '':
                props.append(('AUTHOR', self.format_authors_bibtex(self.authors)))
            if self.editors != '':
                props.append(('TITLE', self.editors))
            if self.year is not None:
                props.append(('YEAR', self.year))
            for cr in self.cross_references.all():
                if cr.source == 'URL':
                    props.append(('URL', cr.xid))
                    break
            
        props.append(('created', set_time_zone(self.created_date).isoformat()))
        props.append(('lastUpdated', set_time_zone(self.last_updated_date).isoformat()))
        
        tmp = []
        for prop in props:
            if prop[0] == 'TITLE':
                tmp.append('\n\t%s = "{%s}",' % prop)
            else:
                tmp.append('\n\t%s = {%s},' % prop)
        return '@%s{%s%s\n}' % (type, self.wid, ''.join(tmp))
        
    def format_authors_bibtex(self, authors):
        authors = authors.split(", ")
        for idx in range(len(authors)):
            author = authors[idx]
            
            names = author.split(" ")
            for fNameIdx in range(len(names)):
                if names[fNameIdx].upper() == names[fNameIdx]:
                    break    
            
            tmpFirstName = names[fNameIdx]
            firstName = ""
            for i in range(len(tmpFirstName)):
                firstName += tmpFirstName[i] + ". "        
            firstName = firstName.strip()
            
            lastName = " ".join(names[0:fNameIdx])

            suffix = " ".join(names[fNameIdx + 1:len(names)])
            
            authors[idx] = lastName        
            if firstName != '':
                authors[idx] += ", " + firstName
                if suffix != '':
                    authors[idx] += " " + suffix
            
        return " and ".join(authors)
    
    #meta information
    class Meta:        
        concrete_entry_model = True
        fieldsets = [
            ('Type', {'fields': ['model_type']}),
            ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), 
            ('Classification', {'fields': ['type']}),            
            ('Citation', {'fields': [{'verbose_name': 'Citation', 'name': 'citation'}]}), 
            ('Cited by', {'fields': [{'verbose_name': 'Cited by', 'name': 'referenced_entries'}]}),
            ('Comments', {'fields': ['comments', 'publication_references']}),
            ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]}),
            ]            
        field_list = [
            'id', 'wid', 'name', 'synonyms', 'cross_references',
            'type',             
            'authors', 'editors', 'year', 'title', 'publication', 'publisher', 'volume', 'issue', 'pages', 
            'comments',
            'publication_references', 
            'created_user', 'created_date', 'last_updated_user', 'last_updated_date', 
            ]
        facet_fields = ['type', 'year', 'publication']
        verbose_name='Publication Reference'
        verbose_name_plural = 'Publication References'

''' END: specific data types'''

#http://isoelectric.ovh.org/files/practise-isoelectric-point.html
def calculate_nucleic_acid_pi(seq):
    numA = float(seq.count('A'))
    numC = float(seq.count('C'))
    numG = float(seq.count('G'))
    numT = float(seq.count('T'))
    numU = float(seq.count('U'))

    pH = 6.5             #starting point pI = 6.5 - theoretically it should be 7, but average protein pI is 6.5 so we increase the probability
    pHprev = 0.0         #of finding the solution
    pHnext = 14.0        #0-14 is possible pH range
    E = 0.01             #epsilon means precision [pI = pH +/- E]

    #the infinite loop
    while True:

        # http://www.steve.gb.com/science/nucleic_acids.html
        NQ = 0.
        #NQ = NQ - 1. / (1. + 10.**(3.65 - pH))     #3' charge
        NQ = NQ - numA / (1. + 10.**(3.5 - pH))   #A charge
        NQ = NQ - numC / (1. + 10.**(4.2 - pH))   #C charge
        NQ = NQ + numG / (1. + 10.**(pH - 9.2))   #G charge
        NQ = NQ - numG / (1. + 10.**(1.6 - pH))   #G charge
        NQ = NQ + numT / (1. + 10.**(pH - 9.7))   #T charge
        NQ = NQ + numU / (1. + 10.**(pH - 9.2))   #U charge
        #NQ = NQ + 1. / (1. + 10.**(pH - 8.2))      #5' charge

        if pH >= 14.0:
            raise

        #%%%%%%%%%%%%%%%%%%%%%%%%%%   BISECTION   %%%%%%%%%%%%%%%%%%%%%%%%

        #we are out of range, thus the new pH value must be smaller
        if NQ < 0.:
            temp = pH
            pH = pH - ((pH - pHprev) / 2.)
            pHnext = temp

        #we used to small pH value, so we have to increase it
        else:
            temp = pH
            pH = pH + ((pHnext - pH) / 2.)
            pHprev = temp

        #terminal condition, finding isoelectric point with given precision
        if (pH - pHprev < E) and (pHnext - pH < E):
            break

    return pH

def format_list_html(val, comma_separated=False, numbered=False, separator=None, force_list=False, vertical_spacing=False, default_items=50):
    if val is None or len(val) == 0:
        return
    if len(val) == 1 and not force_list:
        return val[0]
    class_name = ''
    
    if vertical_spacing:
        class_name = ' class="vertical_spacing"'
    
    if comma_separated:
        separator = ', '        
    
    if separator is not None:    
        if len(val) <= default_items + 1:
            return separator.join(val)
        return    separator.join(val[:default_items]) + separator\
            + ('<span><span class="button"> ... %s <a href="javascript:void(0)" onclick="$(this).parent().hide(); $(this).parent().parent().find(\'.content\').show();">more</a></span><span class="content" style="display:none;">' % (len(val) - default_items))\
            + separator.join(val[default_items:])\
            + '</span></span>'
    
    sepHeadInvisible = '<li style="display:none"><p>'
    sepHead = '<li><p>'
    sepTail = '</p></li>'
    if numbered:
        head = '<ol%s>' % class_name
        tail = '</ol>'        
    else:
        head = '<ul%s>' % class_name
        tail = '</ul>'
        
    if len(val) <= default_items + 1:
        return head + sepHead + (sepTail + sepHead).join(val) + sepTail + tail
    else:
        return head \
            + sepHead + (sepTail + sepHead).join(val[:default_items]) + sepTail \
            + sepHead + ('... %s <a href="javascript:void(0);" onclick="$(this).parent().parent().parent().find(\'li\').css(\'display\', \'list-item\'); $(this).parent().parent().hide();">more</a>' % (len(val) - default_items)) + sepTail\
            + sepHeadInvisible + (sepTail + sepHeadInvisible).join(val[default_items:]) + sepTail \
            + tail

def format_with_evidence(obj = None, txt = None, list_item = False):
    if isinstance(obj, EvidencedEntryData):
        evidence = obj.evidence.all()
    else:
        evidence = EmptyQuerySet()
        for tmp in obj:
            evidence = tmp.evidence.all() | evidence
            
    if len(evidence) == 0:
        return txt
    
    style = ''
    if list_item:
        style = 'style="margin-left:-4px;"'
    
    return '<div %s>%s <a href="javascript:void(0)" onclick="toggleEvidence(this);">(Show evidence)</a><div class="evidence" style="display:none;">%s</div></div>' % (style, txt, format_evidence(evidence), )
            
def format_evidence(evidence):
    txt = []
    for ev in evidence:
        if ev.units is None and ev.units != '':
            tmp = 'Value: %s (%s)' % (ev.value, ev.units)
        else:
            tmp = 'Value: %s' % ev.value
        
        tmp2 = []
        if ev.species is not None and ev.species != '':
            tmp2.append('Species: <i>%s</i>' % ev.species)            
        if ev.media is not None and ev.media != '':
            tmp2.append('Media: %s' % ev.media)            
        if ev.pH is not None:
            tmp2.append('pH: %s' % ev.pH)            
        if ev.temperature is not None:
            tmp2.append('Temperature (C): %s' % ev.temperature)
        
        if ev.is_experimentally_constrained:
            tmp += '<div style="margin-top:4px;"><i>Measured: Yes</i><br/>%s' % ', '.join(tmp2)
        else:
            tmp += '<div style="margin-top:4px;"><i>Measured: No<i><br/>%s' % ', '.join(tmp2)
        
        if ev.comments is not None and ev.comments != '':
            tmp += '<div style="margin-top:4px;"><i>Comments</i><br/>%s</div>' % ev.comments
            
        if len(ev.references.all()) > 0:
            tmp2 = []
            for ref in ev.references.all():
                tmp2.append(ref.get_citation(True))
            tmp += '<div style="margin-top:4px;"><i>References</i><ol style="margin-top:0px"><li style="margin-bottom:0px">%s</li></ol></div>' % '</li><li style="margin-bottom:0px">'.join(tmp2)
            
        txt.append(tmp)
        
    return '<ul style="margin-top:4px;"><li style="margin-top:4px;">%s</li></ul>' % '</li><li style="margin-top:4px;">'.join(txt)
            
def sub_rate_law(species_wid):
    def inner_method(match):
        if match.group(0) == 'Vmax':
            return '<i>V</i><sub>max</sub>'
        if match.group(0)[0:2] == 'Km':
            return '<i>K</i><sub>m%s</sub>' % match.group(0)[2:]
        try:
            obj = SpeciesComponent.objects.get(species__wid=species_wid, wid=match.group(0))
            return '[<a href="%s">%s</a>]' % (obj.get_absolute_url(species), obj.wid)
        except:
            return match.group(0)
    return inner_method