"""
Cyanofactory data model. Contains the database layout.

Author: Jonathan Karr, jkarr@stanford.edu
Affiliation: Covert Lab, Department of Bioengineering, Stanford University
Last updated: 2012-07-17

Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from __future__ import unicode_literals

import itertools
import math
import re
from django.db.models.fields import NullBooleanField
import subprocess
import sys
from datetime import datetime
import StringIO

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqRecord, SeqIO

from django.contrib.auth.models import User
from django.contrib.auth.models import Group
from django.core import validators
from django.core.exceptions import ValidationError, ObjectDoesNotExist
from django.core.urlresolvers import reverse
from django.db.models import F, Model, OneToOneField, CharField, IntegerField, URLField, PositiveIntegerField,\
    FloatField, BooleanField, SlugField, TextField, DateTimeField, options, permalink, SET_NULL, Min, signals
from django.db.models.query import EmptyQuerySet, QuerySet
from django.template import loader, Context
from django.dispatch.dispatcher import receiver
from django.db.models.signals import m2m_changed

from .templatetags.templatetags import set_time_zone
from .cache import Cache
from .history import HistoryForeignKey as ForeignKey
from .history import HistoryManyToManyField as ManyToManyField
from Bio.SeqFeature import SeqFeature, FeatureLocation
from copy import deepcopy
from django.conf import settings

from model_utils.managers import PassThroughManager

from south.modelsinspector import add_introspection_rules
add_introspection_rules([], ["^cyano\.history\.HistoryForeignKey"])
add_introspection_rules([], ["^cyano\.history\.HistoryManyToManyField"])

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
options.DEFAULT_NAMES = options.DEFAULT_NAMES + ('listing', 'concrete_entry_model', 'fieldsets', 'field_list',
                                                 'facet_fields', 'clean', 'validate_unique', 'wid_unique',
                                                 'group_field')

''' BEGIN: validators '''
def validate_dna_sequence(seq):
    validators.RegexValidator(regex=r'^[ACGT]+$', message='Enter a valid DNA sequence consisting of only the letters A, C, G, and T')(seq)

def validate_protein_sequence(seq):
    validators.RegexValidator(regex=r'^[ACDEFGHIKLMNPQRSTVWY]+$', message='Enter a valid Protein sequence consisting of only the letters A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W and Y')(seq)

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
                usedKmMax = max(usedKmMax, 1)
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

    pre = ''
    blocks = []
    posts = []
    pattern = '%s'
    begin = 0
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
    def get_instance(permission):
        if isinstance(permission, basestring):
            permission = Permission.get_by_name(permission)
        elif isinstance(permission, int):
            permission = Permission.get_by_pk(permission)
        elif isinstance(permission, Permission):
            permission = permission
        else:
            raise ValueError("Must be string, int or Permission type")

        if not isinstance(permission, Permission):
            raise ValueError("Invalid permission argument")

        return permission

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
    entry = ForeignKey("Entry", related_name = "user_permissions", auto_created = True)
    user = ForeignKey("UserProfile", related_name = 'permissions')
    allow = ManyToManyField(Permission, verbose_name = 'Allowed Permissions', related_name = 'user_permission_allow')
    deny = ManyToManyField(Permission, verbose_name = 'Denied Permissions', related_name = 'user_permission_deny')

    def __unicode__(self):
        perm_allow = [str(1) if Permission.get_by_pk(x) in self.allow.all() else str(0) for x in range(1, 9)]
        perm_deny = [str(1) if Permission.get_by_pk(x) in self.deny.all() else str(0) for x in range(1, 9)]

        return "[{}, {}]".format(
            "".join(perm_allow), "".join(perm_deny)
        )

    class Meta:
        unique_together = ('entry', 'user')

class GroupPermission(Model):
    entry = ForeignKey("Entry", related_name = "group_permissions", auto_created = True)
    group = ForeignKey("GroupProfile", related_name = 'permissions')
    allow = ManyToManyField(Permission, verbose_name = 'Allowed Permissions', related_name = 'group_permission_allow')
    deny = ManyToManyField(Permission, verbose_name = 'Denied Permissions', related_name = 'group_permission_deny')

    def __unicode__(self):
        perm_allow = [str(1) if Permission.get_by_pk(x) in self.allow.all() else str(0) for x in range(1, 9)]
        perm_deny = [str(1) if Permission.get_by_pk(x) in self.deny.all() else str(0) for x in range(1, 9)]

        return "[{}, {}]".format(
            "".join(perm_allow), "".join(perm_deny)
        )

    class Meta:
        unique_together = ('entry', 'group')

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

    def add_allow_permission(self, entry, permission):
        self._update_permission(entry, permission, allow=True, add=True)

    def delete_allow_permission(self, entry, permission):
        self._update_permission(entry, permission, allow=True, delete=True)

    def add_deny_permission(self, entry, permission):
        self._update_permission(entry, permission, deny=True, add=True)

    def delete_deny_permission(self, entry, permission):
        self._update_permission(entry, permission, deny=True, delete=True)

    def _handle_permission_list(self, entry, perm_list, allow):
        """
        Internal Api used by permission view
        """
        for i, perm in enumerate(perm_list):
            if allow:
                if perm == 0:
                    self.delete_allow_permission(entry, i + 1)
                else:
                    self.add_allow_permission(entry, i + 1)
            else:
                if perm == 0:
                    self.delete_deny_permission(entry, i + 1)
                else:
                    self.add_deny_permission(entry, i + 1)

    def _update_permission(self, permission, perm_get, allow=False, deny=False, add=False, delete=False):
        permission = Permission.get_instance(permission)

        if allow ^ deny or add ^ delete:
            ValueError("Invalid args")

        if add:
            op = "add"
            perm = perm_get()
        else:
            op = "remove"
            try:
                perm = perm_get()
            except ObjectDoesNotExist:
                return

        if allow:
            field = perm.allow
        else:
            field = perm.deny

        getattr(field, op)(permission)

        # Test if both fields are empty now, in that case the permission can be deleted
        if op == "remove" and perm.allow.count() + perm.deny.count() == 0:
            perm.delete()

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

    def _update_permission(self, entry, permission, allow=False, deny=False, add=False, delete=False):
        if allow ^ deny or add ^ delete:
            ValueError("Invalid args")

        if add:
            perm_get = lambda: GroupPermission.objects.get_or_create(entry = entry, group = self)[0]
        else:
            perm_get = lambda: GroupPermission.objects.get(entry = entry, group = self)

        super(GroupProfile, self)._update_permission(permission, perm_get, allow, deny, add, delete)

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
        if allow is None:
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
    force_password_change = BooleanField(default=False)

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

    def _update_permission(self, entry, permission, allow=False, deny=False, add=False, delete=False):
        if allow ^ deny or add ^ delete:
            ValueError("Invalid args")

        if add:
            perm_get = lambda: UserPermission.objects.get_or_create(entry = entry, user = self)[0]
        else:
            perm_get = lambda: UserPermission.objects.get(entry = entry, user = self)

        super(UserProfile, self)._update_permission(permission, perm_get, allow, deny, add, delete)

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

def create_user_profile_signal(sender, instance, created, **kwargs):
    if created:
        UserProfile.objects.create(user=instance)

def password_change_signal(sender, instance, **kwargs):
    try:
        user = User.objects.get(username=instance.username)
        if not user.password == instance.password:
            profile = user.profile
            profile.force_password_change = False
            profile.save()
    except User.DoesNotExist:
        pass

signals.pre_save.connect(password_change_signal, sender=User, dispatch_uid='cyano.models')
signals.post_save.connect(create_user_profile_signal, sender=User, dispatch_uid='cyano.models')

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

    @staticmethod
    def get_by_id(pk):
        model_type_key = "model/model_type/pk/" + str(pk)
        cache_model_type = Cache.try_get(model_type_key, lambda: TableMeta.objects.get(pk = pk))
        return cache_model_type

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
    current = ForeignKey("Entry", verbose_name = "Current version of entry", related_name = 'revisions', editable = False, auto_created = True)
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
    current = ForeignKey("Entry", verbose_name = "Current version of entry points to source of m2m", related_name = 'revisions_m2m', editable = False, auto_created = True)
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

class EntryBasicTextData(EntryData):
    value = TextField(blank=False, default='', verbose_name='Value')

    class Meta:
        verbose_name = "Entry basic text data"
        verbose_name_plural = "Entry basic text data"

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

        if self.vmax_unit.name  == 'U/mg':
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
        if action == "pre_add":
            # target ist species?
            if TableMetaManyToMany.get_by_m2m_model_name(sender._meta.object_name).target_table.model_name == "Species":
                species = model.objects.get(pk=list(pk_set)[0])
                if species.species_components.filter(model_type = instance.model_type, wid = instance.wid).exists():
                    raise ValidationError("Species {}: Wid {} already in use for model_type .{}".format(species.wid, instance, instance._meta.object_name))

        if action == "post_add":
            #print "m2m-add:",TableMetaManyToMany.get_by_m2m_model_name(sender._meta.object_name),instance,pk_set

            table_meta = TableMetaManyToMany.get_by_m2m_model_name(sender._meta.object_name)

            save_list = []

            for pk in pk_set:
                # FIXME: detail_id is the current one
                save_list.append(RevisionManyToMany(current_id = instance.pk, detail_id = instance.detail_id, action = "I", table = table_meta, new_value = pk))

            if len(save_list) > 0:
                #print instance.wid + ": revisioning", len(save_list), "items"
                # Don't update the detail when actually nothing changed for that entry
                RevisionManyToMany.objects.bulk_create(save_list)
        elif action == "post_delete":
            #print "m2m-del:",TableMetaManyToMany.get_by_m2m_model_name(sender._meta.object_name),instance,pk_set

            table_meta = TableMetaManyToMany.get_by_m2m_model_name(sender._meta.object_name)

            save_list = []

            for pk in pk_set:
                # FIXME: detail_id is the current one
                save_list.append(RevisionManyToMany(current_id = instance.pk, detail_id = instance.detail_id, action = "D", table = table_meta, new_value = pk))

            if len(save_list) > 0:
                #print instance.wid + ": deleting", len(save_list), "items"
                # Don't update the detail when actually nothing changed for that entry
                RevisionManyToMany.objects.bulk_create(save_list)

class EntryQuerySet(QuerySet):
    def with_permission(self, permission):
        perm = Permission.get_by_name(permission)

        # Todo

    def for_wid(self, wid, get=True, create=False, creation_status=False):
        if get:
            try:
                if not creation_status:
                    return self.get(wid = wid)

                return self.get(wid = wid), False
            except ObjectDoesNotExist:
                if not create:
                    raise

                if not creation_status:
                    return self.model(wid = wid)

                return self.model(wid = wid), True

        if create:
            raise ValueError("Create only compatible with get")

        return self.filter(wid = wid)

class SpeciesComponentQuerySet(EntryQuerySet):
    def for_species(self, species):
        return self.filter(species__pk = species.pk)

class AbstractEntry(Model):
    objects = PassThroughManager.for_queryset_class(EntryQuerySet)()

    class Meta:
        abstract = True

class Entry(AbstractEntry):
    """Base class for all knowledge base objects.
    
    :Columns:
        * ``model_type``: Specifies the child type of the item (used to find the child table for inheritance)
        * ``wid``: Warehouse Identifier. Unique identifier in the warehouse, used in the URLs e.g.
        * ``name``: Short, human readable name of that entry
        * ``synonyms``: Other names for that entry
        * ``comments``: Detailed notes about this entry
    """
    model_type = ForeignKey(TableMeta)
    wid = SlugField(max_length=150, verbose_name='WID', validators=[validators.validate_slug])
    name = CharField(max_length=255, blank=True, default='', verbose_name='Name')
    synonyms = ManyToManyField(Synonym, blank=True, null=True, related_name='entry', verbose_name='Synonyms')
    comments = TextField(blank=True, default='', verbose_name='Comments')

    created_detail = ForeignKey(RevisionDetail, verbose_name='Entry created revision', related_name='entry_created_detail', editable = False)
    detail = ForeignKey(RevisionDetail, verbose_name='Last Revision entry', related_name='entry_detail', editable = False)

    def __unicode__(self):
        return self.wid

    def natural_key(self):
        return self.wid

    def save(self, revision_detail, *args, **kwargs):
        # Optimized to reduce number of database accesses to a minimum

        from cyano.helpers import slugify

        if self.wid != slugify(self.wid):
            # Slugify the WID
            raise ValidationError("Wid must be slug!")

        if self.pk is not None:
            # can this be optimized? Probably not
            old_item = self._meta.concrete_model.objects.get(pk = self.pk)
            super(Entry, self).save(*args, **kwargs)
        else:
            # New entry (no primary key)
            # The latest entry is not revisioned to save space (and time)
            cache_model_type = TableMeta.get_by_model_name(self._meta.object_name)
            self.created_detail = revision_detail
            self.detail = revision_detail
            self.model_type = cache_model_type

            if self._meta.concrete_model._meta.wid_unique and self._meta.concrete_model.objects.for_wid(self.wid, get = False).exists():
                raise ValidationError("Wid {} for model_type {} is shared and must be unique for all species".format(self.wid, self.model_type.model_name))

            ##print "CREATE: No revision needed for", str(self.wid)
            super(Entry, self).save(*args, **kwargs)
            return

        save_list = []

        fields = self._meta.fields
        # Remove primary keys and some fields that don't need revisioning:
        fields = itertools.ifilter(lambda x: (not x.primary_key) and
                         (not x.name in ["detail", "created_detail"]), fields)

        for field in fields:
            if hasattr(self, field.name + "_id"):
                new_value = getattr(self, field.name + "_id")
            else:
                new_value = getattr(self, field.name)

            if old_item is None:
                old_value = None
            else:
                if hasattr(old_item, field.name + "_id"):
                    old_value = getattr(old_item, field.name + "_id")
                else:
                    old_value = getattr(old_item, field.name)

            if old_value != new_value:
                column = TableMetaColumn.get_by_field(field)

                if new_value is None:
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
            self.detail = revision_detail
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
        detail = self.created_detail
        if is_user_anonymous:
            return '%s<br>%s' % (detail.date.strftime("%Y-%m-%d %H:%M:%S"), detail.reason)
        else:
            user = detail.user.user
            return '<a href="%s">%s %s</a> on %s<br>%s' % (user.get_absolute_url(), user.first_name, user.last_name, detail.date.strftime("%Y-%m-%d %H:%M:%S"), detail.reason)

    def get_as_html_last_updated_user(self, species, is_user_anonymous):
        detail = self.detail
        if is_user_anonymous:
            return '%s<br>%s' % (detail.date.strftime("%Y-%m-%d %H:%M:%S"), detail.reason)
        else:
            user = detail.user.user
            return '<a href="%s">%s %s</a> on %s<br>%s' % (user.get_absolute_url(), user.first_name, user.last_name, detail.date.strftime("%Y-%m-%d %H:%M:%S"), detail.reason)

    def get_as_fasta(self, species):
        raise NotImplementedError("FASTA export not supported for %s" % self._meta.verbose_name_plural)

    def get_fasta_header(self):
        return ">" + self.wid + "|" + self.name

    def get_as_genbank(self, species):
        raise NotImplementedError("GenBank export not supported for %s" % self._meta.verbose_name_plural)
    
    def get_as_sbml(self, species):
        raise NotImplementedError("SBML export not supported for %s" % self._meta.verbose_name_plural)

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
            'id', 'wid', 'name', 'synonyms', 'cross_references', 'comments', 'created_detail', 'detail'
            ]
        listing = ['wid', 'name']
        facet_fields = []
        ordering = ['wid']
        get_latest_by = 'createdDate'
        verbose_name = 'Entry'
        verbose_name_plural = 'Entries'
        wid_unique = False

class CrossReferenceMeta(Model):
    name = CharField(max_length=255, verbose_name = "Identifier for crossreference (DB name or URL)", editable = False)

class AbstractSpeciesComponent(Entry):
    objects = PassThroughManager.for_queryset_class(SpeciesComponentQuerySet)()

    class Meta:
        abstract = True

class SpeciesComponent(AbstractSpeciesComponent):
    '''
    Contains all components that belong to an organism.
    Improves the lookup speed and allows inheritance.
    
    * ``cross_references``: Databases referencing that entry
    * ``publication_references``: Publications referencing that entry
    '''
    species = ManyToManyField("Species", verbose_name = "Organism containing that entry", related_name = "species_components")
    type = ManyToManyField('Type', blank=True, null=True, related_name='members', verbose_name='Type')
    cross_references = ManyToManyField("CrossReference", blank=True, null=True, related_name='cross_referenced_components', verbose_name='Cross references')
    publication_references = ManyToManyField("PublicationReference", blank=True, null=True, related_name='publication_referenced_components', verbose_name='Publications')
    parent = ForeignKey('self', blank=True, null=True, on_delete=SET_NULL, related_name='children', verbose_name='Parent')

    #getters

    def delete(self, species, using=None):
        """
        Delete is referenced counted and only detached when at least one more item is left
        """
        #if self.species.count() > 1:
            # Detach
        #    # Todo revisioning
        self.species.remove(species.pk)
        #else:
            # Delete
        #    # Todo revisioning
        #    super(SpeciesComponent, self).delete(using=using)

    #@permalink
    def get_absolute_url(self, species, history_id = None):
        if history_id is None:
            return "%s/%s/%s/%s" % (settings.ROOT_URL, species.wid, TableMeta.get_by_id(self.model_type_id).model_name, self.wid)
        else:
            return "%s/%s/%s/%s/history/%s" % (settings.ROOT_URL, species.wid, TableMeta.get_by_id(self.model_type_id).model_name, self.wid, history_id)

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
            lambda match: '[' + ', '.join(['<a href="%s">%s</a>' % (reverse('cyano.views.detail', kwargs={'species_wid':species.wid, 'model_type': 'PublicationReference', 'wid': x}), x, ) for x in match.group(0)[1:-1].split(', ')]) + ']',
            txt)

    def get_as_html_publication_references(self, species, is_user_anonymous):
        results = {}
        for r in self.get_all_references():
            key = r.authors + ' ' + r.editors
            results[key] = r.get_citation(species, cross_references=True)

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
            'id', 'wid', 'name', 'synonyms', 'cross_references', 'type',  'comments', 'publication_references',  'created_detail', 'detail',
            ]
        facet_fields = ['type']
        verbose_name = 'Species component'
        verbose_name_plural = 'Species components'
        wid_unique = False

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
            'id', 'wid', 'name', 'synonyms', 'cross_references', 'type',  'comments', 'publication_references', 'created_detail', 'detail'
            ]
        facet_fields = ['type']
        verbose_name = 'Molecule'
        verbose_name_plural = 'Molecules'
        wid_unique = False

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
            'id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'prosthetic_groups', 'chaperones', 'dna_footprint', 'regulatory_rule', 'comments', 'publication_references', 'created_detail', 'detail'
            ]
        facet_fields = ['type', 'chaperones', 'dna_footprint__binding', 'dna_footprint__region']
        verbose_name='Protein'
        verbose_name_plural = 'Proteins'
        wid_unique = False

        def clean(self, obj_data, all_obj_data=None, all_obj_data_by_model=None):
            #regulatory rule
            if obj_data['regulatory_rule'] is not None and obj_data['regulatory_rule']['value'] is not None and obj_data['regulatory_rule']['value'] != '':
                parse_regulatory_rule(obj_data['regulatory_rule']['value'], all_obj_data, obj_data['species'])

            #DNA footprint
            if obj_data['dna_footprint'] is not None:
                chr_lens = []

                if all_obj_data_by_model is None:
                    species_id = Species.objects.values('id').get(wid=obj_data['species'])['id']
                    for obj in Genome.objects.values('length').filter(species__id=species_id):
                        chr_lens.append(obj['length'])
                else:
                    for obj in all_obj_data_by_model['Genome']:
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

class Genome(Molecule):
    #parent pointer
    parent_ptr_molecule = OneToOneField(Molecule, related_name='child_ptr_genome', parent_link=True, verbose_name='Molecule')

    #additional fields
    sequence = TextField(blank=True, default='', verbose_name='Sequence', validators=[validate_dna_sequence])
    length = PositiveIntegerField(verbose_name='Length (nt)')

    #getters
    def get_sequence(self):
        return self.sequence
    
    def get_length(self):
        return self.length

    def get_gc_content(self):
        seq = self.sequence
        return float(seq.count('G') + seq.count('C')) / float(len(seq))

    def get_transcription_units(self):
        all_tu_pks = self.genes.all().values_list("transcription_units", flat=True)
        all_tu = TranscriptionUnit.objects.filter(pk__in = all_tu_pks).distinct()

        return all_tu

    #http://www.owczarzy.net/extinct.htm
    def get_extinction_coefficient(self):
        from cyano.helpers import ExtinctionCoefficient

        seq = self.sequence

        value = 0
        for i in range(len(seq) - 1):
            value += ExtinctionCoefficient.pairwise_dna[seq[i]][seq[i+1]]
        value += ExtinctionCoefficient.pairwise_dna[seq[-1]][seq[0]]
        return value / 2

    def get_pi(self):
        return calculate_nucleic_acid_pi(self.sequence)

    #html formatting
    def get_as_html_sequence(self, species, is_user_anonymous):
        from cyano.helpers import format_sequence_as_html
        return format_sequence_as_html(species, self.sequence, show_protein_seq=True)

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
        from .helpers import shift, overlaps
        from collections import namedtuple
        from cyano.string_template import Loader as StringTemplateLoader

        ntPerSegment = 1e4
        segmentHeight = 27
        geneHeight = 10
        featureHeight = 5
        nSegments = int(math.ceil(self.length / ntPerSegment))
        W = 636
        segmentLeft = 35
        segmentW = W - 4 - segmentLeft
        oneNtW = segmentW / ntPerSegment
        segmentWLast = (self.length - ntPerSegment * (nSegments - 1)) * oneNtW
        chrTop = -12
        feature_draw = [[] for x in range(nSegments)]
        row_offset = []
        gene_template = StringTemplateLoader().load_template("cyano/genome/draw_gene.tmpl")[0]
        feature_template_rect = StringTemplateLoader().load_template("cyano/genome/draw_feature_rect.tmpl")[0]
        feature_template_triangle = StringTemplateLoader().load_template("cyano/genome/draw_feature_triangle.tmpl")[0]
        fake_gene = Gene(model_type=TableMeta.get_by_model_name("Gene"))
        fake_cf = ChromosomeFeature(model_type=TableMeta.get_by_model_name("ChromosomeFeature"))

        #style
        colors = ['3d80b3', '3db34a', 'cd0a0a', 'e78f08', 'b33da6', '0acdcd', '0860e7']
        style = ''
        for i in range(len(colors)):
            style += '.color-%s{fill:#%s; stroke: #%s;}' % (i, colors[i], colors[i], )

        #chromosome
        chrStyle = '\
            .chr text{fill:#222; text-anchor:end; alignment-baseline:middle; font-size:10px}\
            .chr line{stroke:#666; stroke-width:0.5px;}\
        '

        #genes
        geneStyle = '\
            .genes g polygon{stroke-width:1px; fill-opacity:0.5;}\
            .genes g text{text-anchor:middle; alignment-baseline:middle; font-size:8px; fill: #222}\
        '

        genes = StringIO.StringIO()

        GeneTuple = namedtuple("GeneTuple", "name wid coordinate direction length transcription_units transcription_units__wid transcription_units__name")

        genesList = map(GeneTuple._make, self.genes.values_list(
            "name", "wid", "coordinate", "direction", "length",
            "transcription_units", "transcription_units__wid", "transcription_units__name").all())

        nTus = 0
        iTUs = {}
        tus = []

        all_feature_type_pks = [None] + list(self.features.prefetch_related("chromosome_feature").values_list("chromosome_feature__type", flat=True).distinct())

        def draw_gene(gene, tu):
            iSegment = math.floor((gene.coordinate - 1) / ntPerSegment)

            tip_title = gene.name or gene.wid
            tip_text = 'Transcription unit: %s' % (tu or "(None)")

            fake_gene.wid = gene.wid
            gene_abs_url = fake_gene.get_absolute_url(species)

            tip_title = tip_title.replace("'", "\'")
            tip_text = tip_text.replace("'", "\'")

            context_dict = {'title': tip_title,
                            'text': tip_text,
                            'url': gene_abs_url,
                            'color': iTu % len(colors)}

            ret = []

            segments = shift(range(nSegments), int(iSegment))

            w_drawn = 0

            for i in itertools.cycle(segments):
                context_dict.update({'direction': gene.direction})

                if i == iSegment:
                    x1 = segmentLeft + ((gene.coordinate - 1) % ntPerSegment) / ntPerSegment * segmentW
                else:
                    x1 = segmentLeft

                w_needed = gene.length / ntPerSegment * segmentW - w_drawn

                row = i + 1
                if row == nSegments:
                    w_space = segmentLeft + segmentWLast - x1
                else:
                    w_space = segmentLeft + segmentW - x1

                y2 = chrTop + row_offset[i] + geneHeight - 2
                y1 = y2 - geneHeight

                if w_space < w_needed:
                    # Not enough space left on line

                    w = max(1, w_space)
                    w_drawn += w
                    x2 = x1 + w

                    if math.fabs(x2 - x1) > len(gene.wid) * 5:
                        label = gene.wid
                    else:
                        label = ''

                    if gene.direction == "f" or iSegment == i:
                        if gene.direction == "r" and x1 + 5 > x2:
                            arrow_size = x2 - x1
                        else:
                            arrow_size = 5

                        if gene.direction == "f":
                            x1_arrow = x1
                            x2_arrow = x2
                        else:
                            x1_arrow = x1 + arrow_size
                            x2_arrow = x2
                    else:
                        x1_arrow, x2_arrow = x1, x2

                    context_dict.update({'x1': x1,
                                         'x2': x2,
                                         'y1': y1,
                                         'y2': y2,
                                         'x1_arrow': x1_arrow,
                                         'x2_arrow': x2_arrow,
                                         'x_middle': (x1 + x2) / 2,
                                         'y_middle': (y1 + y2) / 2,
                                         'label': label,
                                         'color': iTu % len(colors)
                                         })
                    ret.append(gene_template.render(Context(context_dict)))
                else:
                    w = w_needed
                    x2 = x1 + w
                    if math.fabs(x2 - x1) > len(gene.wid) * 5:
                        label = gene.wid
                    else:
                        label = ''

                    if gene.direction == "f" or iSegment == i:
                        if x2 - 5 < x1:
                            arrow_size = abs(x1 - x2)
                        else:
                            arrow_size = 5

                        if gene.direction == "f":
                            x1_arrow = x1
                            x2_arrow = x2 - arrow_size
                        else:
                            x1_arrow = x1 + arrow_size
                            x2_arrow = x2
                    else:
                        x1_arrow, x2_arrow = x1, x2

                    context_dict.update({'x1': x1,
                                         'x2': x2,
                                         'y1': y1,
                                         'y2': y2,
                                         'x1_arrow': x1_arrow,
                                         'x2_arrow': x2_arrow,
                                         'x_middle': (x1 + x2) / 2,
                                         'y_middle': (y1 + y2) / 2,
                                         'label': label,
                                         })
                    ret.append(gene_template.render(Context(context_dict)))
                    break

            return ret

        #promoters
        promoterStyle = '.promoters rect{fill:#%s; opacity:0.5}' % (colors[0], )
        tfSiteStyle = '.tfSites rect{fill:#%s; opacity:0.5}' % (colors[1], )

        promoters = StringIO.StringIO()
        tfSites = StringIO.StringIO()

        def preprocess_draw_segment(coordinate, length, tip_title, typ_pk, typ, url):
            iSegment = math.floor((coordinate - 1) / ntPerSegment)

            segments = shift(range(nSegments), int(iSegment))

            w_drawn = 0
            done = False
            for i in itertools.cycle(segments):
                if i == iSegment:
                    x = (coordinate - 1) % ntPerSegment
                else:
                    x = 0

                w_needed = length - w_drawn

                row = i + 1
                if row == nSegments:
                    w_space = ntPerSegment - x
                else:
                    w_space = ntPerSegment - x

                y = i
                if w_space < w_needed:
                    # Not enough space left on line
                    w = max(1, w_space)
                    w_drawn += w
                else:
                    w = w_needed
                    done = True

                new_item = [x, x + w, tip_title, typ_pk, typ, url]
                inserted = False
                for row in feature_draw[i]:
                    if any(overlaps(item, new_item, 50) for item in row):
                        continue
                    # Space left -> insert
                    row.append(new_item)
                    inserted = True
                    break

                if not inserted:
                    # Create new row
                    feature_draw[i].append([new_item])

                if done:
                    break

        def draw_segment(i, row):
            ret = []

            for j, subrow in enumerate(row):
                for item in subrow:
                    start, end, tip_title, typ_pk, typ, url = item
                    tip_text = typ

                    tip_title = tip_title.replace("'", "\'")

                    color = "#" + str(colors[all_feature_type_pks.index(typ_pk) % len(colors)])
                    context_dict = {'h': featureHeight,
                                    'title': tip_title,
                                    'text': tip_text,
                                    'url': url,
                                    'color': color,
                                    'opacity': 1}

                    x = segmentLeft + start * oneNtW
                    w = max(3, (end - start) * oneNtW)
                    y = row_offset[i] + j * (featureHeight + 1)

                    context_dict.update({'x': x,
                                         'x_middle': 0,
                                         'y': y,
                                         'w': w})

                    if w <= 3:
                        feature_template = feature_template_triangle
                        context_dict.update({
                            'h': y + featureHeight,
                            'w': x + 3,
                            'x': x - 3,
                            'x_middle': x,
                        })
                    else:
                        feature_template = feature_template_rect

                    ret.append(feature_template.render(Context(context_dict)))

            return ret

        for tu in tus:
            if tu.promoter_35_coordinate is not None:
                tu_coordinate = tu.get_coordinate() + tu.promoter_35_coordinate
                tu_length = tu.promoter_35_length

                tip_title = tu.name or tu.wid
                tip_text = 'Promoter -35 box'
                url = tu.get_absolute_url(species)

                [promoters.write(item) for item in draw_segment(tu_coordinate, tu_length)]

            if tu.promoter_10_coordinate is not None:
                tu_coordinate = tu.get_coordinate() + tu.promoter_10_coordinate
                tu_length = tu.promoter_10_length

                tip_title = tu.name or tu.wid
                tip_text = 'Promoter -10 box'
                url = tu.get_absolute_url(species)

                [promoters.write(item) for item in draw_segment(tu_coordinate, tu_length)]

            for tr in tu.transcriptional_regulations.all():
                if tr.binding_site is not None:
                    tr_coordinate = tr.binding_site.coordinate
                    tr_length = tr.binding_site.length

                    tip_title = tr.name or tr.wid
                    tip_text = 'Transcription factor binding site'
                    url = tr.get_absolute_url(species)

                    [tfSites.write(item) for item in draw_segment(tr_coordinate, tr_length)]

        #features
        featureStyle = '.features rect{fill:#%s;}' % (colors[2], )
        features = StringIO.StringIO()


        feature_values = self.features.all().order_by("coordinate").\
            values("coordinate", "length",
                   "chromosome_feature__name", "chromosome_feature__wid",
                   "chromosome_feature__type", "chromosome_feature__type__name", "chromosome_feature__type__wid")

        for feature in feature_values:
            coordinate = feature["coordinate"]
            length = feature["length"]

            tip_title = feature["chromosome_feature__name"] or feature["chromosome_feature__wid"]

            if feature["chromosome_feature__type"]:
                typ = feature["chromosome_feature__type__name"] or feature["chromosome_feature__type__wid"]
            else:
                typ = None
            fake_cf.wid = feature["chromosome_feature__wid"]
            url = fake_cf.get_absolute_url(species)

            preprocess_draw_segment(coordinate, length, tip_title, feature["chromosome_feature__type"], typ, url)

        for i, row in enumerate(feature_draw):
            if i == 0:
                row_offset.append(chrTop + segmentHeight + 2)
            else:
                row_offset.append(chrTop + segmentHeight + 2 + (len(feature_draw[i - 1]) or 1) * (featureHeight + 1) + row_offset[-1])

        for i, row in enumerate(feature_draw):
            [features.write(item) for item in draw_segment(i, row)]

        chro = StringIO.StringIO()

        for i in range(nSegments):
            x1 = segmentLeft
            x2 = segmentLeft + ((min(self.length, (i+1) * ntPerSegment) - 1) % ntPerSegment) / ntPerSegment * segmentW
            y = chrTop + row_offset[i] + geneHeight
            chro.write('<text x="%s" y="%s">%d</text>' % (segmentLeft - 2, y, i * ntPerSegment + 1))
            chro.write('<line x1="%s" x2="%s" y1="%s" y2="%s"/>' % (x1, x2, y, y))

        for i, gene in enumerate(genesList):
            if gene.transcription_units:
                tu_wid = gene.transcription_units__wid
                tu_name = gene.transcription_units__name

                tu = tu_name or tu_wid

                if tu_wid in iTUs:
                    iTu = iTUs[tu_wid]
                else:
                    iTu = nTus
                    iTUs[tu_wid] = iTu
                    nTus += 1
            else:
                tu = None
                iTu = nTus
                nTus += 1

            [genes.write(item) for item in draw_gene(gene, tu)]

        H = row_offset[-1] + len(feature_draw[-1]) * (featureHeight + 2)

        return '<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="%s" height="%s" viewport="0 0 %s %s"><style>%s%s%s%s%s%s</style><g class="chr">%s</g><g class="genes">%s</g><g class="promoters">%s</g><g class="tfSites">%s</g><g class="features">%s</g></svg>' % (
            W, H, W, H, style, chrStyle, geneStyle, promoterStyle, tfSiteStyle, featureStyle, chro.getvalue(), genes.getvalue(), promoters.getvalue(), tfSites.getvalue(), features.getvalue())

    def get_as_html_structure_local(self, species, is_user_anonymous, start_coordinate = None, end_coordinate = None, highlight_wid = None):
        from .helpers import overlaps
        from cyano.string_template import Loader as StringTemplateLoader

        length = end_coordinate - start_coordinate + 1

        W = 636
        geneHeight = 20
        featureHeight = 10

        geneY = 2
        chrY = geneY + geneHeight + 4
        promoterY = chrY + 1 + 2
        featureY = promoterY

        feature_draw = []
        feature_template_rect = StringTemplateLoader().load_template("cyano/genome/draw_feature_rect.tmpl")[0]
        feature_template_triangle = StringTemplateLoader().load_template("cyano/genome/draw_feature_triangle.tmpl")[0]

        #style
        colors = ['3d80b3', '3db34a', 'cd0a0a', 'e78f08', 'b33da6', '0acdcd', '0860e7']
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

        chro = '<text x="%s" y="%s" style="text-anchor:end;">%s</text><line x1="%s" x2="%s" y1="%s" y2="%s" /><text x="%s" y="%s" style="text-anchor:start;">%s</text>' % (
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

        all_feature_type_pks = [None] + list(self.features.prefetch_related("chromosome_feature").values_list("chromosome_feature__type", flat=True).distinct())

        template = loader.get_template("cyano/genome/draw_gene.html")

        for i, gene in enumerate(genesList):
            if len(gene.transcription_units.all()[:1]) == 1:
                tu = gene.transcription_units.all()[:1][0]
                if iTUs.has_key(tu.wid):
                    iTu = iTUs[tu.wid]
                else:
                    tus.append(tu)
                    iTu = nTus
                    iTUs[tu.wid] = iTu
                    nTus += 1
            else:
                tu = None
                iTu = nTus
                nTus += 1

            x1 = chrL + float(gene.coordinate - start_coordinate) / length * chrW
            x2 = chrL + float(gene.coordinate + gene.length - 1 - start_coordinate) / length * chrW

            if gene.direction == "r" and x1 + 5 > x2:
                arrow_size = x2 - x1
            else:
                arrow_size = 5

            arrow_size = -arrow_size if gene.direction == "f" else arrow_size

            x1 = max(chrL, min(chrR, x1))
            x2 = max(chrL, min(chrR, x2))

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

            if math.fabs(x2 - x1) > len(gene.wid) * 5:
                label = gene.wid
            else:
                label = ''

            tip_title = gene.name or gene.wid

            tip_content = 'Transcription unit: %s' % ("(None)" if tu is None else tu.name or tu.wid)
            tip_title = tip_title.replace("'", "\'")
            tip_content = tip_content.replace("'", "\'")

            gene_abs_url = gene.get_absolute_url(species)

            context_dict = {'title': tip_title,
                            'text': tip_content,
                            'url': gene_abs_url,
                            'color': iTu % len(colors),
                            'x1': x1,
                            'x2': x2,
                            'y1': y1,
                            'y2': y2,
                            'label': label,
                            'direction': gene.direction,
                            'arrow': True,
                            'arrow_size': arrow_size,
                            'fill_opacity': fillOpacity,
                            'stroke_opacity': strokeOpacity,
                            'stroke_width': strokeWidth
                            }

            genes.write(template.render(Context(context_dict)))

        #promoters
        promoterStyle = '.promoters rect{fill:#%s; opacity:0.5}' % (colors[0], )
        tfSiteStyle = '.tfSites rect{fill:#%s; opacity:0.5}' % (colors[1], )

        promoters = StringIO.StringIO()
        tfSites = StringIO.StringIO()

        def draw_segment(wid, coordinate, item_length, tip_title, typ, url):
            if highlight_wid is None or wid in highlight_wid:
                opacity = 1
            else:
                opacity = 0.25

            x = chrL + float(coordinate - start_coordinate) / length * chrW
            w = chrL + float(coordinate + item_length - 1 - start_coordinate) / length * chrW

            x = max(chrL, min(chrR, x))
            w = max(chrL, min(chrR, w)) - x

            tip_title = tip_title.replace("'", "\'")
            tip_text = (typ.name or typ.wid) if typ else ''

            ret = []

            context_dict = {'h': featureHeight,
                            'title': tip_title,
                            'text': tip_text,
                            'url': url,
                            'x': x,
                            'w': w,
                            'x_middle': 0,
                            'coordinate': coordinate,
                            'length': item_length,
                            'opacity': opacity,
                            'color': "#" + str(colors[all_feature_type_pks.index(typ.pk if typ else None) % len(colors)])}

            new_item = [x, x + w]
            inserted = False
            y = 0
            for i, row in enumerate(feature_draw):
                if any(overlaps(item, new_item, 5) for item in row):
                    continue
                # Space left -> insert
                row.append(new_item)
                inserted = True
                y = featureY + i * (featureHeight + 2)
                context_dict.update({
                    'y': y
                })
                break

            if not inserted:
                # Create new row
                y = featureY + len(feature_draw) * (featureHeight + 2)
                context_dict.update({
                    'y': y
                })
                feature_draw.append([new_item])

            if w <= 3:
                feature_template = feature_template_triangle
                context_dict.update({
                    'h': y + featureHeight,
                    'w': x + 3,
                    'x': x - 3,
                    'x_middle': x,
                })
            else:
                feature_template = feature_template_rect

            return feature_template.render(Context(context_dict))

        for tu in tus:
            if tu.promoter_35_coordinate is not None:
                tu_coordinate = tu.get_coordinate() + tu.promoter_35_coordinate
                tu_length = tu.promoter_35_length

                if not (tu_coordinate > end_coordinate or tu_coordinate + tu_length - 1 < start_coordinate):
                    tip_title = tu.name or tu.wid
                    url = tu.get_absolute_url(species)

                    promoters.write(draw_segment(tu.wid, tu_coordinate, tu_length, tip_title, 'Promoter -35 box', url))

            if tu.promoter_10_coordinate is not None:
                tu_coordinate = tu.get_coordinate() + tu.promoter_10_coordinate
                tu_length = tu.promoter_10_length

                if not (tu_coordinate > end_coordinate or tu_coordinate + tu_length - 1 < start_coordinate):
                    tip_title = tu.name or tu.wid
                    url = tu.get_absolute_url(species)

                    promoters.write(draw_segment(tu.wid, tu_coordinate, tu_length, tip_title, 'Promoter -10 box', url))

            for tr in tu.transcriptional_regulations.all():
                if tr.binding_site is not None and not (tr.binding_site.coordinate > end_coordinate or tr.binding_site.coordinate + tr.binding_site.length - 1 < start_coordinate):
                    tip_title = tu.name or tu.wid
                    url = tu.get_absolute_url(species)

                    tfSites.write(draw_segment(tu.wid, tr.binding_site.coordinate, tr.binding_site.length, tip_title, 'Transcription factor binding site', url))

        #features
        featureStyle = '.features rect{fill:#%s;}' % (colors[2], )
        features = StringIO.StringIO()

        #feature_values = self.features.all().order_by("coordinate").\
        #    values("coordinate", "length",
        #           "chromosome_feature__name", "chromosome_feature__wid",
        #           "chromosome_feature__type", "chromosome_feature__type__name", "chromosome_feature__type__wid")

        feature_values = self.features.filter(coordinate__lte=end_coordinate, coordinate__gte=start_coordinate + 1 - F("length")).prefetch_related('chromosome_feature', 'chromosome_feature__type').all()

        for feature in feature_values:
            if feature.coordinate > end_coordinate or feature.coordinate + feature.length - 1 < start_coordinate:
                continue

            tip_title = feature.chromosome_feature.name or feature.chromosome_feature.wid
            url = feature.chromosome_feature.get_absolute_url(species)

            if feature.chromosome_feature.type.all().count() > 0:
                type_ = feature.chromosome_feature.type.all()[0]
            else:
                type_ = None

            features.write(draw_segment(feature.chromosome_feature.wid, feature.coordinate, feature.length, tip_title, type_, url))

        H = 2 + geneHeight + 2 + 4 + 1 * (2 + len(feature_draw) * (featureHeight + 2)) + 2

        return '<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="%s" height="%s" viewport="0 0 %s %s"><style>%s%s%s%s%s%s</style><g class="chr">%s</g><g class="genes">%s</g><g class="promoters">%s</g><g class="tfSites">%s</g><g class="features">%s</g></svg>' % (
            W, H, W, H, style, chrStyle, geneStyle, promoterStyle, tfSiteStyle, featureStyle, chro, genes.getvalue(), promoters.getvalue(), tfSites.getvalue(), features.getvalue())

    def get_as_html_new_structure(self, species, is_user_anonymous):
        return "xxx"

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

    def get_as_fasta(self, species):
        return self.get_fasta_header() + "\r\n" + re.sub(r"(.{70})", r"\1\r\n", self.get_sequence()) + "\r\n"
    
    def get_as_genbank(self, species):
        genbank = StringIO.StringIO()
        genes = Gene.objects.filter(species=species, chromosome_id=self.pk).prefetch_related("cross_references", "protein_monomers")
        record = SeqRecord.SeqRecord(Seq(self.sequence, IUPAC.IUPACUnambiguousDNA()))

        record.description = self.name
        record.name = self.wid
        accession = self.cross_references.filter(source="RefSeq")
        if len(accession) > 0:
            record.annotations["accession"] = accession[0].xid

        record.annotations["date"] = self.detail.date.strftime("%d-%b-%Y").upper()
        record.annotations["source"] = species.name
        record.annotations["organism"] = species.name
        record.annotations["comment"] = self.comments
        
        features = record.features

        source = SeqFeature(FeatureLocation(0, self.length - 1), type="source")
        source.qualifiers["organism"] = species.name

        features += [source]

        for item in genes:
            features += item.get_as_seqfeature(species, record.seq)
        
        SeqIO.write(record, genbank, "genbank")
        
        return genbank.getvalue()

    #meta information
    class Meta:
        concrete_entry_model = True
        fieldsets = [
            ('Type', {'fields': ['model_type']}),
            ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}),
            ('Classification', {'fields': ['type']}),
            ('Sequence', {'fields': [
                {'verbose_name': 'Structure', 'name': 'structure'},
                {'verbose_name': 'New structure', 'name': 'new_structure'},
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
            'id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'sequence', 'length', 'comments', 'publication_references', 'created_detail', 'detail'
            ]
        facet_fields = ['type']
        verbose_name = 'Genome'
        verbose_name_plural = 'Genome'
        wid_unique = False

        def clean(self, obj_data, all_obj_data=None, all_obj_data_by_model=None):
            if obj_data['sequence'] is not None and obj_data['sequence'] != '' and len(obj_data['sequence']) != obj_data['length']:
                raise ValidationError({'length': 'Length of sequence property must match length property'})

class Chromosome(Genome):
    #parent pointer
    parent_ptr_genome = OneToOneField(Genome, related_name='child_ptr_chromosome', parent_link=True, verbose_name='Genome')

    #meta information
    class Meta:
        fieldsets = Genome._meta.fieldsets
        field_list = Genome._meta.field_list
        facet_fields = Genome._meta.facet_fields

        concrete_entry_model = True
        verbose_name = 'Chromosome'
        verbose_name_plural = 'Chromosomes'
        wid_unique = False

class Plasmid(Genome):
    #parent pointer
    parent_ptr_genome = OneToOneField(Genome, related_name='child_ptr_plasmid', parent_link=True, verbose_name='Genome')

    #meta information
    class Meta:
        fieldsets = Genome._meta.fieldsets
        field_list = Genome._meta.field_list
        facet_fields = Genome._meta.facet_fields

        concrete_entry_model = True
        verbose_name = 'Plasmid'
        verbose_name_plural = 'Plasmids'
        wid_unique = False

class ChromosomeFeature(SpeciesComponent):
    #parent pointer
    parent_ptr_species_component = OneToOneField(SpeciesComponent, related_name='child_ptr_chromosome_feature',
                                                 parent_link=True, verbose_name='Species component')

    #additional fields
    #positions = reverse relation

    #meta information
    class Meta:
        concrete_entry_model = True
        fieldsets = [
            ('Type', {'fields': ['model_type']}),
            ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}),
            ('Classification', {'fields': ['type']}),
            {'inline': 'positions'},
            ('Comments', {'fields': ['comments', 'publication_references']}),
            ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'},
                                     {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]}),
        ]
        field_list = [
            'id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'comments', 'publication_references', 'created_detail', 'detail'
        ]
        facet_fields = ['type',] #ToDo 'positions.chromosome', 'positions.direction']
        verbose_name = 'Chromosome feature'
        verbose_name_plural = 'Chromosome features'
        wid_unique = False

        def clean(self, obj_data, all_obj_data=None, all_obj_data_by_model=None):
            if all_obj_data is None:
                chro = Genome.objects.get(species__wid=obj_data['species'], wid=obj_data['genome'])
            else:
                chro = all_obj_data[obj_data['genome']]

            if isinstance(chro, Entry):
                chr_len = chro.length
            else:
                chr_len = chro['length']

            if obj_data['coordinate'] > chr_len:
                raise ValidationError({'coordinate': 'Coordinate must be less then chromosome length.'})
            if obj_data['length'] > chr_len:
                raise ValidationError({'length': 'Length must be less then chromosome length.'})

class FeaturePosition(EntryData):
    """
    Stores position data of a feature
    """
    chromosome_feature = ForeignKey(ChromosomeFeature, related_name = "positions", verbose_name="")
    chromosome = ForeignKey(Genome, related_name='features', verbose_name='Chromosome or Plasmid')
    coordinate = PositiveIntegerField(verbose_name='Coordinate (nt)')
    length = PositiveIntegerField(verbose_name='Length (nt)')
    direction = CharField(max_length=10, choices=CHOICES_DIRECTION, verbose_name='Direction')

    def get_sequence(self, cache=False):
        if cache:
            chromosome_key = "chromosome/%s" % self.chromosome_id
            chromosome = Cache.try_get(chromosome_key, lambda: self.chromosome, 60)
        else:
            chromosome = self.chromosome

        seq = chromosome.sequence[self.coordinate - 1:self.coordinate - 1 + self.length]
        if self.direction == 'r':
            seq = unicode(Seq(seq, IUPAC.unambiguous_dna).reverse_complement())
        return seq

    def get_genes(self):
        genes = []
        for g in self.chromosome.genes.all():
            if (
                        g.coordinate <= self.coordinate <= 1 - g.length + g.coordinate
            ) or (
                        g.coordinate <= self.coordinate + self.length - 1 <= 1 - g.length + g.coordinate
            ) or (
                        g.coordinate >= self.coordinate and g.coordinate + g.length - 1 <= self.coordinate + self.length - 1
            ):
                genes.append(g)
        return genes

    def get_transcription_units(self):
        tus = []
        chr_tus = self.chromosome.get_transcription_units().prefetch_related("genes")

        for tu in chr_tus:
            coordinate = tu.get_coordinate()
            length = tu.get_length()
            if (
                        coordinate <= self.coordinate <= 1 - length + coordinate
            ) or (
                        coordinate <= self.coordinate + self.length - 1 <= 1 - length + coordinate
            ) or (
                        coordinate >= self.coordinate and coordinate + length - 1 <= self.coordinate + self.length - 1
            ):
                tus.append(tu)
        return tus

    #html formatting
    def get_as_html_structure(self, species, is_user_anonymous):
        return self.chromosome.get_as_html_structure(species, is_user_anonymous,
                                                     zoom=1,
                                                     start_coordinate=self.coordinate - 500,
                                                     end_coordinate=self.coordinate + self.length + 500,
                                                     highlight_wid=[self.chromosome_feature.wid])

    def get_as_html_sequence(self, species, is_user_anonymous):
        from cyano.helpers import format_sequence_as_html

        direction = CHOICES_DIRECTION[[x[0] for x in CHOICES_DIRECTION].index(self.direction)][1]

        return '%s: <a href="%s">%s</a>, Coordinate: %s (nt), Length: %s (nt), Direction: %s, Sequence: %s' % (
            self.chromosome.model_type.model_name,
            self.chromosome.get_absolute_url(species), self.chromosome.wid,
            self.coordinate, self.length, direction,
            format_sequence_as_html(species, self.get_sequence(), show_protein_seq=True))

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

    def get_as_fasta(self, species):
        return self.get_fasta_header() + "\r\n" + re.sub(r"(.{70})", r"\1\r\n", self.get_sequence(cache=True)) + "\r\n"

    class Meta:
        verbose_name = "Feature Position"
        verbose_name_plural = "Feature Positions"
        fieldsets = [
            ('Structure', {'fields': [
                {'verbose_name': 'Structure', 'name': 'structure'},
                {'verbose_name': 'Sequence', 'name': 'sequence'},
                {'verbose_name': 'Genes', 'name': 'genes'},
                {'verbose_name': 'Transcription units', 'name': 'transcription_units'},
            ]}),
        ]
        field_list = [
            'id', 'chromosome_feature', 'chromosome', 'coordinate', 'length',
        ]

class Compartment(SpeciesComponent):
    #parent pointer
    parent_ptr_species_component = OneToOneField(SpeciesComponent, related_name='child_ptr_compartment', parent_link=True, verbose_name='Species component')

    #additional fields
    def get_protein_complexes(self, species):
        arr = []
        for obj in ProteinComplex.objects.for_species(species):
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
        for p in self.get_protein_complexes(species):
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
            'id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'comments', 'publication_references', 'created_detail', 'detail'
            ]
        facet_fields = ['type']
        verbose_name='Compartment'
        verbose_name_plural = 'Compartments'
        wid_unique = False

class Gene(Molecule):
    #parent pointer
    parent_ptr_molecule = OneToOneField(Molecule, related_name='child_ptr_gene', parent_link=True, verbose_name='Molecule')

    #additional fields
    symbol = CharField(max_length=255, blank=True, default='', verbose_name='Symbol')
    chromosome = ForeignKey(Genome, related_name='genes', verbose_name='Chromosome or Plasmid')
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
    def get_sequence(self, cache=False):
        if cache:
            chromosome_key = "chromosome/%s" % self.chromosome_id
            chromosome = Cache.try_get(chromosome_key, lambda: self.chromosome, 60)
        else:
            chromosome = self.chromosome
        seq = chromosome.sequence[self.coordinate - 1:self.coordinate - 1 + self.length]

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
            + Metabolite.objects.for_species(self.species).for_wid('AMP').get_empirical_formula() * seq.count('A') \
            + Metabolite.objects.for_species(self.species).for_wid('GMP').get_empirical_formula() * seq.count('C') \
            + Metabolite.objects.for_species(self.species).for_wid('GMP').get_empirical_formula() * seq.count('G') \
            + Metabolite.objects.for_species(self.species).for_wid('UMP').get_empirical_formula() * seq.count('T') \
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

        return '%s: <a href="%s">%s</a>, Coordinate: %s (nt), Length: %s (nt), Direction: %s, G/C content: %.1f%%, Sequence: %s' % (
            self.chromosome.model_type.model_name,
            self.chromosome.get_absolute_url(species), self.chromosome.wid,
            self.coordinate, self.length, direction,
            self.get_gc_content() * 100,
            format_sequence_as_html(species, self.get_sequence(), show_protein_seq=True))

    def get_as_fasta(self, species):
        return self.get_fasta_header() + "\r\n" + re.sub(r"(.{70})", r"\1\r\n", self.get_sequence(cache=True)) + "\r\n"
    
    def get_as_seqfeature(self, species, sequence):
        gene = SeqFeature(FeatureLocation(self.coordinate - 1,
                                          self.coordinate + self.length - (0 if self.direction == "f" else 1),
                                          strand=1 if self.direction == "f" else -1), type="gene")
        gene.qualifiers["locus_tag"] = [self.wid]
        if self.name:
            gene.qualifiers["gene"] = [self.name]

        db_xref = []
        EC_number = []
        for reference in self.cross_references.all():
            if reference.source == "EC":
                EC_number.append(reference.xid)
            db_xref.append(":".join([reference.source, reference.xid]))
        
        if len(db_xref) > 0:
            gene.qualifiers["db_xref"] = db_xref
            
        cds = deepcopy(gene)
        if self.comments:
            cds.qualifiers["note"] = [self.comments]
        if len(EC_number) > 0:
            gene.qualifiers["EC_number"] = EC_number
        cds.qualifiers["transl_table"] = [species.genetic_code]
        cds.qualifiers["codon_start"] = [1]

        typ = list(self.type.all()[:1])
        cds.type = typ[0].wid
        if cds.type == "mRNA":
            cds.type = "CDS"
            s = gene.extract(sequence)
            cds.qualifiers["translation"] = [s.translate(table=species.genetic_code)[:-1]]

        #monomer = ProteinMonomer.objects.values_list('wid', 'name').filter(species = species).get(gene_id = self.pk)
        monomer = list(self.protein_monomers.all()[:1])
        if len(monomer) > 0:
            cds.qualifiers["protein_id"] = monomer[0].wid
            cds.qualifiers["product"] = monomer[0].name

        return gene, cds
        
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
        wid_unique = False

        #chromosome coordinate, length
        def clean(self, obj_data, all_obj_data=None, all_obj_data_by_model=None):
            if all_obj_data is None:
                chro = Genome.objects.get(species__wid=obj_data['species'], wid=obj_data['chromosome'])
            else:
                chro = all_obj_data[obj_data['chromosome']]

            if isinstance(chro, Entry):
                chr_len = chro.length
            else:
                chr_len = chro['length']

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
        iupac_name, traditional_name, log_p, log_d, charge, pi, volume = sys.stdout[1]('\t')[1:]
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
            'created_detail', 'detail'
            ]
        facet_fields = ['type', 'charge', 'is_hydrophobic']
        verbose_name='Metabolite'
        verbose_name_plural = 'Metabolites'
        wid_unique = False

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
            'created_detail', 'detail'
            ]
        facet_fields = ['type']
        verbose_name='Note'
        verbose_name_plural = 'Notes'
        wid_unique = False

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
            'created_detail', 'detail'
            ]
        facet_fields = ['type', 'reactions', 'molecules', 'state', 'process']
        verbose_name='Misc. parameter'
        verbose_name_plural = 'Misc. parameters'
        wid_unique = False

class Pathway(SpeciesComponent):
    #parent pointer
    parent_ptr_species_component = OneToOneField(SpeciesComponent, related_name='child_ptr_pathway', parent_link=True, verbose_name='Species component')

    #helpers
    @staticmethod
    def add_boehringer_pathway(species, revdetail):
        wid = "Boehringer"
        try:
            x = Pathway.objects.for_species(species).for_wid(wid)
            return
        except ObjectDoesNotExist:
            # Create new boehringer (if it was never created yet) or assign to one
            x, created = Pathway.objects.for_wid(wid, create=True, creation_status=True)
            if created:
                typ = Type.objects.for_wid("Main-Pathway", create=True)
                typ.save(revdetail)
                typ.species.add(species)
                x.name = "Biochemical Pathways"
                x.save(revdetail)
                x.type.add(typ)

        x.species.add(species)

    @staticmethod
    def add_kegg_pathway(species, revdetail):
        from kegg.models import Map

        crs = CrossReference.objects.filter(cross_referenced_components__species = species, source = "EC").values_list('xid', flat=True)
        maps = Map.objects.filter(ec_numbers__name__in=crs).distinct()

        for map_ in maps:
            x = Pathway.objects.for_wid(map_.name, create=True)
            x.name = map_.title
            x.save(revdetail)
            if map_.overview:
                typ = Type.objects.for_wid("Main-Pathway", create=True)
                typ.save(revdetail)
                typ.species.add(species)
                x.type.add(typ)
            x.species.add(species)

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
        species_ecs = CrossReference.objects.filter(cross_referenced_components__species = species, source = "EC").values_list('xid', flat=True)
        enzymes = bmodels.Enzyme.objects.filter(ec__in = species_ecs).order_by('title')

        species_metabolites = Metabolite.objects.for_species(species).values_list('name', flat=True)
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

        # Wordaround broken PIL installation
        # see: https://code.djangoproject.com/ticket/6054
        try:
            from PIL import Image
        except ImportError:
            import Image

        from xml.etree.ElementTree import ElementTree, Element

        # All Pathways of the organism
        species_ecs = CrossReference.objects.filter(cross_referenced_components__species = species, source = "EC").values_list('xid', flat=True)

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
                color_component = "red"
                fill_opacity = "0.0"
                fill_color = "green" # not displayed with 0.0

                # Found objects in the organism show a background color

                # Repoint pathway links to our versions
                if "show_pathway?" in url:
                    pathway_name = url[url.index("?") + 1:]
                    color_component = "blue"

                    try:
                        pw_obj = Pathway.objects.for_species(species).for_wid(pathway_name)
                        elem.set("xlink:href", pw_obj.get_absolute_url(species))
                        fill_opacity = "0.2"
                        fill_color = "blue"
                    except ObjectDoesNotExist:
                        elem.set("xlink:href", url)
                        elem.set("target", "_blank")
                else:
                    elem.set("xlink:href", url)
                    elem.set("target", "_blank")

                elem.set("xlink:title", title)
                root.remove(area)
                elem.append(area)

                ecs = self.uniqify(self.extract_ecs(title))

                if len(ecs) > 0:
                    color_component = "green"
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

            template = loader.get_template("cyano/pathway/sidebar.html")

            out = StringIO.StringIO()
            out.write(template.render(Context()))

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
            'id', 'name', 'synonyms', 'cross_references',
            'comments',
            'publication_references',
            'created_detail', 'detail'
            ]
        group_field = 'type'
        facet_fields = ['type']
        verbose_name='Pathway'
        verbose_name_plural = 'Pathways'
        wid_unique = True

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
        for complexe in self.formed_complexes.all():
            results.append('<a href="%s">%s</a><br/>%s' % (complexe.get_absolute_url(species), complexe.name, complexe.get_as_html_biosynthesis(is_user_anonymous)))
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
            'created_detail', 'detail'
            ]
        facet_fields = ['type']
        verbose_name='Process'
        verbose_name_plural = 'Processes'
        wid_unique = False

        @staticmethod
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
        from cyano.helpers import getEntry
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
            return Compartment.objects.for_species(self.species).for_wid('m')

        if len(localizations) == 3 and len(set(['c', 'm', 'e']) & set([x.wid for x in localizations])) == 3:
            return Compartment.objects.for_species(self.species).for_wid('m')

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
            'created_detail', 'detail'
            ]
        facet_fields = ['type', 'dna_footprint__binding', 'dna_footprint__region', 'formation_process', 'chaperones']
        verbose_name='Protein complex'
        verbose_name_plural = 'Protein complexes'
        wid_unique = False

        def clean(self, obj_data, all_obj_data=None, all_obj_data_by_model=None):
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
    def get_sequence(self, species, cache=False):
        return unicode(Seq(self.gene.get_sequence(cache=cache), IUPAC.unambiguous_dna).translate(table=species.genetic_code))

    def get_length(self, species):
        return len(self.get_sequence(species))

    def get_neg_aa(self, species):
        seq = self.get_sequence(species)
        return seq.count('E') + seq.count('D')

    def get_pos_aa(self, species):
        seq = self.get_sequence(species)
        return seq.count('R') + seq.count('H') + seq.count('K')

    def get_n_terminal_aa(self, species):
        return self.get_sequence(species)[0]

    def get_empirical_formula(self, species):
        from cyano.helpers import EmpiricalFormula

        seq = self.get_sequence(species)
        return \
            + Metabolite.objects.for_species(self.species).for_wid('ALA').get_empirical_formula() * seq.count('A') \
            + Metabolite.objects.for_species(self.species).for_wid('ARG').get_empirical_formula() * seq.count('R') \
            + Metabolite.objects.for_species(self.species).for_wid('ASN').get_empirical_formula() * seq.count('N') \
            + Metabolite.objects.for_species(self.species).for_wid('ASP').get_empirical_formula() * seq.count('D') \
            + Metabolite.objects.for_species(self.species).for_wid('CYS').get_empirical_formula() * seq.count('C') \
            + Metabolite.objects.for_species(self.species).for_wid('GLU').get_empirical_formula() * seq.count('E') \
            + Metabolite.objects.for_species(self.species).for_wid('GLN').get_empirical_formula() * seq.count('Q') \
            + Metabolite.objects.for_species(self.species).for_wid('GLY').get_empirical_formula() * seq.count('G') \
            + Metabolite.objects.for_species(self.species).for_wid('HIS').get_empirical_formula() * seq.count('H') \
            + Metabolite.objects.for_species(self.species).for_wid('ILE').get_empirical_formula() * seq.count('I') \
            + Metabolite.objects.for_species(self.species).for_wid('LEU').get_empirical_formula() * seq.count('L') \
            + Metabolite.objects.for_species(self.species).for_wid('LYS').get_empirical_formula() * seq.count('K') \
            + Metabolite.objects.for_species(self.species).for_wid('MET').get_empirical_formula() * (seq.count('M')-1) + Metabolite.objects.for_species(self.species).for_wid('FMET').get_empirical_formula() * (1) \
            + Metabolite.objects.for_species(self.species).for_wid('PHE').get_empirical_formula() * seq.count('F') \
            + Metabolite.objects.for_species(self.species).for_wid('PRO').get_empirical_formula() * seq.count('P') \
            + Metabolite.objects.for_species(self.species).for_wid('SER').get_empirical_formula() * seq.count('S') \
            + Metabolite.objects.for_species(self.species).for_wid('THR').get_empirical_formula() * seq.count('T') \
            + Metabolite.objects.for_species(self.species).for_wid('TRP').get_empirical_formula() * seq.count('W') \
            + Metabolite.objects.for_species(self.species).for_wid('TYR').get_empirical_formula() * seq.count('Y') \
            + Metabolite.objects.for_species(self.species).for_wid('VAL').get_empirical_formula() * seq.count('V') \
            - EmpiricalFormula(H=2, O=1) * (len(seq)-1)

    def get_pi(self, species):
        seq = self.get_sequence(species)

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
        return value

    def get_half_life(self):
        return 20 * 60

    #http://ca.expasy.org/tools/protparam-doc.html
    def get_instability(self, species):
        from cyano.helpers import DipeptideInstabilityWeight

        seq = self.get_sequence(species)
        value = 0.
        for i in range(len(seq)-1):
            if seq[i] != '*' and seq[i+1] != '*':
                value += DipeptideInstabilityWeight.value[seq[i]][seq[i+1]]
        return 10. / float(len(seq)) * value

    #http://ca.expasy.org/tools/protparam-doc.html
    def get_is_stable(self, species):
        return self.get_instability(species) < 40.

    #http://ca.expasy.org/tools/protparam-doc.html
    def get_aliphatic(self, species):
        seq = self.get_sequence(species)
        return 100. * (
            + 1.0 * float(seq.count('A'))
            + 2.9 * float(seq.count('V'))
            + 3.9 * float(seq.count('I'))
            + 3.9 * float(seq.count('L'))
                      ) / float(len(seq))

    #http://ca.expasy.org/tools/protparam-doc.html
    def get_gravy(self, species):
        seq = self.get_sequence(species)
        return \
            (
                + 1.8 * float(seq.count('A'))
                - 4.5 * float(seq.count('R'))
                - 3.5 * float(seq.count('N'))
                - 3.5 * float(seq.count('D'))
                + 2.5 * float(seq.count('C'))
                - 3.5 * float(seq.count('Q'))
                - 3.5 * float(seq.count('E'))
                - 0.4 * float(seq.count('G'))
                - 3.2 * float(seq.count('H'))
                + 4.5 * float(seq.count('I'))
                + 3.8 * float(seq.count('L'))
                - 3.9 * float(seq.count('K'))
                + 1.9 * float(seq.count('M'))
                + 2.8 * float(seq.count('F'))
                - 1.6 * float(seq.count('P'))
                - 0.8 * float(seq.count('S'))
                - 0.7 * float(seq.count('T'))
                - 0.9 * float(seq.count('W'))
                - 1.3 * float(seq.count('Y'))
                + 4.2 * float(seq.count('V'))
            ) / float(len(seq))

    #Source: http://ca.expasy.org/tools/protparam-doc.html
    def get_extinction_coefficient(self, species):
        seq = self.get_sequence(species)
        return \
            + seq.count('W') * 5500 \
            + seq.count('Y') * 1490 \
            + seq.count('C') * 125

    #html formatting
    def get_as_html_sequence(self, species, is_user_anonymous):
        from cyano.helpers import format_sequence_as_html
        return format_sequence_as_html(species, self.get_sequence(species))

    def get_as_html_signal_sequence(self, species, is_user_anonymous):
        ss = self.signal_sequence
        if ss is None:
            return
        return format_with_evidence(obj = ss, txt = 'Type: %s, Location: %s, Length: %s (nt)' % (ss.type, ss.location, ss.length))

    def get_as_html_disulfide_bonds(self, species, is_user_anonymous):
        results = []
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
    
    def get_as_fasta(self, species):
        return self.get_fasta_header() + "\r\n" + re.sub(r"(.{70})", r"\1\r\n", self.get_sequence(species, cache=True)) + "\r\n"

    def get_as_genbank(self, species):
        genbank = StringIO.StringIO()

        record = SeqRecord.SeqRecord(Seq(self.get_sequence(species)[:-1], IUPAC.IUPACProtein()))
        record.description = self.name
        record.name = self.wid
        accession = self.cross_references.filter(source="RefSeq")
        if len(accession) > 0:
            record.annotations["accession"] = accession[0].xid

        record.annotations["date"] = self.detail.date.strftime("%d-%b-%Y").upper()
        record.annotations["source"] = species.name
        record.annotations["organism"] = species.name
        record.annotations["comment"] = self.comments

        source = SeqFeature(FeatureLocation(0, self.get_length(species) - 1), type="source")
        source.qualifiers["organism"] = species.name

        record.features += [source]
        record.features += self.get_as_seqfeature(species)

        SeqIO.write(record, genbank, "genbank")

        return genbank.getvalue()

    def get_as_seqfeature(self, species):
        gene = SeqFeature(FeatureLocation(0, self.get_length(species) - 1), type="Protein")
        gene.qualifiers["locus_tag"] = [self.wid]
        if self.name:
            gene.qualifiers["gene"] = [self.name]

        db_xref = []
        EC_number = []
        for reference in self.cross_references.all():
            if reference.source == "EC":
                EC_number.append(reference.xid)
            db_xref.append(":".join([reference.source, reference.xid]))

        if len(db_xref) > 0:
            gene.qualifiers["db_xref"] = db_xref

        cds = deepcopy(gene)
        if self.comments:
            cds.qualifiers["note"] = [self.comments]
        if len(EC_number) > 0:
            gene.qualifiers["EC_number"] = EC_number
        cds.qualifiers["transl_table"] = [species.genetic_code]
        cds.qualifiers["codon_start"] = [1]

        cds.type = "CDS"
        #if cds.type == "mRNA":
        #    cds.type = "CDS"
            #s = self.get_sequence(species)
            #cds.qualifiers["translation"] = [s.translate(table=species.genetic_code)[:-1]]

        return gene, cds

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
            'created_detail', 'detail'
            ]
        facet_fields = ['type', 'is_n_terminal_methionine_cleaved__value', 'signal_sequence__type', 'signal_sequence__location', 'dna_footprint__binding', 'dna_footprint__region', 'localization', 'chaperones']
        verbose_name='Protein monomer'
        verbose_name_plural = 'Protein monomers'
        wid_unique = False

        def clean(self, obj_data, all_obj_data=None, all_obj_data_by_model=None):
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
            law = '<i>v</i> = %s;' % re.sub(r'([a-z0-9_]+)', sub_rate_law(self.species), law, flags=re.I)
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
            law = '<i>v</i> = %s;' % re.sub(r'([a-z0-9_]+)', sub_rate_law(self.species), law, flags=re.I)
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
            'created_detail', 'detail'
            ]
        facet_fields = ['type', 'direction', 'enzyme__protein', 'coenzymes__metabolite', 'is_spontaneous', 'pathways', 'processes', 'states']
        verbose_name='Reaction'
        verbose_name_plural = 'Reactions'
        wid_unique = False

        def clean(self, obj_data, all_obj_data=None, all_obj_data_by_model=None):
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
    def get_absolute_url(self, species = None, history_id = None):
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
            'id', 'wid', 'name', 'synonyms', 'publication_references', 'cross_references', 'genetic_code',
            'comments'
            ]
        facet_fields = []
        verbose_name = 'Species'
        verbose_name_plural = 'Species'
        wid_unique = True

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
            'created_detail', 'detail'
            ]
        facet_fields = ['type']
        verbose_name='State'
        verbose_name_plural = 'States'
        wid_unique = False

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
            'created_detail', 'detail'
            ]
        facet_fields = ['type', 'value__units']
        verbose_name='Stimulus'
        verbose_name_plural = 'Stimuli'
        wid_unique = False

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
    def get_chromosome(self, cache=False):
        chro = list(set([g.chromosome_id for g in self.genes.all()]))
        if len(chro) != 1:
            raise

        if cache:
            chromosome_key = "chromosome/%s" % chro[0]
            chromosome = Cache.try_get(chromosome_key, lambda: Genome.objects.get(pk=chro[0]), 60)
        else:
            chromosome = Genome.objects.get(pk=chro[0])
        return chromosome

    def get_coordinate(self):
        # Callee can use prefetch_related for speedup
        return min(map(lambda gene: gene.coordinate, self.genes.all()))

    def get_length(self):
        # Callee can use prefetch_related for speedup
        return max(map(lambda gene: gene.coordinate + gene.length - 1 - self.get_coordinate() + 1, self.genes.all()))

    def get_direction(self):
        direction = list(set([g[0] for g in self.genes.values_list('direction').all()]))
        if len(direction) != 1:
            raise
        return direction[0]

    def get_sequence(self, cache=False):
        seq = self.get_chromosome(cache=cache).sequence[self.get_coordinate() - 1:self.get_coordinate() - 1 + self.get_length()]
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

        return '%s: <a href="%s">%s</a>, Coordinate: %s, Length: %s, Direction: %s, Sequence: %s' % (
            self.get_chromosome().model_type.model_name,
            self.get_chromosome().get_absolute_url(species), self.get_chromosome().wid,
            self.get_coordinate(), self.get_length(), direction,
            format_sequence_as_html(species, self.get_sequence(), show_protein_seq=True))

    def get_as_html_transcriptional_regulations(self, species, is_user_anonymous):
        results = []
        for r in self.transcriptional_regulations.all():
            results.append('<a href="%s">%s</a>: <a href="%s">%s</a>' % (r.get_absolute_url(species), r.wid, r.transcription_factor.get_absolute_url(species), r.transcription_factor.wid))
        return format_list_html(results)
    
    def get_as_fasta(self, species):
        return self.get_fasta_header() + "\r\n" + re.sub(r"(.{70})", r"\1\r\n", self.get_sequence(cache=True)) + "\r\n"

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
            'created_detail', 'detail'
            ]
        facet_fields = ['type']
        verbose_name='Transcription unit'
        verbose_name_plural = 'Transcription units'
        wid_unique = False

        def clean(self, obj_data, all_obj_data=None, all_obj_data_by_model=None):
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

        chro = self.transcription_unit.get_chromosome()

        structure = chr.get_as_html_structure(species, is_user_anonymous,
                zoom = 1,
                start_coordinate = bs.coordinate - 500,
                end_coordinate = bs.coordinate + bs.length + 500,
                highlight_wid = [self.wid])

        txt = '%s<br/>%s: <a href="%s">%s</a>, Coordinate: %s (nt), Length: %s (nt), Direction: %s, Sequence: %s' % (
            chro.model_type.model_name,
            structure, chro.get_absolute_url(species), chro.wid,
            bs.coordinate, bs.length, direction,
            format_sequence_as_html(species, self.get_binding_site_sequence(), show_protein_seq=True))

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
            'created_detail', 'detail'
            ]
        facet_fields = ['type', 'transcription_unit', 'transcription_factor']
        verbose_name='Transcriptional regulation'
        verbose_name_plural = 'Transcriptional regulation'
        wid_unique = False

        def clean(self, obj_data, all_obj_data=None, all_obj_data_by_model=None):
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
                chro = Genome.objects.get(species__wid=obj_data['species'], wid=chr_wid)
            else:
                chro = all_obj_data[chr_wid]
            if isinstance(chro, Entry):
                chr_len = chro.length
            else:
                chr_len = chro['length']

            #error check binding site coordinate, length
            if obj_data['binding_site'] is not None:
                if obj_data['binding_site']['coordinate'] > chr_len:
                    raise ValidationError({'binding_site': 'Coordinate must be less then chromosome length.'})
                if obj_data['binding_site']['length'] > chr_len:
                    raise ValidationError({'binding_site': 'Length must be less then chromosome length.'})

class Type(SpeciesComponent):
    #parent pointer
    parent_ptr_species_component = OneToOneField(SpeciesComponent, related_name='child_ptr_type', parent_link=True, verbose_name='Species component')

    #getters
    def get_all_members(self, species):
        members = []
        for m in self.members.for_species(species):
            members.append(m)
        for c in self.children.for_species(species):
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
        for c in self.children.for_species(species):
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
            'created_detail', 'detail'
            ]
        facet_fields = ['type', 'parent']
        verbose_name='Type'
        verbose_name_plural = 'Types'
        wid_unique = True

class CrossReference(EntryData):
    xid = CharField(max_length=255, verbose_name='External ID')
    source = CharField(max_length=20, choices=CHOICES_CROSS_REFERENCE_SOURCES, verbose_name='Source')

    class Meta:
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
        cite_type = None
        props = []
        if self.type.all()[0].wid == 'article':
            cite_type = 'ARTICLE'
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
            cite_type = 'BOOK'
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
            cite_type = 'THESIS'
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
            cite_type = 'MISC'
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
        return '@%s{%s%s\n}' % (cite_type, self.wid, ''.join(tmp))

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
            'created_detail', 'detail'
            ]
        facet_fields = ['type', 'year', 'publication']
        verbose_name='Publication Reference'
        verbose_name_plural = 'Publication References'
        wid_unique = True

class MassSpectrometryJob(SpeciesComponent):
    #parent pointer
    parent_ptr_species_component = OneToOneField(SpeciesComponent, related_name='child_ptr_mass_spectrometry', parent_link=True, verbose_name='Species component')

    #additional fields
    #

    #getters

    #html formatting
    def get_as_html_target_peptide(self, species, is_user_anonymous):
        results = []
        try:
            target_type = Type.objects.for_wid("Target-Peptide")
        except ObjectDoesNotExist:
            return ""

        for g in self.children.for_species(species).filter(type=target_type).order_by("wid"):
            results.append('<a href="%s">%s</a>' % (g.get_absolute_url(species), g.wid))
        return format_list_html(results, comma_separated=True)

    def get_as_html_decoy_peptide(self, species, is_user_anonymous):
        results = []
        try:
            decoy_type = Type.objects.for_wid("Decoy-Peptide")
        except ObjectDoesNotExist:
            return ""

        for g in self.children.for_species(species).filter(type=decoy_type).order_by("wid"):
            results.append('<a href="%s">%s</a>' % (g.get_absolute_url(species), g.wid))
        return format_list_html(results, comma_separated=True)

    def get_as_html_related_proteins(self, species, is_user_anonymous):
        results = []

        for g in self.children.for_species(species).filter(model_type=TableMeta.get_by_model_name("MassSpectrometryProtein")).order_by("wid"):
            results.append('<a href="%s">%s</a>' % (g.get_absolute_url(species), g.wid))
        return format_list_html(results, comma_separated=True)

    #meta information
    class Meta:
        concrete_entry_model = True
        fieldsets = [
            ('Type', {'fields': ['model_type']}),
            ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}),
            ('Classification', {'fields': ['type']}),
            ('Mass Spectrometry', {'fields': [
                {'verbose_name': 'Target Peptides', 'name': 'target_peptide'},
                {'verbose_name': 'Decoy Peptides', 'name': 'decoy_peptide'},
                {'verbose_name': 'Related Proteins', 'name': 'related_proteins'}
            ]}),
            ('Comments', {'fields': ['comments', 'publication_references']}),
            ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]}),
        ]
        field_list = [
            'id', 'wid', 'name', 'synonyms', 'cross_references', 'type',  'comments', 'publication_references',  'created_detail', 'detail',
        ]
        facet_fields = ['type']
        verbose_name = 'Mass Spectrometry Job'
        verbose_name_plural = 'Mass Spectrometry Jobs'
        wid_unique = False

class Peptide(Protein):
    #parent pointer
    parent_ptr_species_component = OneToOneField(SpeciesComponent, related_name='child_ptr_peptide', parent_link=True, verbose_name='Species component')

    #additional fields
    sequence = TextField(blank=True, default='', verbose_name='Sequence', validators=[validate_protein_sequence])
    length = PositiveIntegerField(verbose_name='Length (nt)')
    proteotypic = NullBooleanField(null=True, verbose_name='Proteotypic')
    charge = IntegerField(verbose_name='Charge')
    mass = FloatField(verbose_name='m/z')
    zscore = FloatField(verbose_name='z-score')
    retention_time = FloatField(verbose_name='Retention Time')
    proteins = ManyToManyField(EntryBasicTextData, verbose_name='Proteins belonging to the Peptide', related_name='peptides')

    # Matched Proteins -> FK Protein?
    # Target or Decoy via type

    # SeqPTM via ChromosomeFeature

    #getters

    #html formatting
    def get_as_html_sequence(self, species, is_user_anonymous):
        from cyano.helpers import format_sequence_as_html
        return format_sequence_as_html(species, self.sequence, show_protein_seq=True)

    def get_as_html_matched_proteins(self, species, is_user_anonymous):
        results = []

        for matched_protein in self.proteins.all():
            try:
                res = Protein.objects.for_species(species).get(
                    parent=MassSpectrometryJob.objects.all()[0], wid__startswith=matched_protein.value + "_")
                results.append('<a href="%s">%s</a>' % (res.get_absolute_url(species), res.wid))
            except ObjectDoesNotExist:
                results.append(matched_protein.value)

        return format_list_html(results, comma_separated=True)

    #meta information
    class Meta:
        concrete_entry_model = True
        fieldsets = [
            ('Type', {'fields': ['model_type']}),
            ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}),
            ('Classification', {'fields': ['type']}),
            ('Related', {'fields': [
                'parent',
                {'verbose_name': 'Matched Proteins', 'name': 'matched_proteins'},
            ]}),
            ('Structure', {'fields': [
                'prosthetic_groups', 'chaperones', 'dna_footprint',
                {'verbose_name': 'Sequence', 'name': 'sequence'},
            ]}),
            ('Regulation', {'fields': ['regulatory_rule']}),
            ('Function', {'fields': [
                {'verbose_name': 'Enzyme', 'name': 'enzyme_participants'},
                {'verbose_name': 'Transcriptional regulation', 'name': 'transcriptional_regulations'},
                {'verbose_name': 'Protein folding substrates', 'name': 'chaperone_substrates'},
                {'verbose_name': 'Reaction participant', 'name':'reaction_stoichiometry_participants'},
                {'verbose_name': 'Complex subunit', 'name':'protein_complex_biosythesis_participants'},
                ]}),
            ('Statistics', {'fields': [
                'proteotypic',
                'charge',
                'mass',
                'zscore',
                'retention_time',
            ]}),
            ('Parameters', {'fields': ['parameters']}),
            ('Comments', {'fields': ['comments', 'publication_references']}),
            ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]}),
            ]
        field_list = [
            'id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'prosthetic_groups', 'chaperones', 'dna_footprint', 'regulatory_rule', 'comments', 'publication_references', 'created_detail', 'detail'
            ]
        facet_fields = ['type', 'chaperones', 'dna_footprint__binding', 'dna_footprint__region']
        verbose_name = 'Peptide'
        verbose_name_plural = 'Peptides'
        wid_unique = False


class MassSpectrometryProtein(Protein):
    parent_ptr_species_component = OneToOneField(SpeciesComponent, related_name='child_ptr_ms_protein', parent_link=True, verbose_name='Species component')

    score = FloatField(verbose_name="Protein Score")
    coverage = FloatField(verbose_name="% Coverage")
    sequence = TextField(verbose_name='Sequence', validators=[validate_protein_sequence])
    length = PositiveIntegerField(verbose_name='Length (nt)')
    ambiguous = ManyToManyField(EntryBasicTextData, verbose_name='Ambiguous Proteins', related_name='ambiguous')
    sub = ManyToManyField(EntryBasicTextData, verbose_name='Sub-Proteins', related_name='sub')
    pi = FloatField(verbose_name="Protein PI")
    mass = FloatField(verbose_name="Protein Mass (Da)")

    #html formatting
    def get_as_html_sequence(self, species, is_user_anonymous):
        from cyano.helpers import format_sequence_as_html
        return format_sequence_as_html(species, self.sequence)

    def get_as_html_structure(self, species, is_user_anonymous):
        from .helpers import overlaps, create_detail_fieldset

        length = self.length

        W = 636
        geneHeight = 20
        featureHeight = 10

        geneY = 2
        chrY = geneY + geneHeight + 4
        featureY = chrY + 1 + 2

        start_coordinate = 0
        end_coordinate = length

        feature_draw = []

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

        chro = '<text x="%s" y="%s" style="text-anchor:end;">%s</text><line x1="%s" x2="%s" y1="%s" y2="%s" /><text x="%s" y="%s" style="text-anchor:start;">%s</text>' % (
            chrL - 4, chrY, start_coordinate,
            chrL, chrR, chrY, chrY,
            chrR + 2, chrY, end_coordinate)

        def draw_segment(feature):
            opacity = 1
            coordinate = feature.coordinate
            item_length = feature.length

            x = chrL + float(coordinate - start_coordinate) / length * chrW
            w = chrL + float(coordinate + item_length - 1 - start_coordinate) / length * chrW

            x = max(chrL, min(chrR, x))
            w = max(chrL, min(chrR, w)) - x

            fieldsets = deepcopy(feature._meta.fieldsets)
            fielddata = create_detail_fieldset(species, feature, fieldsets, False)

            tip_title = fielddata[0][0]
            tip_text = StringIO.StringIO()

            for dataset in fielddata:
                for item in dataset[1]["fields"]:
                    tip_text.write("<br><b>" + item["verbose_name"] + "</b>: "+str(item["data"]))

            template = loader.get_template("cyano/genome/draw_feature.html")

            context_dict = {'h': featureHeight,
                            'title': tip_title,
                            'text': tip_text.getvalue(),
                            'url': feature.protein.get_absolute_url(species),
                            'x': x,
                            'w': w,
                            'opacity': opacity}

            new_item = [x, x + w]
            inserted = False
            for i, row in enumerate(feature_draw):
                if any(overlaps(item, new_item, 5) for item in row):
                    continue
                # Space left -> insert
                row.append(new_item)
                inserted = True
                context_dict.update({
                    'y': featureY + i * (featureHeight + 2)
                })
                break

            if not inserted:
                # Create new row
                context_dict.update({
                    'y': featureY + len(feature_draw) * (featureHeight + 2)
                })
                feature_draw.append([new_item])

            return template.render(Context(context_dict))


        #features
        featureStyle = '.features rect{fill:#%s;}' % (colors[2], )
        features = StringIO.StringIO()

        for feature in self.protein_details.all():
            features.write(draw_segment(feature))

        H = 2 + geneHeight + 2 + 4 + 1 * (2 + len(feature_draw) * (featureHeight + 2)) + 2

        return '<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="%s" height="%s" viewport="0 0 %s %s"><style>%s%s%s</style><g class="chr">%s</g><g class="features">%s</g></svg>' % (
            W, H, W, H, style, chrStyle, featureStyle, chro, features.getvalue())

    def get_as_html_ambiguous_proteins(self, species, is_user_anonymous):
        results = []

        for ambiguous_protein in self.ambiguous.all():
            try:
                res = Protein.objects.for_species(species).get(
                    parent=MassSpectrometryJob.objects.all()[0], wid__startswith=ambiguous_protein.value + "_")
                results.append('<a href="%s">%s</a>' % (res.get_absolute_url(species), res.wid))
            except ObjectDoesNotExist:
                results.append(ambiguous_protein.value)

        return format_list_html(results, comma_separated=True)

    def get_as_html_sub_proteins(self, species, is_user_anonymous):
        results = []

        for sub_protein in self.sub.all():
            try:
                res = Protein.objects.for_species(species).get(
                    parent=MassSpectrometryJob.objects.all()[0], wid__startswith=sub_protein.value + "_")
                results.append('<a href="%s">%s</a>' % (res.get_absolute_url(species), res.wid))
            except ObjectDoesNotExist:
                results.append(sub_protein.value)

        return format_list_html(results, comma_separated=True)

    class Meta:
        concrete_entry_model = True
        fieldsets = [
            ('Type', {'fields': ['model_type']}),
            ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}),
            ('Classification', {'fields': ['type']}),
            ('Related', {'fields': [
                'parent',
                {'verbose_name': 'Ambiguous Proteins', 'name': 'ambiguous_proteins'},
                {'verbose_name': 'Sun-Proteins', 'name': 'sub_proteins'},
            ]}),
            ('Structure', {'fields': [
                'prosthetic_groups', 'chaperones', 'dna_footprint',
                {'verbose_name': 'Sequence', 'name': 'sequence'},
                {'verbose_name': 'Structure', 'name': 'structure'},
            ]}),
            ('Regulation', {'fields': ['regulatory_rule']}),
            ('Function', {'fields': [
                {'verbose_name': 'Enzyme', 'name': 'enzyme_participants'},
                {'verbose_name': 'Transcriptional regulation', 'name': 'transcriptional_regulations'},
                {'verbose_name': 'Protein folding substrates', 'name': 'chaperone_substrates'},
                {'verbose_name': 'Reaction participant', 'name':'reaction_stoichiometry_participants'},
                {'verbose_name': 'Complex subunit', 'name':'protein_complex_biosythesis_participants'},
            ]}),
            ('Statistics', {'fields': [
                'score',
                'coverage',
                'pi',
                'mass'
            ]}),
            {'inline': 'protein_details'},
            ('Parameters', {'fields': ['parameters']}),
            ('Comments', {'fields': ['comments', 'publication_references']}),
            ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]}),
            ]
        field_list = [
            'id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'prosthetic_groups', 'chaperones', 'dna_footprint', 'regulatory_rule', 'score', 'coverage', 'pi', 'mass', 'sequence', 'comments', 'publication_references', 'created_detail', 'detail'
            ]
        facet_fields = ['type', 'chaperones', 'dna_footprint__binding', 'dna_footprint__region']
        verbose_name = 'Mass Spectrometry Protein'
        verbose_name_plural = 'Mass Spectrometry Proteins'
        wid_unique = False


class MassSpectrometryProteinDetail(EntryData):
    protein = ForeignKey(MassSpectrometryProtein, related_name="protein_details")
    sequence = TextField(verbose_name='Sequence')
    sequence_ptm = TextField(verbose_name='Sequence + PTMs')
    coordinate = PositiveIntegerField(verbose_name='Coordinate (nt)')
    length = PositiveIntegerField(verbose_name='Length (nt)')
    proteotypic = BooleanField(verbose_name='Proteotypic')
    zscore = FloatField(verbose_name='zscore')
    delta_mass = FloatField(verbose_name='Delta Mass (ppm)')
    mass = FloatField(verbose_name='Experimental Mass (m/z)')
    charge = IntegerField(verbose_name='Charge')
    retention_time = IntegerField(verbose_name='Retention Time (min)')
    theoretical_mass = FloatField(verbose_name='Theoretical Mass (Da)')
    missed_cleavages = IntegerField(verbose_name='Missed Cleavages')

    class Meta:
        verbose_name = "Mass Spectrometry Protein Detail"
        verbose_name_plural = "Mass Spectrometry Protein Details"
        fieldsets = [
            ('Mass Spectrometry Details', {'fields': [
                'sequence',
                'sequence_ptm',
                'proteotypic',
                'zscore',
                'delta_mass',
                'mass',
                'charge',
                'retention_time',
                'theoretical_mass',
                'missed_cleavages'
            ]}),
        ]
        field_list = [
            'id', 'chromosome_feature', 'chromosome', 'coordinate', 'length',
        ]



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
                tmp2.append(ref.get_citation(Species.objects.for_wid(ev.species), cross_references=True))
            tmp += '<div style="margin-top:4px;"><i>References</i><ol style="margin-top:0px"><li style="margin-bottom:0px">%s</li></ol></div>' % '</li><li style="margin-bottom:0px">'.join(tmp2)

        txt.append(tmp)

    return '<ul style="margin-top:4px;"><li style="margin-top:4px;">%s</li></ul>' % '</li><li style="margin-top:4px;">'.join(txt)

def sub_rate_law(species):
    def inner_method(match):
        if match.group(0) == 'Vmax':
            return '<i>V</i><sub>max</sub>'
        if match.group(0)[0:2] == 'Km':
            return '<i>K</i><sub>m%s</sub>' % match.group(0)[2:]
        try:
            obj = SpeciesComponent.objects.for_species(species).for_wid(match.group(0))
            return '[<a href="%s">%s</a>]' % (obj.get_absolute_url(species), obj.wid)
        except:
            return match.group(0)
    return inner_method
