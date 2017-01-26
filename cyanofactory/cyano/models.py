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
import json
import math
import re
from typing import List

from django.contrib.contenttypes.fields import GenericForeignKey
from django.contrib.contenttypes.models import ContentType
from django.db.models.aggregates import Count
from django.db.models.fields import NullBooleanField
import subprocess
import sys
from datetime import datetime
from io import StringIO, BytesIO

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqRecord, SeqIO

from django.contrib.auth.models import User, Permission
from django.contrib.auth.models import Group
from django.core import validators
from django.core.exceptions import ValidationError, ObjectDoesNotExist
from django.core.urlresolvers import reverse
from django.db.models import F, Model, OneToOneField, CharField, IntegerField, URLField, PositiveIntegerField,\
    FloatField, BooleanField, SlugField, TextField, DateTimeField, options, permalink, SET_NULL, signals
from django.db.models.fields.related import ForeignKey, ManyToManyField
from django.db.models.query import QuerySet, Prefetch
from django.template import loader, Context
from django.dispatch.dispatcher import receiver
from django.db.models.signals import m2m_changed
from model_utils.managers import InheritanceManager

from .templatetags.templatetags import set_time_zone
from .cache import Cache
from Bio.SeqFeature import SeqFeature, FeatureLocation
from copy import deepcopy
from django.conf import settings

from django.db import models

from guardian.models import UserObjectPermissionBase
from guardian.models import GroupObjectPermissionBase
from guardian.shortcuts import assign_perm, remove_perm

def enum(**enums):
    return type(str('Enum'), (), enums)

PermissionEnum = enum(
    READ_NORMAL="view_normal",
    READ_DELETE="view_delete",
    READ_PERMISSION="view_permission",
    READ_HISTORY="view_history",
    WRITE_NORMAL="change_entry",
    WRITE_DELETE="delete_entry",
    WRITE_PERMISSION="edit_permission",
    ACCESS_SPECIES="access_species",
    CREATE_SPECIES="create_mutant",
    ACCESS_SBGN="access_sbgn"
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

NUCLEOTIDE_SUBSTITUTION = {
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'T': 'T',
    'U': 'U',
    'R': 'GA',
    'Y': 'TC',
    'K': 'GT',
    'M': 'AC',
    'S': 'GC',
    'W': 'AT',
    'B': 'GTC',
    'D': 'GAT',
    'H': 'ACT',
    'V': 'GCA',
    'N': 'AGCT',
    '-': ''
}

''' END: CHOICES '''

# add model options
options.DEFAULT_NAMES = options.DEFAULT_NAMES + ('listing', 'concrete_entry_model', 'fieldsets', 'field_list',
                                                 'facet_fields', 'clean', 'validate_unique', 'wid_unique',
                                                 'group_field')


def get_anonymous_user_instance(User):
    return User.objects.get(username='guest')

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


class GlobalPermission(Model):
    class Meta:
        permissions = (
            ("access_species", "Can access any species"),
            ("create_mutant", "Can create new species or mutants"),
            ("access_sbgn", "Can access SBGN map"),
        )


class GroupProfile(Model):
    """
    Additional information for groups
    """
    group = OneToOneField(Group, related_name="profile")
    description = CharField(max_length=255, blank=True, default='', verbose_name="Group description")

    def assign_perm(self, permission, obj=None):
        if obj is None:
            if "." in permission:
                app_label, codename = permission.split(1)
                perm = Permission.objects.filter(content_type__app_label=app_label, codename=codename)
            else:
                codename = permission
                perm = Permission.objects.filter(codename=codename)
                if perm.count() > 1:
                    raise ValueError("Permission is ambiguous: {}".format(codename))

            perm = perm.first()
            assign_perm("{}.{}".format(perm.content_type.app_label, perm.codename), self.group)
        else:
            old_cls = obj.__class__
            obj.__class__ = Entry

            try:
                assign_perm(permission, self.group, obj)
            finally:
                obj.__class__ = old_cls

    def has_perm(self, permission, obj=None):
        """
        Improved perm check respecting all groups
        """
        if obj is None:
            return self.group.permissions.filter(codename=permission).exists()
        else:
            return self.has_perms(permission, Entry.objects.filter(pk=obj.pk)).count()

    def has_perms(self, permission, objs):
        """
        Mass perm check on objs.
        Returns objects group has perms on.
        Only works for entrys
        """
        from django.contrib.auth.models import Permission
        from django.contrib.contenttypes.models import ContentType

        if len(objs) == 0:
            return Permission.objects.none()

        ct = ContentType.objects.get_for_model(Entry)
        try:
            perm = Permission.objects.get(content_type=ct, codename=permission)
        except ObjectDoesNotExist:
            return Permission.objects.none()

        obj_ids = objs.values_list("pk", flat=True)
        all_perms = EntryGroupObjectPermission.objects.filter(permission=perm, group=self.group, content_object_id__in=obj_ids).values_list("content_object_id", flat=True)

        return objs.model.objects.filter(pk__in=all_perms)

    def __str__(self):
        return self.group.name


class UserProfile(Model):
    """
    Additional information for users
    """
    user = OneToOneField(User, related_name="profile")
    affiliation = CharField(max_length=255, blank=True, default='', verbose_name='Affiliation')
    website = URLField(max_length=255, blank=True, default='', verbose_name='Website')
    phone = CharField(max_length=255, blank=True, default='', verbose_name='Phone')
    address = CharField(max_length=255, blank=True, default='', verbose_name='Address')
    city = CharField(max_length=255, blank=True, default='', verbose_name='City')
    state = CharField(max_length=255, blank=True, default='', verbose_name='State')
    zip = CharField(max_length=255, blank=True, default='', verbose_name='Zip')
    country = CharField(max_length=255, blank=True, default='', verbose_name='Country')
    force_password_change = BooleanField(default=False)

    @staticmethod
    def get_profile(user: User):
        if not user.is_authenticated():
            # Special handling for guests
            return UserProfile.objects.get(user__username="guest")
        else:
            return user.profile

    def get_name(self):
        if len(self.user.first_name) > 0:
            name = self.user.first_name
            if len(self.user.last_name) > 0:
                name += " " + self.user.last_name
        elif len(self.user.last_name) > 0:
            name = self.user.last_name
        else:
            name = self.user.username
        return name

    def get_website(self):
        if not self.website.startswith("http://") and not self.website.startswith("https://"):
            return "http://" + self.website
        return self.website

    def get_global_perms(self):
        from django.contrib.auth.models import Permission

        if self.user.is_superuser:
            return Permission.objects.all()

        return self.user.user_permissions.all() |\
               Permission.objects.filter(pk__in=self.get_groups().values_list("group__permissions", flat=True))

    def get_groups(self):
        """
        Returns all groups where the user is in, including "Everybody" and
        (if the user isn't a guest) "Registred"
        
        :returns: GroupProfile List
        """
        groups = self.user.groups.values_list("pk", flat=True)
        profiles = GroupProfile.objects.filter(pk__in=groups)

        if self.user.username != "guest":
            profiles |= GroupProfile.objects.filter(group__name="Registred")

        profiles |= GroupProfile.objects.filter(group__name="Everybody")

        return profiles

    def assign_perm(self, permission, obj=None):
        if obj is None:
            if "." in permission:
                app_label, codename = permission.split(1)
                perm = Permission.objects.filter(content_type__app_label=app_label, codename=codename)
            else:
                codename = permission
                perm = Permission.objects.filter(codename=codename)
                if perm.count() > 1:
                    raise ValueError("Permission is ambiguous: {}".format(codename))

            perm = perm.first()
            assign_perm("{}.{}".format(perm.content_type.app_label, perm.codename), self.user)
        else:
            old_cls = obj.__class__
            obj.__class__ = Entry

            try:
                assign_perm(permission, self.user, obj)
            finally:
                obj.__class__ = old_cls

    def has_perm(self, permission: str, obj=None) -> bool:
        """
        Improved perm check respecting all groups
        """
        from django.contrib.auth.models import Group

        if obj is None:
            if self.user.is_superuser:
                return True
            return self.get_global_perms().filter(codename=permission).exists()
        else:
            return self.has_perms(permission, Entry.objects.filter(pk=obj.pk)).count() == 1

    def has_perms(self, permission: str, objs):
        """
        Mass perm check on objs.
        Returns objects user has perms on.
        Only works for entrys
        """
        from django.contrib.auth.models import Permission
        from django.contrib.contenttypes.models import ContentType

        if len(objs) == 0:
            return Permission.objects.none()

        if self.user.is_superuser:
            return objs

        ct = ContentType.objects.get_for_model(Entry)
        try:
            perm = Permission.objects.get(content_type=ct, codename=permission)
        except ObjectDoesNotExist:
            return Permission.objects.none()

        obj_ids = objs.values_list("pk", flat=True)
        all_user_perms = EntryUserObjectPermission.objects.filter(permission=perm, user=self.user, content_object_id__in=obj_ids).values_list("content_object_id", flat=True)
        obj_with_user_perms = objs.model.objects.filter(pk__in=all_user_perms)

        groups = self.get_groups().values_list("pk", flat=True)

        all_group_perms = EntryGroupObjectPermission.objects.filter(group__in=groups, permission=perm, content_object_id__in=obj_ids).values_list("content_object_id", flat=True)
        obj_with_group_perms = objs.model.objects.filter(pk__in=all_group_perms)

        return obj_with_user_perms | obj_with_group_perms

    def __str__(self):
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

    def __str__(self):
        return self.table_name + " - " + self.model_name

class RevisionDetail(Model):
    user = ForeignKey(UserProfile, verbose_name = "Modified by", related_name = '+', editable = False)
    date = DateTimeField(default=datetime.now, verbose_name = "Modification date")
    reason = TextField(blank=True, default='', verbose_name='Reason for edit')

class Revision(Model):
    """To allow reverting of edits all edit operations are stored in this table.
    
    Always only the new value of a single cell is stored to save memory.
    
    :Columns:
        * ``current``: Reference to the entry
        * ``detail``: Contains additional information about the edit
        * ``action``: Type of the operation (Update, Insert, Delete)
        * ``column``: Table and column where the modification occured
        * ``new_value``: New value in the cell of the column
    """
    current = GenericForeignKey()
    content_type = ForeignKey(ContentType)
    object_id = IntegerField(db_index=True, verbose_name="Current version primary key")
    detail = ForeignKey(RevisionDetail, verbose_name='Details about this revision', related_name='revisions', editable=False)
    action = CharField(max_length=1, choices=CHOICES_DBOPERATION)
    new_data = TextField(blank=True)

    def __str__(self):
        return str(self.current)

    class Meta:
        ordering = ['detail__date']

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

    def __str__(self):
        arr = []
        for field in self.__class__._meta.fields:
            if not field.auto_created and getattr(self, field.name) is not None and not (isinstance(getattr(self, field.name), (str, unicode, )) and getattr(self, field.name) == ''):
                arr.append('%s: %s' % (field.name, getattr(self, field.name)))
        for field in self.__class__._meta.many_to_many:
            arr.append('%s: %s' % (field.name, getattr(self, field.name).all()))
        return ', '.join(arr)

    class Meta:
        ordering = ['value', 'units']
        verbose_name = 'Evidence'
        verbose_name_plural = 'Evidence'

class EntryDataQuerySet(QuerySet):
    def get_or_create_with_revision(self, revision, *args, **kwargs):
        try:
            return self.get(*args, **kwargs)
        except ObjectDoesNotExist:
            obj = self.model(*args, **kwargs)
            obj.save(revision)
            return obj

class EntryData(Model):
    objects = EntryDataQuerySet.as_manager()

    def __str__(self):
        arr = []
        txt = ''
        nFields = 0
        for field in self.__class__._meta.fields:
            if not field.auto_created:
                nFields += 1
                try:
                    txt = getattr(self, field.name)
                    arr.append('%s: %s' % (field.name, txt))
                except ObjectDoesNotExist:
                    pass
        for field in self.__class__._meta.many_to_many:
            nFields += 1
            txt = getattr(self, field.name).all()
            arr.append('%s: %s' % (field.name, txt))

        if nFields == 1:
            return txt
        else:
            return ', '.join(arr)

    def save(self, revision_detail, force_revision=True, *args, **kwargs):
        if self.pk is not None:
            # Update item
            action = "U"

            # can this be optimized? Probably not
            old_item = self._meta.concrete_model.objects.get(pk=self.pk)
        else:
            # New entry (no primary key)
            action = "I"

            old_item = None

        super(EntryData, self).save(*args, **kwargs)

        save_data = {}

        fields = self._meta.fields
        # Remove primary keys and some fields that don't need revisioning:
        fields = filter(lambda x: not x.primary_key, fields)

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
                save_data[field.name] = new_value

        if force_revision or len(save_data) > 0:
            ##print str(self) + ": revisioning", len(save_data), "items"
            # Don't update the detail when actually nothing changed for that entry

            # Try fetching an existing revision
            try:
                r = Revision.objects.get(object_id=self.pk,
                                         content_type_id=ContentType.objects.get_for_model(self._meta.concrete_model).pk,
                                         detail_id=revision_detail.pk,
                                         action=action)
                new_data = json.loads(r.new_data)
                new_data.update(save_data)
                r.new_data = json.dumps(new_data)
            except ObjectDoesNotExist:
                r = Revision(current=self,
                             object_id=self.pk,
                             detail_id=revision_detail.pk,
                             action=action,
                             new_data=json.dumps(save_data))
            ##print save_data
            r.save()

    def delete(self, revision_detail, using=None):
        r = Revision(current=self,
                     object_id=self.pk,
                     detail_id=revision_detail.pk,
                     action="D",
                     new_data="{}")
        r.save()

        super(EntryData, self).delete(using=using)

    class Meta:
        abstract = True

class EvidencedEntryData(EntryData):
    evidence = ManyToManyField(Evidence, blank=True, verbose_name='Evidence')

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
        print("WARN: Reverse support not implemented:",sender,instance,pk_set)
        return

    if not pk_set is None and len(pk_set) > 0 and isinstance(instance, Entry):
        if action == "pre_add":
            # target ist species?
            model_name, field_name = sender._meta.object_name.split("_", 1)
            if model_name == "Species":
                species = model.objects.get(pk=list(pk_set)[0])
                if species.species_components.filter(model_type = instance.model_type, wid = instance.wid).exists():
                    raise ValidationError("Species {}: Wid {} already in use for model_type .{}".format(species.wid, instance, instance._meta.object_name))

        if action == "post_add":
            #print "m2m-add:",TableMetaManyToMany.get_by_m2m_model_name(sender._meta.object_name),instance,pk_set

            # Try extracting the field name from object_name...
            model_name, field_name = sender._meta.object_name.split("_", 1)

            # Test if that field is correct
            field = getattr(instance, field_name)

            # Fetch latest revision for instance
            revision = Revision.objects.filter(
                object_id=instance.pk,
                content_type=ContentType.objects.get_for_model(instance._meta.concrete_model)).last()

            new_data = json.loads(revision.new_data)
            item = new_data.get(field_name, [])
            item += list(pk_set)
            new_data[field_name] = item
            revision.new_data = json.dumps(new_data)
            revision.save()

        elif action == "post_remove":
            #print "m2m-del:",TableMetaManyToMany.get_by_m2m_model_name(sender._meta.object_name),instance,pk_set

            # Try extracting the field name from object_name...
            model_name, field_name = sender._meta.object_name.split("_", 1)

            # Test if that field is correct
            field = getattr(instance, field_name)

            # Fetch latest revision for instance
            revision = Revision.objects.filter(
                object_id=instance.pk,
                content_type=ContentType.objects.get_for_model(instance._meta.concrete_model)).last()

            new_data = json.loads(revision.new_data)
            item = new_data.get(field_name, [])

            for x in list(pk_set):
                try:
                    item.remove(x)
                except ValueError:
                    #print "fixme?"
                    pass
            new_data[field_name] = item
            revision.new_data = json.dumps(new_data)
            revision.save()


class EntryQuerySet(QuerySet):
    def for_wid(self, wid, get=True, create=False, creation_status=False):
        if get:
            try:
                if not creation_status:
                    return self.get(wid=wid)

                return self.get(wid=wid), False
            except ObjectDoesNotExist:
                if not create:
                    raise

                if not creation_status:
                    return self.model(wid=wid)

                return self.model(wid=wid), True

        if create:
            raise ValueError("Create only compatible with get")

        return self.filter(wid=wid)

    def for_permission(self, permission, user):
        return self.filter(pk__in=user.has_perms(permission, self).values_list("pk", flat=True))


class SpeciesComponentQuerySet(EntryQuerySet):
    def for_species(self, species):
        return self.filter(species=species.pk)

    def for_permission(self, permission, user):
        return self.filter(pk__in=user.has_perms(permission, self).values_list("pk", flat=True))


class AbstractEntry(Model):
    objects = EntryQuerySet.as_manager()
    child_objects = InheritanceManager()

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
    synonyms = ManyToManyField(Synonym, blank=True, related_name='entry', verbose_name='Synonyms')
    comments = TextField(blank=True, default='', verbose_name='Comments')

    def __str__(self):
        return self.wid

    def natural_key(self):
        return self.wid
    
    def get_model(self):
        return TableMeta.get_by_id(self.model_type_id)

    def get_name_or_wid(self):
        return self.name or self.wid

    def last_revision(self):
        return Revision.objects.filter(
            content_type_id=ContentType.objects.get_for_model(self).pk,
            object_id=self.pk).prefetch_related('detail').last()

    def first_revision(self):
        return Revision.objects.filter(
            content_type_id=ContentType.objects.get_for_model(self).pk,
            object_id=self.pk).prefetch_related('detail').first()

    def get_permissions(self):
        return [
            EntryUserObjectPermission.objects.filter(content_object_id=self.pk),
            EntryGroupObjectPermission.objects.filter(content_object_id=self.pk)
        ]

    def save(self, revision_detail, force_revision=True, *args, **kwargs):
        from cyano.helpers import slugify

        if not self.wid or self.wid != slugify(self.wid):
            # Slugify the WID
            raise ValidationError("Wid must be slug!")

        if self.pk is not None:
            # Update item
            action = "U"

            # can this be optimized? Probably not
            old_item = self._meta.concrete_model.objects.get(pk = self.pk)
            super(Entry, self).save(*args, **kwargs)
        else:
            # New entry (no primary key)
            # The latest entry is not revisioned to save space (and time)
            action = "I"
            cache_model_type = TableMeta.get_by_model_name(self._meta.object_name)
            self.model_type = cache_model_type

            old_item = None

            super(Entry, self).save(*args, **kwargs)

        save_data = {}

        fields = self._meta.fields
        # Remove primary keys and some fields that don't need revisioning:
        fields = filter(lambda x: (not x.primary_key) and
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
                save_data[field.name] = new_value

        if force_revision or len(save_data) > 0:
            ##print self.wid + ": revisioning", len(save_list), "items"
            # Don't update the detail when actually nothing changed for that entry
            #self.detail = revision_detail

            # Try fetching an existing revision
            try:
                r = Revision.objects.get(object_id=self.pk,
                                         content_type_id=ContentType.objects.get_for_model(self._meta.concrete_model).pk,
                                         detail_id=revision_detail.pk,
                                         action=action)
                new_data = json.loads(r.new_data)
                new_data.update(save_data)
                r.new_data = json.dumps(new_data)
            except ObjectDoesNotExist:
                r = Revision(current=self,
                             object_id=self.pk,
                             detail_id=revision_detail.pk,
                             action=action,
                             new_data=json.dumps(save_data))
            ##print save_data
            r.save()

    def delete(self, revision_detail, using=None):
        r = Revision(current=self,
                     object_id=self.pk,
                     detail_id=revision_detail.pk,
                     action="D",
                     new_data="{}")
        r.save()

        super(Entry, self).delete(using=using)

    #html formatting
    def get_as_html_synonyms(self, is_user_anonymous):
        return format_list_html([x.name for x in self.synonyms.all()], comma_separated=True)

    def get_as_html_cross_references(self, is_user_anonymous):
        results = []
        for cr in self.cross_references.all():
            results.append('%s: <a href="%s">%s</a>' % (cr.source, reverse("db_xref.views.dbxref", kwargs={"source" : cr.source, "xid" : cr.xid}), cr.xid ))
        return format_list_html(results, separator=', ')

    def get_as_html_created_user(self, is_user_anonymous):
        revision = self.first_revision()
        detail = revision.detail
        if is_user_anonymous:
            return '%s<br>%s' % (detail.date.strftime("%Y-%m-%d %H:%M:%S"), detail.reason)
        else:
            user = detail.user.user
            return '<a href="%s">%s %s</a> on %s<br>%s' % (user.get_absolute_url(), user.first_name, user.last_name, detail.date.strftime("%Y-%m-%d %H:%M:%S"), detail.reason)

    def get_as_html_last_updated_user(self, is_user_anonymous):
        revision = self.last_revision()
        detail = revision.detail
        if is_user_anonymous:
            return '%s<br>%s' % (detail.date.strftime("%Y-%m-%d %H:%M:%S"), detail.reason)
        else:
            user = detail.user.user
            return '<a href="%s">%s %s</a> on %s<br>%s' % (user.get_absolute_url(), user.first_name, user.last_name, detail.date.strftime("%Y-%m-%d %H:%M:%S"), detail.reason)

    def get_as_fasta(self):
        raise NotImplementedError("FASTA export not supported for %s" % self._meta.verbose_name_plural)

    def get_fasta_header(self):
        return ">" + self.wid + "|" + self.name

    def get_as_genbank(self):
        raise NotImplementedError("GenBank export not supported for %s" % self._meta.verbose_name_plural)
    
    def get_as_sbml(self):
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
            'id', 'wid', 'name', 'synonyms', 'comments'
            ]
        listing = ['wid', 'name']
        facet_fields = []
        ordering = ['wid']
        get_latest_by = 'createdDate'
        verbose_name = 'Entry'
        verbose_name_plural = 'Entries'
        wid_unique = False

        permissions = (
            ('view_normal', 'View entry'),
            ('view_delete', 'View deleted revisions of entry'),
            ('view_permission', 'View permissions of entry'),
            ('view_history', 'View older (not deleted) revisions of entry'),
            ('edit_permission', 'Allow modifying of permissions'),
        )


class EntryUserObjectPermission(UserObjectPermissionBase):
    content_object = models.ForeignKey(Entry)


class EntryGroupObjectPermission(GroupObjectPermissionBase):
    content_object = models.ForeignKey(Entry)


class AbstractSpeciesComponent(Entry):
    objects = SpeciesComponentQuerySet.as_manager()

    class Meta:
        abstract = True


class SpeciesComponent(AbstractSpeciesComponent):
    '''
    Contains all components that belong to an organism.
    Improves the lookup speed and allows inheritance.
    
    * ``cross_references``: Databases referencing that entry
    * ``publication_references``: Publications referencing that entry
    '''
    species = ForeignKey('Species', related_name='components', verbose_name='Species')
    type = ManyToManyField('Type', blank=True, related_name='members', verbose_name='Type')
    cross_references = ManyToManyField("CrossReference", blank=True, related_name='cross_referenced_components', verbose_name='Cross references')
    publication_references = ManyToManyField("PublicationReference", blank=True, related_name='publication_referenced_components', verbose_name='Publications')
    parent = ForeignKey('self', blank=True, null=True, on_delete=SET_NULL, related_name='children', verbose_name='Parent')

    #getters

    def delete(self, revision_detail, using=None):
        super(SpeciesComponent, self).delete(revision_detail, using=using)

    #@permalink
    def get_absolute_url(self, history_id=None, api=False):
        species_key = "species/%s" % self.species_id
        species = Cache.try_get(species_key, lambda: self.species, 60)

        return "%s/%s%s/%s/%s%s" % (
            settings.ROOT_URL,
            "api/" if api else "",
            species.wid,
            TableMeta.get_by_id(self.model_type_id).model_name,
            "%s/" % history_id if history_id is not None else "",
            self.wid
        )

    def get_absolute_url_api(self):
        return self.get_absolute_url(api=True)

    @classmethod
    def get_model_url(cls, species):
        return "%s/%s/%s" % (settings.ROOT_URL, species.wid, cls._meta.object_name)

    def get_all_references(self):
        return self.publication_references.all() | PublicationReference.objects.filter(evidence__species_component__id = self.id)

    @classmethod
    def get_statistics(cls, species):
        # API:
        # [[0, "<a href="Chromosome">, 123],

        amount = cls.objects.for_species(species).count()
        url = cls.get_model_url(species)
        types = Type.objects.filter(members__in=cls.objects.for_species(species)).order_by("wid").values("wid").annotate(count=Count("wid"))

        if amount == 0:
            return None

        idx = 0

        lst = [[idx, cls._meta.verbose_name_plural, amount, None, url]]

        for typ in types:
            lst.append([idx+1, typ["wid"], typ["count"], None, "{}?type={}".format(url, typ["wid"])])

        return lst

    #html formatting
    def get_as_html_parameters(self, is_user_anonymous):
        results = []
        for p in self.parameters.all():
            results.append('<a href="%s">%s</a>: <i>%s</i> = %s %s' % (p.get_absolute_url(self.species), p.wid, p.name, p.value.value, p.value.units))
        return format_list_html(results)

    def get_as_html_comments(self, is_user_anonymous):
        txt = self.comments

        #provide links to references
        return re.sub(r'\[(PUB_\d{4,4})(, PUB_\d{4,4})*\]',
            lambda match: '[' + ', '.join(['<a href="%s">%s</a>' % (reverse('cyano.views.detail', kwargs={'species_wid':self.species.wid, 'model_type': 'PublicationReference', 'wid': x}), x, ) for x in match.group(0)[1:-1].split(', ')]) + ']',
            txt)

    def get_as_html_publication_references(self, is_user_anonymous):
        results = {}
        for r in self.get_all_references():
            key = r.authors + ' ' + r.editors
            results[key] = r.get_citation(cross_references=True)

        keys = list(results.keys())
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
            'id', 'wid', 'name', 'synonyms', 'cross_references', 'type',  'comments', 'publication_references'
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
        return self.get_empirical_formula().get_molecular_weight()

    def get_atoms(self):
        return sum(self.get_empirical_formula().values())

    def get_absorbance_factor(self):
        return 1 / self.get_extinction_coefficient()

    #html formatting        
    def get_as_html_empirical_formula(self, is_user_anonymous):
        return ""
        return self.get_empirical_formula().get_as_html()

    def get_as_html_molecular_weight(self, is_user_anonymous):
        return ""
        return self.get_molecular_weight()

    def get_as_html_extinction_coefficient(self, is_user_anonymous):
        return self.get_extinction_coefficient()

    def get_as_html_absorbance_factor(self, is_user_anonymous):
        return self.get_absorbance_factor()

    def get_as_html_pi(self,  is_user_anonymous):
        if hasattr(self, 'pi'):
            return self.pi
        else:
            return self.get_pi()

    def get_as_html_modification_reactions(self, is_user_anonymous):
        return ""# FIXME
        results = []
        for obj in self.modification_reactions.all():
            if len(obj.reactions.all()) == 0:
                continue
            rxn = obj.reactions.all()[0]
            results.append('<a href="%s">%s</a><br/>%s' % (rxn.get_absolute_url(), rxn.name, rxn.get_as_html_stoichiometry(is_user_anonymous)))
        return format_list_html(results, vertical_spacing=True)

    def get_as_html_reaction_stoichiometry_participants(self, is_user_anonymous):
        results = []
        for obj in self.reaction_stoichiometry_participants.all():
            if len(obj.reactions.all()) == 0:
                continue
            rxn = obj.reactions.all()[0]
            results.append('<a href="%s">%s</a><br/>%s' % (rxn.get_absolute_url(), rxn.name, rxn.get_as_html_stoichiometry(is_user_anonymous)))
        return format_list_html(list(set(results)), vertical_spacing=True)

    def get_as_html_protein_complex_biosythesis_participants(self, is_user_anonymous):
        results = []
        for obj in self.protein_complex_biosythesis_participants.all():
            if len(obj.protein_complexes.all()) == 0:
                continue
            pc = obj.protein_complexes.all()[0]
            results.append('<a href="%s">%s</a><br/>%s' % (pc.get_absolute_url(), pc.name, pc.get_as_html_biosynthesis(is_user_anonymous)))
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
            'id', 'wid', 'name', 'synonyms', 'cross_references', 'type',  'comments', 'publication_references'
            ]
        facet_fields = ['type']
        verbose_name = 'Molecule'
        verbose_name_plural = 'Molecules'
        wid_unique = False

class Protein(Molecule):
    #parent pointer
    parent_ptr_molecule = OneToOneField(Molecule, related_name='child_ptr_protein', parent_link=True, verbose_name='Molecule')

    #additional fields
    prosthetic_groups = ManyToManyField(ProstheticGroupParticipant, blank=True, related_name='proteins', verbose_name='Prosthetic groups')
    chaperones = ManyToManyField('self', symmetrical=False, blank=True, related_name='chaperone_substrates', verbose_name='Chaperones')
    dna_footprint = ForeignKey(DNAFootprint, null=True, blank=True, related_name='proteins', verbose_name='DNA footprint')
    regulatory_rule = ForeignKey(EntryCharData, null=True, blank=True, on_delete=SET_NULL, verbose_name='Regulatory rule', related_name='+')

    #html formatting
    def get_as_html_prosthetic_groups(self,  is_user_anonymous):
        results = []
        for p in self.prosthetic_groups.all():
            results.append(format_with_evidence(list_item = True, obj = p, txt = '(%s) <a href="%s">%s</a> [<a href="%s">%s</a>]' % (p.coefficient, p.metabolite.get_absolute_url(), p.metabolite.wid, p.compartment.get_absolute_url(), p.compartment.wid)))
        return format_list_html(results, force_list=True)

    def get_as_html_chaperones(self, is_user_anonymous):
        from cyano.helpers import format_list_html_url
        return format_list_html_url(self.chaperones.all(), self.species)

    def get_as_html_chaperone_substrates(self, is_user_anonymous):
        from cyano.helpers import format_list_html_url
        return format_list_html_url(self.chaperone_substrates.all(), self.species)

    def get_as_html_dna_footprint(self, is_user_anonymous):
        if self.dna_footprint is None:
            return None
        return format_with_evidence(obj = self.dna_footprint, txt = 'Length: %s (nt), Binding: %s, Region: %s' % (self.dna_footprint.length, self.dna_footprint.binding, self.dna_footprint.region))

    def get_as_html_regulatory_rule(self, is_user_anonymous):
        if self.regulatory_rule is not None and self.regulatory_rule.value is not None and self.regulatory_rule.value != '':
            return parse_regulatory_rule(self.regulatory_rule.value, {}, self.species.wid)
        return ''

    def get_as_html_transcriptional_regulations(self, is_user_anonymous):
        results = []
        for r in self.transcriptional_regulations.all():
            results.append('<a href="%s">%s</a>: <a href="%s">%s</a>' % (r.get_absolute_url(), r.wid, r.transcription_unit.get_absolute_url(), r.transcription_unit.wid))
        return format_list_html(results)

    def get_as_html_enzyme_participants(self, is_user_anonymous):
        results = []
        for obj in self.enzyme_participants.all():
            if len(obj.reactions.all()) == 0:
                continue
            rxn = obj.reactions.all()[0]
            results.append('<a href="%s">%s</a><br/>%s' % (rxn.get_absolute_url(), rxn.name, rxn.get_as_html_stoichiometry(is_user_anonymous)))
        return format_list_html(results, vertical_spacing=True)

    def get_as_html_half_life(self, is_user_anonymous):
        return self.get_half_life()

    def clean(self, all_obj_data=None, all_obj_data_by_model=None):
        #regulatory rule
        if self.regulatory_rule is not None and self.regulatory_rule.value is not None and self.regulatory_rule.value != '':
            parse_regulatory_rule(self.regulatory_rule.value, all_obj_data, self.species)

        #DNA footprint
        if self.dna_footprint is not None:
            chr_lens = []

            if all_obj_data_by_model is None:
                species_id = Species.objects.values('id').get(wid=self.species)['id']
                for obj in Genome.objects.values('length').filter(species__id=species_id):
                    chr_lens.append(obj['length'])
            else:
                for obj in all_obj_data_by_model['Genome']:
                    if isinstance(obj, Entry):
                        chr_lens.append(obj.length)
                    else:
                        chr_lens.append(obj['length'])

            if self.dna_footprint.length > max(chr_lens):
                raise ValidationError({'dna_footprint': 'Length must be less than chromosome length'})

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
            'id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'prosthetic_groups', 'chaperones', 'dna_footprint', 'regulatory_rule', 'comments', 'publication_references'
            ]
        facet_fields = ['type', 'chaperones', 'dna_footprint__binding', 'dna_footprint__region']
        verbose_name='Protein'
        verbose_name_plural = 'Proteins'
        wid_unique = False


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

    @classmethod
    def get_statistics(cls, species):
        amount = cls.objects.for_species(species).count()
        url = cls.get_model_url(species)

        if amount == 0:
            return None

        lst = [[0, cls._meta.verbose_name_plural, amount, None, url]]

        amount = Chromosome.objects.for_species(species).count()
        url = Chromosome.get_model_url(species)
        chromosomes = Chromosome.objects.for_species(species)

        if amount > 0:
            lst.append([1, Chromosome._meta.verbose_name_plural, amount, None, url])
            chr_sum = sum(c.length for c in chromosomes)
            lst.append([2, "Length", chr_sum, "nt"])
            lst.append([2, "GC-content", round((sum(c.get_gc_content() for c in chromosomes) / amount) * 100, 1), "%"])

        amount = Plasmid.objects.for_species(species).count()
        url = Plasmid.get_model_url(species)
        chromosomes = Plasmid.objects.for_species(species)

        if amount > 0:
            lst.append([1, Plasmid._meta.verbose_name_plural, amount, None, url])
            chr_sum = sum(c.length for c in chromosomes)
            lst.append([2, "Length", chr_sum, "nt"])
            lst.append([2, "GC-content", round((sum(c.get_gc_content() for c in chromosomes) / amount) * 100, 1), "%"])

        return lst

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
    def get_as_html_sequence(self, is_user_anonymous):
        from cyano.helpers import format_sequence_as_html
        return format_sequence_as_html(self.species, self.sequence, show_protein_seq=True)

    def get_as_html_diff_sequence(self, new_obj, is_user_anonymouse):
            from Bio import pairwise2
            pairwise2.MAX_ALIGNMENTS = 1
            align = pairwise2.align.globalxx(self.sequence[:100], new_obj[:100])
            output = StringIO()

            for o, n in zip(align[0][0], align[0][1]):
                output.write("|" if o == n else " ")

            align_graphic = output.getvalue()
            output.close()
            output = StringIO()

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

    def get_as_html_structure(self, is_user_anonymous, start_coordinate = None, end_coordinate = None, highlight_wid = None, zoom = 0):
        if zoom == 0:
            return self.get_as_html_structure_global(is_user_anonymous)
        else:
            return self.get_as_html_structure_local(is_user_anonymous, start_coordinate = start_coordinate, end_coordinate = end_coordinate, highlight_wid = highlight_wid)

    def get_as_html_structure_global(self, is_user_anonymous):
        from .helpers import shift, overlaps
        from collections import namedtuple

        attrib_class = type(str("StructureAttribClass"), (object,), dict())

        chr_attrib = attrib_class()

        chr_attrib.nt_per_segment = 1e4
        chr_attrib.height = 27
        gene_height = 10
        feature_height = 5
        num_colors = 20
        num_segments = int(math.ceil(self.length / chr_attrib.nt_per_segment))
        chr_attrib.left = len(str(num_segments * chr_attrib.nt_per_segment)) * 5
        chr_attrib.width = 636 - 4 - chr_attrib.left
        chr_attrib.one_nt_width = chr_attrib.width / chr_attrib.nt_per_segment
        chr_attrib.last_width = (self.length - chr_attrib.nt_per_segment * (num_segments - 1)) * chr_attrib.one_nt_width
        chr_attrib.top = -12

        fake_gene = Gene(model_type=TableMeta.get_by_model_name("Gene"))
        fake_cf = ChromosomeFeature(model_type=TableMeta.get_by_model_name("ChromosomeFeature"))

        # data
        chromosomes = []
        genes = [[] for x in range(num_segments)]
        features = [[] for x in range(num_segments)]
        promoters = [[] for x in range(num_segments)]
        tf_sites = [[] for x in range(num_segments)]

        GeneTuple = namedtuple("GeneTuple", "name wid coordinate direction length species transcription_units transcription_units__wid transcription_units__name")

        # TUs
        num_transcription_units = 0
        transcription_units_index = {}
        transcription_units = []

        genes_list = map(GeneTuple._make, self.genes.values_list(
            "name", "wid", "coordinate", "direction", "length", "species",
            "transcription_units", "transcription_units__wid", "transcription_units__name").all())

        all_feature_type_pks = [None] + list(Type.objects.filter(pk__in=self.features.values_list("chromosome_feature__type", flat=True)).distinct().order_by("parent").values_list("pk", flat=True))
        all_features_types_count = dict(map(lambda x: [x, 0], all_feature_type_pks))

        types = []

        for i in range(num_segments):
            attrib = attrib_class()

            attrib.left = chr_attrib.left
            attrib.right = chr_attrib.left + ((min(self.length, (i+1) * chr_attrib.nt_per_segment) - 1) % chr_attrib.nt_per_segment) / chr_attrib.nt_per_segment * chr_attrib.width
            attrib.y = chr_attrib.top + gene_height
            attrib.start = int(i * chr_attrib.nt_per_segment + 1)
            #attrib.end = int(i * chr_attrib.nt_per_segment + chr_attrib.nt_per_segment)
            chromosomes.append(attrib)

        def draw_gene(gene, tu):
            segment_index = int(math.floor((gene.coordinate - 1) / chr_attrib.nt_per_segment))

            #tip_text = 'Transcription unit: %s' % (tu or "(None)")

            fake_gene.wid = gene.wid
            fake_gene.species_id = pk=gene.species
            url = fake_gene.get_absolute_url()

            title = (gene.name or gene.wid).replace("'", "\'")

            ret = []

            segments = shift(list(range(num_segments)), int(segment_index))

            w_drawn = 0

            for i in itertools.cycle(segments):
                gene_attrib = attrib_class()
                gene_attrib.url = url
                gene_attrib.title = title
                gene_attrib.wid = gene.wid

                gene_attrib.direction = gene.direction

                gene_attrib.color = transcription_unit_index % num_colors

                if i == segment_index:
                    gene_attrib.x1 = chr_attrib.left + ((gene.coordinate - 1) % chr_attrib.nt_per_segment) / chr_attrib.nt_per_segment * chr_attrib.width
                else:
                    gene_attrib.x1 = chr_attrib.left

                w_needed = gene.length / chr_attrib.nt_per_segment * chr_attrib.width - w_drawn

                row = i + 1
                if row == num_segments:
                    w_space = chr_attrib.left + chr_attrib.last_width - gene_attrib.x1
                else:
                    w_space = chr_attrib.left + chr_attrib.width - gene_attrib.x1

                gene_attrib.y2 = chr_attrib.top + gene_height - 2  # + row_offset[i]
                gene_attrib.y1 = gene_attrib.y2 - gene_height

                if w_space < w_needed:
                    # Not enough space left on line

                    w = max(1, w_space)
                    w_drawn += w
                    gene_attrib.x2 = gene_attrib.x1 + w

                    if math.fabs(gene_attrib.x2 - gene_attrib.x1) > len(gene.wid) * 5:
                        gene_attrib.label = gene.wid
                    else:
                        gene_attrib.label = ''

                    if gene.direction == "r" and segment_index == i:
                        if gene.direction == "r" and gene_attrib.x1 + 5 > gene_attrib.x2:
                            gene_attrib.arrow_size = gene_attrib.x2 - gene_attrib.x1
                        else:
                            gene_attrib.arrow_size = 5
                        gene_attrib.arrow = True

                    genes[i].append(gene_attrib)
                else:
                    w = w_needed
                    gene_attrib.x2 = gene_attrib.x1 + w
                    if math.fabs(gene_attrib.x2 - gene_attrib.x1) > len(gene.wid) * 5:
                        gene_attrib.label = gene.wid
                    else:
                        gene_attrib.label = ''

                    if gene.direction == "f" or segment_index == i:
                        if gene_attrib.x2 - 5 < gene_attrib.x1:
                            gene_attrib.arrow_size = abs(gene_attrib.x1 - gene_attrib.x2)
                        else:
                            gene_attrib.arrow_size = 5
                        gene_attrib.arrow_size = -gene_attrib.arrow_size if gene_attrib.direction == "f" else gene_attrib.arrow_size
                        gene_attrib.arrow = True

                    genes[i].append(gene_attrib)
                    break

            return ret

        feature_rows = [[] for x in range(num_segments)]

        #promoters
        def add_segment(wid, coordinate, length, typ, url):
            segment_index = math.floor((coordinate - 1) / chr_attrib.nt_per_segment)

            segments = shift(list(range(num_segments)), int(segment_index))

            w_drawn = 0
            done = False
            for i in itertools.cycle(segments):
                feature_attrib = attrib_class()
                feature_attrib.wid = wid
                feature_attrib.coordinate = coordinate
                feature_attrib.length = length
                feature_attrib.url = url
                feature_attrib.title = wid.replace("'", "\'")

                feature_attrib.x = chr_attrib.left
                if i == segment_index:
                    feature_attrib.x += ((feature_attrib.coordinate - 1) % chr_attrib.nt_per_segment) * chr_attrib.one_nt_width
                else:
                    feature_attrib.x += 0

                w_needed = feature_attrib.length - w_drawn

                w_space = chr_attrib.nt_per_segment - feature_attrib.x

                if w_space < w_needed:
                    # Not enough space left on line
                    w = max(1, w_space)
                    w_drawn += w
                else:
                    w = w_needed
                    done = True

                #feature_attrib.x = chr_attrib.left + feature_attrib.x * chr_attrib.one_nt_width
                feature_attrib.y = 0
                feature_attrib.height = feature_height
                #feature_attrib.width = w

                feature_attrib.width = max(3, w * chr_attrib.one_nt_width)
                #y = row_offset[i] + j * (featureHeight + 1)

                if isinstance(typ, str):
                    pass
                else:
                    if typ and all_features_types_count[typ.pk] == 0:
                        types.append(typ)
                        all_features_types_count[typ.pk] = 1
                    else:
                        all_features_types_count[None] = 1
                    feature_attrib.color = all_feature_type_pks.index(typ.pk if typ else None) % num_colors
                    if typ:
                        typ.color = feature_attrib.color
                    feature_attrib.type = typ

                new_item = [feature_attrib.x, feature_attrib.x + feature_attrib.width, feature_attrib]
                inserted = False

                for j, row in enumerate(feature_rows[i]):
                    if any(overlaps(item, new_item, 5) for item in row):
                        continue
                    # Space left -> insert
                    row.append(new_item)
                    features[i].append(feature_attrib)
                    feature_attrib.y += j * (feature_height + 1)
                    inserted = True
                    break

                if not inserted:
                    # Create new row
                    feature_attrib.y += len(feature_rows[i]) * (feature_height + 1)
                    feature_rows[i].append([new_item])
                    features[i].append(feature_attrib)

                if done:
                    break

        for transcription_unit in transcription_units:
            if transcription_unit.promoter_35_coordinate is not None:
                tu_coordinate = transcription_unit.get_coordinate() + transcription_unit.promoter_35_coordinate
                tu_length = transcription_unit.promoter_35_length

                if not tu_coordinate > chr_attrib.end or tu_coordinate + tu_length - 1 < chr_attrib.start:
                    tip_title = transcription_unit.name or transcription_unit.wid
                    url = transcription_unit.get_absolute_url()

                    draw_segment(transcription_unit.wid, tu_coordinate, tu_length, tip_title, 'Promoter -35 box', url)

            if transcription_unit.promoter_10_coordinate is not None:
                tu_coordinate = transcription_unit.get_coordinate() + transcription_unit.promoter_10_coordinate
                tu_length = transcription_unit.promoter_10_length

                if not (tu_coordinate > chr_attrib.end or tu_coordinate + tu_length - 1 < chr_attrib.start):
                    tip_title = transcription_unit.name or transcription_unit.wid
                    url = transcription_unit.get_absolute_url()

                    promoters.append(draw_segment(transcription_unit.wid, tu_coordinate, tu_length, tip_title, 'Promoter -10 box', url))

            for tr in transcription_unit.transcriptional_regulations.all():
                if tr.binding_site is not None and not (tr.binding_site.coordinate > chr_attrib.end or tr.binding_site.coordinate + tr.binding_site.length - 1 < chr_attrib.start):
                    tip_title = transcription_unit.name or transcription_unit.wid
                    url = transcription_unit.get_absolute_url()

                    tf_sites.append(draw_segment(transcription_unit.wid, tr.binding_site.coordinate, tr.binding_site.length, tip_title, 'Transcription factor binding site', url))

        #features

        feature_values = self.features.all().order_by("coordinate").\
            values("coordinate", "length", "chromosome_feature__species",
                   "chromosome_feature__name", "chromosome_feature__wid",
                   "chromosome_feature__type", "chromosome_feature__type__name", "chromosome_feature__type__wid")

        for feature in feature_values:
            coordinate = feature["coordinate"]
            length = feature["length"]

            if feature["chromosome_feature__type"]:
                typ = attrib_class()
                typ.wid = feature["chromosome_feature__type__wid"]
                typ.name = feature["chromosome_feature__type__name"]
                typ.pk = feature["chromosome_feature__type"]
            else:
                typ = None
            fake_cf.wid = feature["chromosome_feature__wid"]
            fake_cf.species_id = feature["chromosome_feature__species"]
            url = fake_cf.get_absolute_url()

            add_segment(fake_cf.wid, coordinate, length, typ, url)

        for i, gene in enumerate(genes_list):
            if gene.transcription_units:
                tu_wid = gene.transcription_units__wid
                tu_name = gene.transcription_units__name

                tu = tu_name or tu_wid

                transcription_unit_index = transcription_units_index.get(tu_wid)
                if transcription_unit_index is None:
                    transcription_units.append(gene.transcription_units)
                    transcription_unit_index = num_transcription_units
                    transcription_units_index[tu_wid] = transcription_unit_index
                    num_transcription_units += 1
            else:
                transcription_unit_index = num_transcription_units
                num_transcription_units += 1
                tu = None

            draw_gene(gene, tu)

        row_offset = [chr_attrib.top + chr_attrib.height + 2]
        # Flatten lists and add offset
        for i, feature in enumerate(feature_rows):
            if i > 0:
                row_offset.append(chr_attrib.top + chr_attrib.height + 2 + (len(feature_rows[i - 1])*(feature_height + 1) + row_offset[-1]))

        for item, offset in zip(chromosomes, row_offset):
            item.y += offset

        for row, offset in zip(genes, row_offset):
            for item in row:
                item.y1 += offset
                item.y2 += offset

        for row, offset in zip(itertools.chain(features, promoters, tf_sites), itertools.cycle(row_offset)):
            for item in row:
                item.y += offset

        H = row_offset[-1] + len(feature_rows[-1]) * (feature_height + 2)

        c = Context({
            'species': self.species,
            'genes': [item for sublist in genes for item in sublist],
            'height': H,
            'chromosomes': chromosomes,
            'features': [item for sublist in features for item in sublist],
            'promoters': [item for sublist in promoters for item in sublist],
            'tf_sites': [item for sublist in tf_sites for item in sublist],
            'types': types
        })

        template = loader.get_template("cyano/fields/structure.html")
        rendered = template.render(c)

        return rendered

    def get_as_html_structure_local(self, is_user_anonymous, start_coordinate=None, end_coordinate=None, highlight_wid=[]):
        from .helpers import overlaps

        attrib_class = type(str("StructureAttribClass"), (object,), dict())

        gene_y = 2
        gene_height = 20
        feature_height = 10
        num_colors = 20

        # Chromosome
        chr_attrib = attrib_class()

        chr_attrib.start = start_coordinate
        chr_attrib.end = end_coordinate
        chr_attrib.y = gene_y + gene_height + 4
        chr_attrib.left = 4.5 * len(str(chr_attrib.start)) + 4
        chr_attrib.right = 636 - 4.5 * len(str(chr_attrib.end)) - 2 - 6
        chr_attrib.width = chr_attrib.right - chr_attrib.left
        chr_attrib.length = chr_attrib.end - chr_attrib.start + 1

        promoter_y = chr_attrib.y + 1 + 2
        feature_y = promoter_y

        # data
        genes = []
        features = []
        promoters = []
        tf_sites = []
        feature_draw = []

        # TUs
        num_transcription_units = 0
        transcription_units_index = {}
        transcription_units = []

        genes_list = self.genes.filter(
            coordinate__lte=chr_attrib.end,
            coordinate__gte=chr_attrib.start + 1 - F("length")
        ).prefetch_related('transcription_units', 'transcription_units__transcriptional_regulations').all()

        all_feature_type_pks = [None] + list(
            Type.objects.filter(pk__in=self.features.values_list("chromosome_feature__type", flat=True)).distinct()
            .order_by("parent").values_list("pk", flat=True)
        )

        all_features_types_count = dict(map(lambda x: [x, 0], all_feature_type_pks))

        types = []

        for i, gene in enumerate(genes_list):
            if len(gene.transcription_units.all()[:1]) == 1:
                transcription_unit = gene.transcription_units.all()[:1][0]

                transcription_unit_index = transcription_units_index.get(transcription_unit.wid)
                if transcription_unit_index is None:
                    transcription_units.append(transcription_unit)
                    transcription_unit_index = num_transcription_units
                    transcription_units_index[transcription_unit.wid] = transcription_unit_index
                    num_transcription_units += 1
            else:
                transcription_unit_index = num_transcription_units
                num_transcription_units += 1

            gene_attrib = attrib_class()
            gene_attrib.wid = gene.wid
            gene_attrib.x1 = chr_attrib.left + float(gene.coordinate - chr_attrib.start) / chr_attrib.length * chr_attrib.width
            gene_attrib.x2 = chr_attrib.left + float(gene.coordinate + gene.length - 1 - chr_attrib.start) / chr_attrib.length * chr_attrib.width

            gene_attrib.arrow = True

            gene_attrib.direction = gene.direction

            gene_attrib.arrow_size = 0
            if gene.direction == "r" and gene_attrib.x1 + 5 > gene_attrib.x2:
                gene_attrib.arrow_size = gene_attrib.x2 - gene_attrib.x1
            else:
                gene_attrib.arrow_size = 5

            gene_attrib.arrow_size = -gene_attrib.arrow_size if gene_attrib.direction == "f" else gene_attrib.arrow_size

            gene_attrib.x1 = max(chr_attrib.left, min(chr_attrib.right, gene_attrib.x1))
            gene_attrib.x2 = max(chr_attrib.left, min(chr_attrib.right, gene_attrib.x2))

            gene_attrib.y1 = gene_y
            gene_attrib.y2 = gene_y + gene_height

            if math.fabs(gene_attrib.x2 - gene_attrib.x1) > len(gene.wid) * 5:
                gene_attrib.label = gene.wid
            else:
                gene_attrib.label = ''

            gene_attrib.title = (gene.name or gene.wid).replace("'", "\'")

            gene_attrib.url = gene.get_absolute_url()

            gene_attrib.color = transcription_unit_index % num_colors

            genes.append(gene_attrib)

        #promoters
        def draw_segment(wid, coordinate, item_length, tip_title, typ, url):
            feature_attrib = attrib_class()

            feature_attrib.wid = wid
            feature_attrib.coordinate = coordinate
            feature_attrib.length = item_length

            feature_attrib.x = chr_attrib.left + float(feature_attrib.coordinate - chr_attrib.start) / chr_attrib.length * chr_attrib.width
            feature_attrib.width = chr_attrib.left + float(feature_attrib.coordinate + feature_attrib.length - 1 - chr_attrib.start) / chr_attrib.length * chr_attrib.width

            feature_attrib.x = max(chr_attrib.left, min(chr_attrib.right, feature_attrib.x))
            feature_attrib.width = max(chr_attrib.left, min(chr_attrib.right, feature_attrib.width)) - feature_attrib.x

            feature_attrib.title = tip_title.replace("'", "\'")

            if isinstance(typ, str):
                pass
            else:
                if typ and all_features_types_count[typ.id] == 0:
                    types.append(typ)
                    all_features_types_count[typ.id] = 1
                else:
                    all_features_types_count[None] = 1
                feature_attrib.color = all_feature_type_pks.index(typ.pk if typ else None) % num_colors
                if typ:
                    typ.color = feature_attrib.color
                feature_attrib.type = typ

            feature_attrib.url = url
            new_item = [feature_attrib.x, feature_attrib.x + feature_attrib.width]
            inserted = False
            y = 0

            for i, row in enumerate(feature_draw):
                if any(overlaps(item, new_item, 5) for item in row):
                    continue
                # Space left -> insert
                row.append(new_item)
                inserted = True
                feature_attrib.y = feature_y + i * (feature_height + 2)
                break

            if not inserted:
                # Create new row
                feature_attrib.y = feature_y + len(feature_draw) * (feature_height + 2)
                feature_draw.append([new_item])

            feature_attrib.height = y + feature_height

            return feature_attrib

        for transcription_unit in transcription_units:
            if transcription_unit.promoter_35_coordinate is not None:
                tu_coordinate = transcription_unit.get_coordinate() + transcription_unit.promoter_35_coordinate
                tu_length = transcription_unit.promoter_35_length

                if not tu_coordinate > chr_attrib.end or tu_coordinate + tu_length - 1 < chr_attrib.start:
                    tip_title = transcription_unit.name or transcription_unit.wid
                    url = transcription_unit.get_absolute_url()

                    promoters.append(draw_segment(transcription_unit.wid, tu_coordinate, tu_length, tip_title, 'Promoter -35 box', url))

            if transcription_unit.promoter_10_coordinate is not None:
                tu_coordinate = transcription_unit.get_coordinate() + transcription_unit.promoter_10_coordinate
                tu_length = transcription_unit.promoter_10_length

                if not (tu_coordinate > chr_attrib.end or tu_coordinate + tu_length - 1 < chr_attrib.start):
                    tip_title = transcription_unit.name or transcription_unit.wid
                    url = transcription_unit.get_absolute_url()

                    promoters.append(draw_segment(transcription_unit.wid, tu_coordinate, tu_length, tip_title, 'Promoter -10 box', url))

            for tr in transcription_unit.transcriptional_regulations.all():
                if tr.binding_site is not None and not (tr.binding_site.coordinate > chr_attrib.end or tr.binding_site.coordinate + tr.binding_site.length - 1 < chr_attrib.start):
                    tip_title = transcription_unit.name or transcription_unit.wid
                    url = transcription_unit.get_absolute_url()

                    tf_sites.append(draw_segment(transcription_unit.wid, tr.binding_site.coordinate, tr.binding_site.length, tip_title, 'Transcription factor binding site', url))

        feature_values = self.features.filter(coordinate__lte=chr_attrib.end, coordinate__gte=chr_attrib.start + 1 - F("length")).prefetch_related('chromosome_feature', 'chromosome_feature__type').distinct()

        for feature in feature_values:
            if feature.coordinate > chr_attrib.end or feature.coordinate + feature.length - 1 < chr_attrib.start:
                continue

            tip_title = feature.chromosome_feature.name or feature.chromosome_feature.wid
            url = feature.chromosome_feature.get_absolute_url()

            if feature.chromosome_feature.type.count() > 0:
                type_ = feature.chromosome_feature.type.all()[0]
            else:
                type_ = None

            features.append(draw_segment(feature.chromosome_feature.wid, feature.coordinate, feature.length, tip_title, type_, url))

        H = 2 + gene_height + 2 + 4 + 1 * (2 + len(feature_draw) * (feature_height + 2)) + 2

        c = Context({
            'species': self.species,
            'genes': genes,
            'height': H,
            'chromosomes': [chr_attrib],
            'features': features,
            'promoters': promoters,
            'tf_sites': tf_sites,
            'highlight_wid': highlight_wid,
            'types': types
        })

        template = loader.get_template("cyano/fields/structure.html")
        rendered = template.render(c)

        return rendered

    def get_as_html_genes(self, is_user_anonymous):
        from cyano.helpers import format_list_html_url
        return format_list_html_url(self.genes.all())

    def get_as_html_features(self, is_user_anonymous):
        from cyano.helpers import format_list_html_url

        features = self.features.all().values_list("chromosome_feature__pk", flat=True)
        
        return format_list_html_url(ChromosomeFeature.objects.for_species(self.species).filter(pk__in=features).distinct())

    def get_as_html_structure_filter(self, is_user_anonymous, start_coordinate=None, end_coordinate=None, zoom=0):
        from itertools import groupby

        if zoom == 0:
            types = Type.objects.filter(pk__in=self.features.values_list("chromosome_feature__type", flat=True)).distinct().order_by("parent")
        else:
            types = Type.objects.filter(pk__in=self.features.filter(coordinate__lte=end_coordinate, coordinate__gte=start_coordinate + 1 - F("length")).values_list("chromosome_feature__type", flat=True)).distinct().order_by("parent")

        if types.count() == 0:
            return ""

        # group by feature parent
        groups = {}
        i = 1
        for k, g in groupby(types, lambda x: x.parent):
            groups[k] = []
            for e in g:
                groups[k].append([i, e])
                i = (i+1)%20

        template = loader.get_template("cyano/fields/structure_filter.html")
        c = Context({'groups': groups})
        rendered = template.render(c)

        return rendered

    def get_as_fasta(self):
        return self.get_fasta_header() + "\r\n" + re.sub(r"(.{70})", r"\1\r\n", self.get_sequence()) + "\r\n"
    
    def get_as_genbank(self):
        genbank = StringIO()
        genes = Gene.objects.filter(species=self.species, chromosome_id=self.pk).prefetch_related("cross_references", "protein_monomers")
        record = SeqRecord.SeqRecord(Seq(self.sequence, IUPAC.IUPACAmbiguousDNA()))

        record.description = self.name
        record.name = self.wid
        accession = self.cross_references.filter(source="RefSeq")
        if len(accession) > 0:
            record.annotations["accession"] = accession[0].xid

        record.annotations["date"] = self.last_revision().detail.date.strftime("%d-%b-%Y").upper()
        record.annotations["source"] = self.species.name
        record.annotations["organism"] = self.species.name
        record.annotations["comment"] = self.comments
        
        features = record.features

        source = SeqFeature(FeatureLocation(0, self.length - 1), type="source")
        source.qualifiers["organism"] = self.species.name

        features += [source]

        for item in genes:
            features += item.get_as_seqfeature(record.seq)
        
        SeqIO.write(record, genbank, "genbank")
        
        return genbank.getvalue()

    def clean(self, all_obj_data=None, all_obj_data_by_model=None):
        if self.sequence is not None and self.sequence != '' and len(self.sequence) != self.length:
            raise ValidationError({'length': 'Length of sequence property must match length property'})

    #meta information
    class Meta:
        concrete_entry_model = True
        fieldsets = [
            ('Type', {'fields': ['model_type']}),
            ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}),
            ('Classification', {'fields': ['type']}),
            ('Sequence', {'fields': [
                {'verbose_name': 'Structure Filter', 'name': 'structure_filter'},
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
            'id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'sequence', 'length', 'comments', 'publication_references'
            ]
        facet_fields = ['type']
        verbose_name = 'Genome'
        verbose_name_plural = 'Genome'
        wid_unique = False

class Chromosome(Genome):
    #parent pointer
    parent_ptr_genome = OneToOneField(Genome, related_name='child_ptr_chromosome', parent_link=True, verbose_name='Genome')

    @classmethod
    def get_statistics(cls, species):
        # Done in Genome
        return None

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

    @classmethod
    def get_statistics(cls, species):
        # Done in Genome
        return None

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

    @classmethod
    def get_statistics(cls, species):
        url = cls.get_model_url(species)
        amount = cls.objects.for_species(species).count()

        if amount == 0:
            return None

        lst = [[0, cls._meta.verbose_name_plural, amount, None, url]]
        return lst

    #additional fields
    #positions = reverse relation
    def get_as_html_tooltip(self, is_user_anonymous):
        return ", ".join(map(lambda x: x.name or x.wid, self.type.filter(species=self.species)))

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
            'id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'comments', 'publication_references'
        ]
        facet_fields = ['type',] #ToDo 'positions.chromosome', 'positions.direction']
        verbose_name = 'Chromosome feature'
        verbose_name_plural = 'Chromosome features'
        wid_unique = False

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
            seq = str(Seq(seq, IUPAC.ambiguous_dna).reverse_complement())
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
    def get_as_html_structure(self, is_user_anonymous):
        return self.chromosome.get_as_html_structure(is_user_anonymous,
                                                     zoom=1,
                                                     start_coordinate=self.coordinate - 500,
                                                     end_coordinate=self.coordinate + self.length + 500,
                                                     highlight_wid=[self.chromosome_feature.wid])

    def get_as_html_sequence(self, is_user_anonymous):
        from cyano.helpers import format_sequence_as_html

        direction = CHOICES_DIRECTION[[x[0] for x in CHOICES_DIRECTION].index(self.direction)][1]

        return '%s: <a href="%s">%s</a>, Coordinate: %s (nt), Length: %s (nt), Direction: %s, Sequence: %s' % (
            self.chromosome.model_type.model_name,
            self.chromosome.get_absolute_url(), self.chromosome.wid,
            self.coordinate, self.length, direction,
            format_sequence_as_html(self.chromosome.species, self.get_sequence(), seq_offset=self.coordinate))

    def get_as_html_structure_filter(self, is_user_anonymous):
        return self.chromosome.get_as_html_structure_filter(is_user_anonymous,
                                                            zoom=1,
                                                            start_coordinate=self.coordinate - 500,
                                                            end_coordinate=self.coordinate + self.length + 500)

    def get_as_html_genes(self, is_user_anonymous):
        from cyano.helpers import format_list_html_url
        return format_list_html_url(self.get_genes())

    def get_as_html_transcription_units(self, is_user_anonymous):
        from cyano.helpers import format_list_html_url
        return format_list_html_url(self.get_transcription_units(), self.chromosome.species)

    def get_as_fasta(self):
        return self.get_fasta_header() + "\r\n" + re.sub(r"(.{70})", r"\1\r\n", self.get_sequence(cache=True)) + "\r\n"

    class Meta:
        verbose_name = "Feature Position"
        verbose_name_plural = "Feature Positions"
        fieldsets = [
            ('Structure', {'fields': [
                {'verbose_name': 'Structure', 'name': 'structure'},
                {'verbose_name': 'Structure Filter', 'name': 'structure_filter'},
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
    def get_protein_complexes(self):
        arr = []
        for obj in ProteinComplex.objects.for_species(self.species):
            if self.pk == obj.get_localization().pk:
                arr.append(obj)
        return arr

    #getters
    def get_as_html_biomass_compositions(self, is_user_anonymous):
        results = []
        for bm in self.biomass_compositions.all():
            if len(bm.metabolites.all()) == 0:
                continue
            m = bm.metabolites.all()[0]
            results.append('<a href="%s">%s</a>: %.4f' % (m.get_absolute_url(), m.name, bm.concentration))
        return format_list_html(results, comma_separated=False)

    def get_as_html_protein_monomers(self, is_user_anonymous):
        from cyano.helpers import format_list_html_url
        return format_list_html_url(self.protein_monomers.all())

    def get_as_html_protein_complexes(self, is_user_anonymous):
        from cyano.helpers import format_list_html_url
        return format_list_html_url(self.get_protein_complexes())

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
            'id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'comments', 'publication_references'
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
    codons = ManyToManyField(Codon, blank=True, related_name='genes', verbose_name='Codons')
    amino_acid = ForeignKey('Metabolite', blank=True, null=True, on_delete=SET_NULL, related_name='genes', verbose_name='Amino acid')
    homologs = ManyToManyField(Homolog, blank=True, related_name='genes', verbose_name='Homologs')

    #getters    
    def get_sequence(self, cache=False):
        if cache:
            chromosome_key = "chromosome/%s" % self.chromosome_id
            chromosome = Cache.try_get(chromosome_key, lambda: self.chromosome, 60)
        else:
            chromosome = self.chromosome
        seq = chromosome.sequence[self.coordinate - 1:self.coordinate - 1 + self.length]

        if self.direction == 'r':
            seq = Seq(seq, IUPAC.ambiguous_dna).reverse_complement()
        return str(seq)

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

        seq = Seq(str(self.get_sequence()), IUPAC.ambiguous_dna).transcribe()

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
    def get_as_html_codons(self, is_user_anonymous):
        result = []
        for codon in self.codons.all():
            result.append(format_with_evidence(list_item = True, obj = codon, txt = '<tt>%s</tt>' % codon.sequence))
        return format_list_html(result, force_list=True)

    def get_as_html_homologs(self, is_user_anonymous):
        results = []
        for h in self.homologs.all():
            results.append(format_with_evidence(list_item = True, obj = h, txt = '%s: <a href="%s">%s</a>' % (h.species, HOMOLOG_SPECIES_URLS[h.species] % h.xid, h.xid)))
        return format_list_html(results, force_list=True)

    def get_as_html_structure(self, is_user_anonymous):
        return self.chromosome.get_as_html_structure(is_user_anonymous,
            zoom = 1,
            start_coordinate = self.coordinate - 2500,
            end_coordinate = self.coordinate + self.length + 2500,
            highlight_wid = [self.wid])

    def get_as_html_sequence(self, is_user_anonymous):
        from cyano.helpers import format_sequence_as_html

        direction = CHOICES_DIRECTION[[x[0] for x in CHOICES_DIRECTION].index(self.direction)][1]

        return '%s: <a href="%s">%s</a>, Coordinate: %s (nt), Length: %s (nt), Direction: %s, G/C content: %.1f%%, Sequence: %s' % (
            self.chromosome.model_type.model_name,
            self.chromosome.get_absolute_url(), self.chromosome.wid,
            self.coordinate, self.length, direction,
            self.get_gc_content() * 100,
            format_sequence_as_html(self.species, self.get_sequence(), seq_offset=self.coordinate, show_protein_seq=True))

    def get_as_html_structure_filter(self, is_user_anonymous):
        return self.chromosome.get_as_html_structure_filter(is_user_anonymous,
            zoom = 1,
            start_coordinate = self.coordinate - 2500,
            end_coordinate = self.coordinate + self.length + 2500,
            )

    def get_as_fasta(self):
        return self.get_fasta_header() + "\r\n" + re.sub(r"(.{70})", r"\1\r\n", self.get_sequence(cache=True)) + "\r\n"
    
    def get_as_seqfeature(self, sequence):
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
        cds.qualifiers["transl_table"] = [self.species.genetic_code]
        cds.qualifiers["codon_start"] = [1]

        typ = list(self.type.all()[:1])
        cds.type = typ[0].wid
        if cds.type == "mRNA":
            cds.type = "CDS"
            s = gene.extract(sequence)
            cds.qualifiers["translation"] = [s.translate(table=self.species.genetic_code)[:-1]]

        #monomer = ProteinMonomer.objects.values_list('wid', 'name').filter(species = species).get(gene_id = self.pk)
        monomer = list(self.protein_monomers.all()[:1])
        if len(monomer) > 0:
            cds.qualifiers["protein_id"] = monomer[0].wid
            cds.qualifiers["product"] = monomer[0].name

        return gene, cds

    def get_as_html_tooltip(self, is_user_anonymous):
        from cyano.helpers import format_field_detail_view
        return "Transcription unit: %s" % (format_field_detail_view(self, "transcription_units", is_user_anonymous))

    def clean(self, all_obj_data=None, all_obj_data_by_model=None):
        if all_obj_data is None:
            chro = Genome.objects.get(species__wid=self.species, wid=self.chromosome)
        else:
            chro = all_obj_data[self.chromosome]

        if isinstance(chro, Entry):
            chr_len = chro.length
        else:
            chr_len = chro['length']

        if self.coordinate > chr_len:
            raise ValidationError({'coordinate': 'Coordinate must be less then chromosome length.'})
        if self.length > chr_len:
            raise ValidationError({'length': 'Length must be less then chromosome length.'})

    #meta information
    class Meta:
        concrete_entry_model = True
        fieldsets = [
            ('Type', {'fields': ['model_type']}),
            ('Name', {'fields': ['wid', 'name', 'symbol', 'synonyms', {'verbose_name': 'Protein product', 'name': 'protein_monomers'}, 'cross_references', 'homologs']}),
            ('Classification', {'fields': ['type']}),
            ('Structure', {'fields': [
                {'verbose_name': 'Structure', 'name': 'structure'},
                {'verbose_name': 'Structure Filter', 'name': 'structure_filter'},
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
    biomass_composition = ManyToManyField(BiomassComposition, blank=True, related_name='metabolites', verbose_name='Biomass composition (SP4 media, <br/>5% CO<sub>2</sub>, 37C; mmol gDCW<sup>-1</sup>)')
    media_composition = ForeignKey(MediaComposition, blank=True, null=True, on_delete=SET_NULL, related_name='metabolites', verbose_name='Media composition (SP4; mM)')
    map_coordinates = ManyToManyField(MetaboliteMapCoordinate, blank=True, related_name='metabolites', verbose_name='Map coordinates')

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
    def get_as_html_structure(self, is_user_anonymous):
        from cyano.helpers import draw_molecule
        return draw_molecule(self.smiles, 'svg', 636, 150)

    def get_as_html_empirical_formula(self, is_user_anonymous):
        from cyano.helpers import EmpiricalFormula
        return EmpiricalFormula(self.empirical_formula).get_as_html()

    def get_as_html_biomass_composition(self, is_user_anonymous):
        results = []
        for b in self.biomass_composition.all():
            results.append(format_with_evidence(list_item = True, obj = b, txt = '%s [<a href="%s">%s</a>]' % (b.concentration, b.compartment.get_absolute_url(), b.compartment.wid)))
        return format_list_html(results, force_list=True)

    def get_as_html_media_composition(self, is_user_anonymous):
        m = self.media_composition
        if m is None:
            return
        if m.is_diffused:
            txt = '%s (diffused)' % (m.concentration, )
        else:
            txt = m.concentration

        return format_with_evidence(obj = m, txt = txt)

    def get_as_html_coenzyme_participants(self, is_user_anonymous):
        results = []
        for obj in self.coenzyme_participants.all():
            if len(obj.reactions.all()) == 0:
                continue
            rxn = obj.reactions.all()[0]
            results.append('<a href="%s">%s</a><br/>%s' % (rxn.get_absolute_url(), rxn.name, rxn.get_as_html_stoichiometry(is_user_anonymous)))
        return format_list_html(results, vertical_spacing=True)

    def get_as_html_prosthetic_group_participants(self, is_user_anonymous):
        results = []
        for obj in self.prosthetic_group_participants.all():
            if len(obj.proteins.all()) == 0:
                continue
            p = obj.proteins.all()[0]
            results.append('<a href="%s">%s</a>' % (p.get_absolute_url(), p.name))
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
            'publication_references'
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
    reactions = ManyToManyField('Reaction', blank=True, related_name='parameters', verbose_name='Reactions')
    molecules = ManyToManyField('Molecule', blank=True, related_name='parameters', verbose_name='Molecules')
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

        # Create new boehringer (if it was never created yet) or assign to one
        x, created = Pathway.objects.for_species(species).for_wid(wid, create=True, creation_status=True)
        if created:
            typ = Type.objects.for_wid("Main-Pathway", create=True)
            typ.species = species
            typ.save(revdetail)
            x.name = "Biochemical Pathways"
            x.species = species
            x.save(revdetail)
            x.type.add(typ)

        x.species = species

    @staticmethod
    def add_kegg_pathway(species, revdetail):
        from kegg.models import Map

        crs = CrossReference.objects.filter(cross_referenced_components__species = species, source = "EC").values_list('xid', flat=True)
        maps = Map.objects.filter(ec_numbers__name__in=crs).distinct()

        for map_ in maps:
            x = Pathway.objects.for_species(species).for_wid(map_.name, create=True)
            x.name = map_.title
            x.species = species
            x.save(revdetail)
            if map_.overview:
                typ = Type.objects.for_wid("Main-Pathway", create=True)
                typ.species = species
                typ.save(revdetail)
                x.type.add(typ)

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

    def get_boehringer_hits(self):
        import boehringer.models as bmodels
        species_ecs = CrossReference.objects.filter(cross_referenced_components__species = self.species, source = "EC").values_list('xid', flat=True)
        enzymes = bmodels.Enzyme.objects.filter(ec__in = species_ecs).order_by('title')

        species_metabolites = Metabolite.objects.for_species(self.species).values_list('name', flat=True)
        metabolites = bmodels.Metabolite.objects.filter(title__in = species_metabolites).order_by('title')

        return enzymes, metabolites

    # TODO: The current logic hardcodes on Boehringer and uses KEGG otherwise
    # Not really nice solution...

    def get_as_html_navigator(self, is_user_anonymous):
        if self.wid != "Boehringer":
            return ""

        enzymes, metabolites = self.get_boehringer_hits()

        template = loader.get_template("cyano/pathway/navigator.html")

        c = Context({'enzymes': enzymes, 'metabolites': metabolites})

        rendered = template.render(c)
        return rendered

    #html formatting
    def get_as_html_reaction_map(self, is_user_anonymous):
        #W = 731
        H = 600

        if self.wid == "Boehringer":
            enzymes, metabolites = self.get_boehringer_hits()
            template = loader.get_template("cyano/pathway/reaction_map_boehringer.html")
            c = Context({'enzymes': enzymes, 'metabolites': metabolites,
                         'width': '100%', 'height': H})
            return template.render(c)

        # Wordaround broken PIL installation
        # see: https://code.djangoproject.com/ticket/6054
        try:
            from PIL import Image
        except ImportError:
            import Image

        from xml.etree.ElementTree import ElementTree, Element

        # All Pathways of the organism
        species_ecs = CrossReference.objects.filter(cross_referenced_components__species = self.species, source = "EC").values_list('xid', flat=True)

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
            root.set("width", str("100%"))
            root.set("height", str(min(iheight,H)))
            root.set("viewport", "0 0 {} {}".format("100%", str(min(iheight,H))))

            script = Element("script")
            script.set("xlink:href", "{}boehringer/js/svg-pan-zoom.js".format(settings.STATIC_URL))
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
                        pw_obj = Pathway.objects.for_species(self.species).for_wid(pathway_name)
                        elem.set("xlink:href", pw_obj.get_absolute_url())
                        fill_opacity = "0.3"
                        fill_color = "blue"
                    except ObjectDoesNotExist:
                        elem.set("xlink:href", reverse("kegg.views.map_view", kwargs={"map_id": pathway_name}))
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
                            fill_opacity = "0.3"

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

            out = StringIO()
            out.write(template.render(Context()))

            out.write('<script type="text/javascript" src="' + settings.STATIC_URL + 'boehringer/js/svg-pan-zoom.js"></script>')
            xmlout = BytesIO()
            tree.write(xmlout)
            out.write(xmlout.getvalue().decode())
            out.write('<script type="text/javascript">var svgPan = svgPanZoom("svg", {minZoom: 0.1});</script>')

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
            'comments',
            'publication_references',
            'type'
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
    def get_as_html_reactions(self, is_user_anonymous):
        results = []
        for reaction in self.reactions.all():
            results.append('<a href="%s">%s</a><br/>%s' % (reaction.get_absolute_url(), reaction.name, reaction.get_as_html_stoichiometry(is_user_anonymous)))
        return format_list_html(results, vertical_spacing=True)

    def get_as_html_formed_complexes(self, is_user_anonymous):
        results = []
        for complexe in self.formed_complexes.all():
            results.append('<a href="%s">%s</a><br/>%s' % (complexe.get_absolute_url(), complexe.name, complexe.get_as_html_biosynthesis(is_user_anonymous)))
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
            ]
        facet_fields = ['type']
        verbose_name='Process'
        verbose_name_plural = 'Processes'
        wid_unique = False

class ProteinComplex(Protein):
    #parent pointer
    parent_ptr_protein = OneToOneField(Protein, related_name='child_ptr_protein_complex', parent_link=True, verbose_name='Protein')

    #additional fields
    biosynthesis = ManyToManyField(ProteinComplexBiosythesisParticipant, related_name='protein_complexes', verbose_name='Biosynthesis')
    disulfide_bonds = ManyToManyField(DisulfideBond, blank=True, related_name='protein_complexes', verbose_name='Disulfide bonds (pH 7.5)')
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

    def get_as_html_num_subunits(self, is_user_anonymous):
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
    def get_as_html_biosynthesis(self, is_user_anonymous):
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
                tmp += '<a href="%s">%s</a>' % (s.molecule.get_absolute_url(), s.molecule.wid)
                if len(compartments) > 1:
                    tmp += '[<a href="%s">%s</a>]' % (s.compartment.get_absolute_url(), s.compartment.wid)
                pos.append(tmp)
            else:
                tmp = ''
                if s.coefficient != 1:
                    tmp += '(%d) ' % s.coefficient
                tmp += '<a href="%s">%s</a>' % (s.molecule.get_absolute_url(), s.molecule.wid)
                if len(compartments) > 1:
                    tmp += '[<a href="%s">%s</a>]' % (s.compartment.get_absolute_url(), s.compartment.wid)
                neg.append(tmp)

        result = ''
        if len(compartments) == 1:
            result += '[<a href="%s">%s</a>]: ' % (compartments[0].get_absolute_url(), compartments[0].wid)
        result += ' + '.join(pos)
        result += ' &rArr; '
        result += ' + '.join(neg)
        return format_with_evidence(obj = self.biosynthesis.all(), txt = result)

    def get_as_html_disulfide_bonds(self, is_user_anonymous):
        results = [];
        for b in self.disulfide_bonds.all():
            results.append(format_with_evidence(list_item = True, obj = b, txt = '<a href="%s">%s</a>: %s-%s' % (b.protein_monomer.get_absolute_url(), b.protein_monomer.wid, b.residue_1, b.residue_2)))
        return format_list_html(results, force_list=True)

    def get_as_html_localization(self, is_user_anonymous):
        localization = self.get_localization()
        return '<a href="%s">%s</a>' % (localization.get_absolute_url(), localization.wid, )

    def clean(self, all_obj_data=None, all_obj_data_by_model=None):
        from cyano.helpers import getModel, getEntry

        #biosynthesis
        coeff = 0
        for b in self.biosynthesis:
            if b['molecule'] == self.wid:
                coeff += b['coefficient']

        if coeff != 1:
            raise ValidationError({'biosynthesis': 'Protein complex must appear on the right side of the biosynthesis reaction'})

        #disulfide bonds
        for dsfb in self.disulfide_bonds:
            mon_wid = dsfb['protein_monomer']

            if all_obj_data is None:
                mon = ProteinMonomer.objects.get(species__wid=self.species, wid=mon_wid)
            else:
                mon = all_obj_data[mon_wid]
            if isinstance(mon, Entry):
                gene_wid = mon.gene.wid
            else:
                gene_wid = mon['gene']

            if all_obj_data is None:
                gene = Gene.objects.get(species_wid=self.species, wid=gene_wid)
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
        for b in self.biosynthesis:
            if all_obj_data is None:
                molecule = getEntry(species_wid=self.species, wid=b['molecule'])
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
                    gene = Gene.objects.get(species__wid=self.species, wid=gene_wid)
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
            ]
        facet_fields = ['type', 'dna_footprint__binding', 'dna_footprint__region', 'formation_process', 'chaperones']
        verbose_name='Protein complex'
        verbose_name_plural = 'Protein complexes'
        wid_unique = False

class ProteinMonomer(Protein):
    #parent pointer
    parent_ptr_protein = OneToOneField(Protein, related_name='child_ptr_protein_monomer', parent_link=True, verbose_name='Protein')

    #additional fields
    gene = ForeignKey(Gene, related_name='protein_monomers', verbose_name='Gene')
    is_n_terminal_methionine_cleaved = ForeignKey(EntryBooleanData, null = True, verbose_name='Is N-terminal methionine cleaved', related_name='+')
    localization = ForeignKey(Compartment, null = True, related_name='protein_monomers', verbose_name='Localization')
    signal_sequence = ForeignKey(SignalSequence, blank=True, null=True, related_name='protein_monomers', on_delete=SET_NULL, verbose_name='Sequence sequence')

    #getters
    def get_sequence(self, cache=False):
        return str(Seq(str(self.gene.get_sequence(cache=cache)), IUPAC.ambiguous_dna).translate(table=self.species.genetic_code))

    def get_length(self):
        return len(self.get_sequence())

    def get_neg_aa(self):
        seq = self.get_sequence
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
        return value

    def get_half_life(self):
        return 20 * 60

    #http://ca.expasy.org/tools/protparam-doc.html
    def get_instability(self):
        from cyano.helpers import DipeptideInstabilityWeight

        seq = self.get_sequence()
        value = 0.
        for i in range(len(seq)-1):
            if seq[i] != '*' and seq[i+1] != '*':
                value += DipeptideInstabilityWeight.value[seq[i]][seq[i+1]]
        return 10. / float(len(seq)) * value

    #http://ca.expasy.org/tools/protparam-doc.html
    def get_is_stable(self):
        return self.get_instability() < 40.

    #http://ca.expasy.org/tools/protparam-doc.html
    def get_aliphatic(self):
        seq = self.get_sequence()
        return 100. * (
            + 1.0 * float(seq.count('A'))
            + 2.9 * float(seq.count('V'))
            + 3.9 * float(seq.count('I'))
            + 3.9 * float(seq.count('L'))
                      ) / float(len(seq))

    #http://ca.expasy.org/tools/protparam-doc.html
    def get_gravy(self):
        seq = self.get_sequence()
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
    def get_extinction_coefficient(self):
        seq = self.get_sequence()
        return \
            + seq.count('W') * 5500 \
            + seq.count('Y') * 1490 \
            + seq.count('C') * 125

    #html formatting
    def get_as_html_sequence(self, is_user_anonymous):
        from cyano.helpers import format_sequence_as_html
        return format_sequence_as_html(self.species, self.get_sequence())

    def get_as_html_signal_sequence(self, is_user_anonymous):
        ss = self.signal_sequence
        if ss is None:
            return
        return format_with_evidence(obj = ss, txt = 'Type: %s, Location: %s, Length: %s (nt)' % (ss.type, ss.location, ss.length))

    def get_as_html_disulfide_bonds(self, is_user_anonymous):
        results = []
        for b in self.disulfide_bonds.all():
            results.append('<a href="%s">%s</a>: %s-%s' % (b.protein_complexes.all()[0].get_absolute_url(), b.protein_complexes.all()[0].wid, b.residue_1, b.residue_2))
        return format_list_html(results)

    def get_as_html_instability(self, is_user_anonymous):
        return self.get_instability()

    def get_as_html_is_stable(self, is_user_anonymous):
        return self.get_is_stable()

    def get_as_html_aliphatic(self, is_user_anonymous):
        return self.get_aliphatic()

    def get_as_html_gravy(self, is_user_anonymous):
        return self.get_gravy()

    def get_as_html_interactions(self, is_user_anonymous):
        from cyanointeraction.models import ProteinsNames as IproteinsNames

        prot_name = IproteinsNames.objects.filter(protein_name=self.gene.wid).first()

        if prot_name:
            return '<a href="{}">Show interactions</a>'.format(
                reverse("cyanointeraction.views.checkInteraction", kwargs={"protID": prot_name.protein_id}))

        return ""

    def get_as_fasta(self):
        return self.get_fasta_header() + "\r\n" + re.sub(r"(.{70})", r"\1\r\n", self.get_sequence(cache=True)) + "\r\n"

    def get_as_genbank(self):
        genbank = StringIO()

        record = SeqRecord.SeqRecord(Seq(self.get_sequence()[:-1], IUPAC.IUPACProtein()))
        record.description = self.name
        record.name = self.wid
        accession = self.cross_references.filter(source="RefSeq")
        if len(accession) > 0:
            record.annotations["accession"] = accession[0].xid

        record.annotations["date"] = self.last_revision().detail.date.strftime("%d-%b-%Y").upper()
        record.annotations["source"] = self.species.name
        record.annotations["organism"] = self.species.name
        record.annotations["comment"] = self.comments

        source = SeqFeature(FeatureLocation(0, self.get_length() - 1), type="source")
        source.qualifiers["organism"] = self.species.name

        record.features += [source]
        record.features += self.get_as_seqfeature()

        SeqIO.write(record, genbank, "genbank")

        return genbank.getvalue()

    def get_as_seqfeature(self):
        gene = SeqFeature(FeatureLocation(0, self.get_length() - 1), type="Protein")
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
        cds.qualifiers["transl_table"] = [self.species.genetic_code]
        cds.qualifiers["codon_start"] = [1]

        cds.type = "CDS"
        #if cds.type == "mRNA":
        #    cds.type = "CDS"
            #s = self.get_sequence(species)
            #cds.qualifiers["translation"] = [s.translate(table=species.genetic_code)[:-1]]

        return gene, cds

    def clean(self, all_obj_data=None, all_obj_data_by_model=None):
        if self.signal_sequence is not None:
            if all_obj_data is None:
                gene = Gene.objects.get(species__wid=self.species, wid=self.gene)
            else:
                gene = all_obj_data[self.gene]
            if isinstance(gene, Entry):
                mon_len = gene.length / 3
            else:
                mon_len = gene['length'] / 3

            if self.signal_sequence.length > mon_len:
                raise ValidationError({'signal_sequence': 'Length must be less than protein length'})

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
            ('Interactions', {'fields': [{'verbose_name': 'Protein/Metabolite interactions', 'name': 'interactions'}]}),
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
            ]
        facet_fields = ['type', 'is_n_terminal_methionine_cleaved__value', 'signal_sequence__type', 'signal_sequence__location', 'dna_footprint__binding', 'dna_footprint__region', 'localization', 'chaperones']
        verbose_name='Protein monomer'
        verbose_name_plural = 'Protein monomers'
        wid_unique = False

class Reaction(SpeciesComponent):
    #parent pointer
    parent_ptr_species_component = OneToOneField(SpeciesComponent, related_name='child_ptr_reaction', parent_link=True, verbose_name='Species component')

    #additional fields
    stoichiometry = ManyToManyField(ReactionStoichiometryParticipant, related_name='reactions', verbose_name='Stoichiometry')
    direction = CharField(max_length=1, choices=CHOICES_REACTION_DIRECTION, verbose_name='Direction')
    modification = ForeignKey(ModificationReaction, blank=True, null=True, on_delete=SET_NULL, related_name='reactions', verbose_name='Modification')
    enzyme = ForeignKey(EnzymeParticipant, blank=True, null=True, on_delete=SET_NULL, related_name='reactions', verbose_name='Enzyme')
    coenzymes = ManyToManyField(CoenzymeParticipant, blank=True, related_name='reactions', verbose_name='Coenzymes')
    is_spontaneous = BooleanField(verbose_name='Is spontaneous (pH 7.5, 25C, <i>I</i> = 0)')
    delta_g = FloatField(blank=True, null=True, verbose_name='&Delta;G (pH 7.5, 25C, <i>I</i> = 0; kJ mol<sup>-1</sup>)')
    keq = ForeignKey(EntryPositiveFloatData, blank=True, null=True, on_delete=SET_NULL, verbose_name='K<sub>eq</sub>', related_name='+')
    kinetics_forward = ForeignKey(Kinetics, blank=True, null=True, on_delete=SET_NULL, related_name='reactions_forward', verbose_name='Forward kinetics')
    kinetics_backward = ForeignKey(Kinetics, blank=True, null=True, on_delete=SET_NULL, related_name='reactions_backward', verbose_name='Backward kinetics')
    optimal_ph = ForeignKey(EntryPositiveFloatData, blank=True, null=True, on_delete=SET_NULL, verbose_name='Optimal pH', related_name='+')
    optimal_temperature = ForeignKey(EntryFloatData, blank=True, null=True, on_delete=SET_NULL, verbose_name='Optimal temperature', related_name='+')
    pathways = ManyToManyField('Pathway', blank=True, related_name='reactions', verbose_name='Pathways')
    processes = ForeignKey('Process', blank=True, null=True, on_delete=SET_NULL, related_name='reactions', verbose_name='Process')
    states = ForeignKey('State', blank=True, null=True, on_delete=SET_NULL, related_name='reactions', verbose_name='State')
    map_coordinates = ManyToManyField(ReactionMapCoordinate, blank=True, related_name='reactions', verbose_name='Map coordinates')

    #getters

    #html formatting
    def get_as_html_stoichiometry(self, is_user_anonymous):
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
                tmp += '<a href="%s">%s</a>' % (s.molecule.get_absolute_url(), s.molecule.name)
                if len(compartments) > 1:
                    tmp += '[<a href="%s">%s</a>]' % (s.compartment.get_absolute_url(), s.compartment.wid)
                pos.append(tmp)
            else:
                tmp = ''
                if s.coefficient != 1:
                    tmp += '(%d) ' % s.coefficient
                tmp += '<a href="%s">%s</a>' % (s.molecule.get_absolute_url(), s.molecule.name)
                if len(compartments) > 1:
                    tmp += '[<a href="%s">%s</a>]' % (s.compartment.get_absolute_url(), s.compartment.wid)
                neg.append(tmp)

        result = ''
        if len(compartments) == 1:
            result += '[<a href="%s">%s</a>]: ' % (compartments[0].get_absolute_url(), compartments[0].wid)
        result += ' + '.join(pos)
        if self.direction == 'f':
            result += ' &rArr; '
        elif self.direction == 'b':
            result += ' &lArr; '
        elif self.direction == 'r':
            result += ' &hArr; '
        result += ' + '.join(neg)
        return format_with_evidence(obj = self.stoichiometry.all(), txt = result)

    def get_as_html_modification(self, is_user_anonymous):
        m = self.modification
        if m is None:
            return
        if m.position is None:
            txt = '<a href="%s">%s</a> [<a href="%s">%s</a>]' % (m.molecule.get_absolute_url(), m.molecule.wid, m.compartment.get_absolute_url(), m.compartment.wid)
        else:
            txt = '(%d) <a href="%s">%s</a> [<a href="%s">%s</a>]' % (m.position, m.molecule.get_absolute_url(), m.molecule.wid, m.compartment.get_absolute_url(), m.compartment.wid)

        return format_with_evidence(obj = m, txt = txt)

    def get_as_html_enzyme(self, is_user_anonymous):
        e = self.enzyme
        if e is None:
            return
        return format_with_evidence(obj = e, txt = '<a href="%s">%s</a> [<a href="%s">%s</a>]' % (e.protein.get_absolute_url(), e.protein.wid, e.compartment.get_absolute_url(), e.compartment.wid))

    def get_as_html_coenzymes(self, is_user_anonymous):
        results = []
        for c in self.coenzymes.all():
            if c.coefficient is None:
                results.append(format_with_evidence(list_item = True, obj = c, txt = '<a href="%s">%s</a> [<a href="%s">%s</a>]' % (c.metabolite.get_absolute_url(), c.metabolite.wid, c.compartment.get_absolute_url(), c.compartment.wid)))
            else:
                results.append(format_with_evidence(list_item = True, obj = c, txt = '(%d) <a href="%s">%s</a> [<a href="%s">%s</a>]' % (c.coefficient, c.metabolite.get_absolute_url(), c.metabolite.wid, c.compartment.get_absolute_url(), c.compartment.wid)))

        return format_list_html(results, force_list=True)

    def get_as_html_kinetics_forward(self, is_user_anonymous):
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

    def get_as_html_kinetics_backward(self, is_user_anonymous):
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

    def clean(self, all_obj_data=None, all_obj_data_by_model=None):
        from cyano.helpers import getEntry, getModel, EmpiricalFormula

        #stoichiometry
        formula = EmpiricalFormula()
        includes_macromolecules = False
        for s in self.stoichiometry:
            if all_obj_data is None:
                molecule = getEntry(species_wid=self.species, wid=s['molecule'])
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
            if self.modification is not None:
                includes_macromolecules = True

        if len(formula) > 0 and not includes_macromolecules: #todo: remove "includes_macromolecules"
            raise ValidationError({'stoichiometry': 'Reaction imbalanced by %s' % formula.get_as_html()})

        #kinetics
        if self.kinetics_forward    is not None:
            validate_kinetics(obj_data, 'f')
        if self.kinetics_backward is not None:
            validate_kinetics(obj_data, 'r')

        #modication
        if self.modification is not None:
            mod = self.modification

            if all_obj_data is None:
                molecule = getEntry(species_wid=self.species, wid=mod['molecule'])
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
                    gene = Gene.objects.get(species__wid=self.species, wid=gene_wid)
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
            ]
        facet_fields = ['type', 'direction', 'enzyme__protein', 'coenzymes__metabolite', 'is_spontaneous', 'pathways', 'processes', 'states']
        verbose_name='Reaction'
        verbose_name_plural = 'Reactions'
        wid_unique = False

class Species(Entry):
    #parent pointer
    parent_ptr_entry = OneToOneField(Entry, related_name='child_ptr_species', parent_link=True, verbose_name='Entry')

    #additional fields
    genetic_code = CharField(max_length=50, verbose_name='Genetic code', choices = CHOICES_GENETIC_CODE)

    cross_references = ManyToManyField("CrossReference", blank=True, related_name='cross_referenced_species', verbose_name='Cross references')
    publication_references = ManyToManyField("PublicationReference", blank=True, related_name='publication_referenced_species', verbose_name='Publication references')

    #getters
    @permalink
    def get_absolute_url(self, history_id = None):
        return ('cyano.views.species', (), {'species_wid': self.wid})

    @permalink
    def get_absolute_url_api(self):
        return ('cyano-api-species', (), {'species_wid': self.wid})

    #html formatting
    def get_as_html_comments(self, is_user_anonymous):
        txt = self.comments

        #provide links to references
        return re.sub(r'\[(PUB_\d{4,4})(, PUB_\d{4,4})*\]',
            lambda match: '[' + ', '.join(['<a href="%s">%s</a>' % (reverse('cyano.views.detail', kwargs={'species_wid':self.wid, 'wid': x}), x, ) for x in match.group(0)[1:-1].split(', ')]) + ']',
            txt)

    def get_as_html_genetic_code(self, is_user_anonymous):
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
            seq = str(Seq(seq, IUPAC.ambiguous_dna).reverse_complement())
        return seq

    def get_gc_content(self):
        seq = self.get_sequence()
        return float(seq.count('G') + seq.count('C')) / float(len(seq))

    #html formatting
    def get_as_html_structure(self, is_user_anonymous):
        return self.get_chromosome().get_as_html_structure(is_user_anonymous,
            zoom = 1,
            start_coordinate = self.get_coordinate() - 500,
            end_coordinate = self.get_coordinate() + self.get_length() + 500,
            highlight_wid = [self.wid] + [g.wid for g in self.genes.all()])

    def get_as_html_genes(self, is_user_anonymous):
        from cyano.helpers import format_list_html_url
        return format_list_html_url(self.genes.all())

    def get_as_html_sequence(self, is_user_anonymous):
        from cyano.helpers import format_sequence_as_html

        direction = CHOICES_DIRECTION[[x[0] for x in CHOICES_DIRECTION].index(self.get_direction())][1]

        return '%s: <a href="%s">%s</a>, Coordinate: %s, Length: %s, Direction: %s, Sequence: %s' % (
            self.get_chromosome().model_type.model_name,
            self.get_chromosome().get_absolute_url(), self.get_chromosome().wid,
            self.get_coordinate(), self.get_length(), direction,
            format_sequence_as_html(self.species, self.get_sequence()))

    def get_as_html_structure_filter(self, is_user_anonymous):
        return self.get_chromosome().get_as_html_structure_filter(is_user_anonymous,
            zoom = 1,
            start_coordinate = self.get_coordinate() - 500,
            end_coordinate = self.get_coordinate() + self.get_length() + 500)

    def get_as_html_transcriptional_regulations(self, is_user_anonymous):
        results = []
        for r in self.transcriptional_regulations.all():
            results.append('<a href="%s">%s</a>: <a href="%s">%s</a>' % (r.get_absolute_url(), r.wid, r.transcription_factor.get_absolute_url(), r.transcription_factor.wid))
        return format_list_html(results)
    
    def get_as_fasta(self):
        return self.get_fasta_header() + "\r\n" + re.sub(r"(.{70})", r"\1\r\n", self.get_sequence(cache=True)) + "\r\n"

    def clean(self, all_obj_data=None, all_obj_data_by_model=None):
        if len(self.genes) == 1:
            return
        if len(self.genes) == 0:
            raise ValidationError({'genes': 'Transcription units most contain at least 1 gene'})

        chr_wids = []
        if all_obj_data is None:
            for gene_wid in self.genes:
                chr_wids.append(Gene.objects.get(species__wid=self.species, wid=gene_wid).chromosome.wid)
        else:
            for gene_wid in self.genes:
                if isinstance(all_obj_data[gene_wid], Entry):
                    chr_wids.append(all_obj_data[gene_wid].chromosome.wid)
                else:
                    chr_wids.append(all_obj_data[gene_wid]['chromosome'])

        if len(set(chr_wids)) > 1:
            raise ValidationError({'genes': 'Genes must all belong to the same chromosome'})

    #meta information
    class Meta:
        concrete_entry_model = True
        fieldsets = [
            ('Type', {'fields': ['model_type']}),
            ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}),
            ('Classification', {'fields': ['type']}),
            ('Structure (Hayflick media, 37C)', {'fields': [
                {'verbose_name': 'Structure', 'name': 'structure'},
                {'verbose_name': 'Structure Filter', 'name': 'structure_filter'},
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
            'publication_references'
            ]
        facet_fields = ['type']
        verbose_name='Transcription unit'
        verbose_name_plural = 'Transcription units'
        wid_unique = False

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
    def get_as_html_binding_site(self, is_user_anonymous):
        from cyano.helpers import format_sequence_as_html

        bs = self.binding_site
        if bs is None:
            return None

        direction = CHOICES_DIRECTION[[x[0] for x in CHOICES_DIRECTION].index(bs.direction)][1]

        chro = self.transcription_unit.get_chromosome()

        structure = chro.get_as_html_structure(is_user_anonymous,
                zoom = 1,
                start_coordinate = bs.coordinate - 500,
                end_coordinate = bs.coordinate + bs.length + 500,
                highlight_wid = [self.wid])

        txt = '%s<br/>%s: <a href="%s">%s</a>, Coordinate: %s (nt), Length: %s (nt), Direction: %s, Sequence: %s' % (
            chro.model_type.model_name,
            structure, chro.get_absolute_url(), chro.wid,
            bs.coordinate, bs.length, direction,
            format_sequence_as_html(self.species, self.get_binding_site_sequence(), show_protein_seq=True))

        return format_with_evidence(obj = bs, txt = txt)

    def clean(self, all_obj_data=None, all_obj_data_by_model=None):
        #gene wid
        if all_obj_data is None:
            tu = TranscriptionUnit.objects.get(species__wid=self.species, wid=self.transcription_unit)
        else:
            tu = all_obj_data[self.transcription_unit]
        if isinstance(tu, Entry):
            gene_wid = tu.genes.all()[0].wid
        else:
            gene_wid = tu['genes'][0]

        #chr wid
        if all_obj_data is None:
            gene = Gene.objects.get(species__wid=self.species, wid=gene_wid)
        else:
            gene = all_obj_data[gene_wid]
        if isinstance(gene, Entry):
            chr_wid = gene.chromosome.wid
        else:
            chr_wid = gene['chromosome']

        #chr length
        if all_obj_data is None:
            chro = Genome.objects.get(species__wid=self.species, wid=chr_wid)
        else:
            chro = all_obj_data[chr_wid]
        if isinstance(chro, Entry):
            chr_len = chro.length
        else:
            chr_len = chro['length']

        #error check binding site coordinate, length
        if self.binding_site is not None:
            if self.binding_site.coordinate > chr_len:
                raise ValidationError({'binding_site': 'Coordinate must be less then chromosome length.'})
            if self.binding_site.length > chr_len:
                raise ValidationError({'binding_site': 'Length must be less then chromosome length.'})

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
            'publication_references'
            ]
        facet_fields = ['type', 'transcription_unit', 'transcription_factor']
        verbose_name='Transcriptional regulation'
        verbose_name_plural = 'Transcriptional regulation'
        wid_unique = False

class Type(SpeciesComponent):
    #parent pointer
    parent_ptr_species_component = OneToOneField(SpeciesComponent, related_name='child_ptr_type', parent_link=True, verbose_name='Species component')

    @classmethod
    def get_statistics(cls, species):
        # API:
        # [[0, "<a href="Chromosome">, 123],

        amount = cls.objects.for_species(species).count()
        url = cls.get_model_url(species)
        types = cls.objects.for_species(species).values("wid").annotate(count=Count("children")).filter(count__gt=0)

        if amount == 0:
            return None

        idx = 0

        lst = [[idx, cls._meta.verbose_name_plural, amount, None, url]]

        for typ in types:
            lst.append([idx+1, typ["wid"], typ["count"], None, "{}?type={}".format(url, typ["wid"])])

        return lst

    #getters
    def get_all_members(self):
        members = []
        for m in self.members.all():
            members.append(m)

        children_pks = self.children.values_list("pk", flat=True)
        t = Type.objects.filter(pk__in=children_pks)
        for c in t:
            members += c.get_all_members()
        return members

    #html formatting    
    def get_as_html_parent(self, is_user_anonymous):
        if self.parent is not None:
            result = '<a href="%s">%s</a>' % (self.parent.get_absolute_url(), self.parent.wid, )
            if self.parent.parent is not None:
                result = self.parent.get_as_html_parent(is_user_anonymous) + ' &#8250; ' + result
            return result

    def get_as_html_children(self, is_user_anonymous):
        from cyano.helpers import format_list_html_url
        return format_list_html_url(self.children.all())

    def get_as_html_members(self, is_user_anonymous):
        from cyano.helpers import format_list_html_url
        return format_list_html_url(self.get_all_members())

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
            'publication_references'
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
        unique_together = ('xid', 'source')
        verbose_name = 'Cross reference'
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
    def get_citation(self, cross_references = False):
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

        cr = self.get_as_html_cross_references(True)
        cr_spacer = ''
        if cr != '':
            cr_spacer = ', '
        return '%s CyanoFactory: <a href="%s">%s</a>%s%s' % (txt, self.get_absolute_url(), self.wid, cr_spacer, cr)

    def get_all_referenced_entries(self):
        entries = []
        for entry in self.publication_referenced_components.filter(species = self.species):
            entries.append(entry)
        for ev in Evidence.objects.filter(references__id=self.id):
            entries.append(ev.species_component)
        return entries

    #html formatting    
    def get_as_html_citation(self, is_user_anonymous):
        return self.get_citation()

    def get_as_html_referenced_entries(self, is_user_anonymous):
        from cyano.helpers import format_list_html_url
        return format_list_html_url(self.get_all_referenced_entries())

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
            'publication_references'
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
    def get_as_html_target_peptide(self, is_user_anonymous):
        results = []
        try:
            target_type = Type.objects.for_wid("Target-Peptide")
        except ObjectDoesNotExist:
            return ""

        from cyano.helpers import format_list_html_url
        return format_list_html_url(self.children.filter(type=target_type).order_by("wid"))

    def get_as_html_decoy_peptide(self, is_user_anonymous):
        results = []
        try:
            decoy_type = Type.objects.for_wid("Decoy-Peptide")
        except ObjectDoesNotExist:
            return ""

        from cyano.helpers import format_list_html_url
        return format_list_html_url(self.children.filter(type=decoy_type).order_by("wid"))

    def get_as_html_related_proteins(self, is_user_anonymous):
        from cyano.helpers import format_list_html_url
        return format_list_html_url(self.children.filter(model_type=TableMeta.get_by_model_name("MassSpectrometryProtein")).order_by("wid"))

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
            'id', 'wid', 'name', 'synonyms', 'cross_references', 'type',  'comments', 'publication_references'
        ]
        facet_fields = ['type']
        verbose_name = 'Mass Spectrometry Job'
        verbose_name_plural = 'Mass Spectrometry Jobs'
        wid_unique = False

class Peptide(Protein):
    #parent pointer
    parent_ptr_protein = OneToOneField(Protein, related_name='child_ptr_peptide', parent_link=True, verbose_name='Species component')

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
    def get_as_html_sequence(self, is_user_anonymous):
        from cyano.helpers import format_sequence_as_html
        return format_sequence_as_html(self.sequence)

    def get_as_html_matched_proteins(self, is_user_anonymous):
        results = []

        for matched_protein in self.proteins.all():
            try:
                res = Protein.objects.for_species(self.species).get(
                    parent=MassSpectrometryJob.objects.all()[0], wid__startswith=matched_protein.value + "_")
                results.append('<a href="%s">%s</a>' % (res.get_absolute_url(), res.wid))
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
            'id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'prosthetic_groups', 'chaperones', 'dna_footprint', 'regulatory_rule', 'comments', 'publication_references'
            ]
        facet_fields = ['type', 'chaperones', 'dna_footprint__binding', 'dna_footprint__region']
        verbose_name = 'Peptide'
        verbose_name_plural = 'Peptides'
        wid_unique = False


class MassSpectrometryProtein(Protein):
    parent_ptr_protein = OneToOneField(Protein, related_name='child_ptr_ms_protein', parent_link=True, verbose_name='Protein')

    score = FloatField(verbose_name="Protein Score")
    coverage = FloatField(verbose_name="% Coverage")
    sequence = TextField(verbose_name='Sequence', validators=[validate_protein_sequence])
    length = PositiveIntegerField(verbose_name='Length (nt)')
    ambiguous = ManyToManyField(EntryBasicTextData, verbose_name='Ambiguous Proteins', related_name='ambiguous')
    sub = ManyToManyField(EntryBasicTextData, verbose_name='Sub-Proteins', related_name='sub')
    pi = FloatField(verbose_name="Protein PI")
    mass = FloatField(verbose_name="Protein Mass (Da)")

    #html formatting
    def get_as_html_sequence(self, is_user_anonymous):
        from cyano.helpers import format_sequence_as_html
        return format_sequence_as_html(self.species, self.sequence)

    def get_as_html_structure(self, is_user_anonymous):
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
            fielddata = create_detail_fieldset(self.species, feature, fieldsets, False)

            tip_title = fielddata[0][0]
            tip_text = StringIO()

            for dataset in fielddata:
                for item in dataset[1]["fields"]:
                    tip_text.write("<br><b>" + item["verbose_name"] + "</b>: "+str(item["data"]))

            template = loader.get_template("cyano/genome/draw_feature.html")

            context_dict = {'h': featureHeight,
                            'title': tip_title,
                            'text': tip_text.getvalue(),
                            'url': feature.protein.get_absolute_url(),
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
        features = StringIO()

        for feature in self.protein_details.all():
            features.write(draw_segment(feature))

        H = 2 + geneHeight + 2 + 4 + 1 * (2 + len(feature_draw) * (featureHeight + 2)) + 2

        return '<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="%s" height="%s" viewport="0 0 %s %s"><style>%s%s%s</style><g class="chr">%s</g><g class="features">%s</g></svg>' % (
            W, H, W, H, style, chrStyle, featureStyle, chro, features.getvalue())

    def get_as_html_ambiguous_proteins(self, is_user_anonymous):
        results = []

        for ambiguous_protein in self.ambiguous.all():
            try:
                res = Protein.objects.for_species(self.species).get(
                    parent=MassSpectrometryJob.objects.all()[0], wid__startswith=ambiguous_protein.value + "_")
                results.append('<a href="%s">%s</a>' % (res.get_absolute_url(), res.wid))
            except ObjectDoesNotExist:
                results.append(ambiguous_protein.value)

        return format_list_html(results, comma_separated=True)

    def get_as_html_sub_proteins(self,  is_user_anonymous):
        results = []

        for sub_protein in self.sub.all():
            try:
                res = Protein.objects.for_species(self.species).get(
                    parent=MassSpectrometryJob.objects.all()[0], wid__startswith=sub_protein.value + "_")
                results.append('<a href="%s">%s</a>' % (res.get_absolute_url(), res.wid))
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
            'id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'prosthetic_groups', 'chaperones', 'dna_footprint', 'regulatory_rule', 'score', 'coverage', 'pi', 'mass', 'sequence', 'comments', 'publication_references'
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

class Basket(Model):
    user = ForeignKey(UserProfile, related_name = 'baskets', verbose_name = "Users baskets")
    name = CharField(max_length=255, blank=False, default='', verbose_name = "Basket name")
    
    def __str__(self):
        return self.name + " " + str(self.components.all())

class BasketComponent(Model):
    basket = ForeignKey(Basket, related_name = "components", verbose_name = "In basket")
    component = ForeignKey(SpeciesComponent, related_name = "+", verbose_name = "component")
    species = ForeignKey(Species, related_name = "+", verbose_name = "Species component belongs to")

    def __str__(self):
        return str(self.component)

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
        return separator.join(val[:default_items]) + separator\
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
        evidence = None
        for tmp in obj:
            if evidence is None:
                evidence = tmp.evidence.all()
            else:
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
            return '[<a href="%s">%s</a>]' % (obj.get_absolute_url(), obj.wid)
        except:
            return match.group(0)
    return inner_method

class ProteinComparison(models.Model):
    wid = models.AutoField(primary_key=True)
    protein_a = models.IntegerField()
    protein_b = models.IntegerField()
    equal_score = models.FloatField()
