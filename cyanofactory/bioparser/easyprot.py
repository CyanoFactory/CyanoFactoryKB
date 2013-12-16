"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from Bio import SeqIO
from django.core.exceptions import ValidationError
from bioparser.easyprot import NoOverwriteDict

import cyano.models as cmodels
from cyano.helpers import slugify
from bioparser import BioParser
from django.db.transaction import commit_on_success


class EasyProt(BioParser):
    class NoOverwriteDict(dict):
        # overwriting still possible with
        def __setitem__(self, key, value):
            if key in self:
                raise ValueError("Key {} already in use".format(key))
            super(NoOverwriteDict, self).__setitem__(key, value)

    def __init__(self, wid, user, reason):
        super(EasyProt, self).__init__(wid, user, reason)

        self.export_parameters = NoOverwriteDict()
        self.job_params = NoOverwriteDict()
        self.target_peptides = []
        self.decoy_peptides = []
        # PTM: Post Translational Modifications

    def parse(self, handle):
        if hasattr(self, "notify_progress"):
            self.notify_progress(current=0, total=1, message="Parsing EasyProt file...")

        from openpyxl.reader.excel import load_workbook
        workbook = load_workbook(handle)

        sheet_export_parameters = workbook.get_sheet_by_name("Export Parameters")
        if sheet_export_parameters is None:
            raise ValidationError("Sheet Export Parameters missing")

        sheet_jobs_params = workbook.get_sheet_by_name("Jobs Params")
        if sheet_jobs_params is None:
            raise ValidationError("Sheet Jobs Params missing")

        sheet_target_peptides = workbook.get_sheet_by_name("Target Peptides")
        if sheet_target_peptides is None:
            raise ValidationError("Sheet Target Peptides missing")

        sheet_decoy_peptides = workbook.get_sheet_by_name("Decoy Peptides")
        if sheet_decoy_peptides is None:
            raise ValidationError("Sheet Decoy Peptides missing")

        self._parse_export_parameters(sheet_export_parameters)
        self._parse_jobs_params(sheet_jobs_params)
        self._parse_target_peptides(sheet_target_peptides)
        self._parse_decoy_peptides(sheet_decoy_peptides)

    @staticmethod
    def _add_to_dict(dictionary, key, value):
        if key in dictionary:
            raise ValidationError("Key {} already in use".format(key))
        dictionary[key] = value

    def _parse_export_parameters(self, sheet_export_parameters):
        """
        :type sheet_export_parameters: openpyxl.worksheet.Worksheet
        """
        row = sheet_export_parameters.get_highest_row()
        parse_jobs = False
        for i in range(row):
            typ = sheet_export_parameters.cell(row=i, column=0).value
            value = sheet_export_parameters.cell(row=i, column=1).value

            typ = typ.lower()

            if parse_jobs:
                if typ.startswith("job"):
                    value2 = sheet_export_parameters.cell(row=i, column=2).value
                    self.export_parameters["jobs"].append([value, value2])
                    continue
                else:
                    parse_jobs = False

            if typ == "easyprot version":
                typ = "version"
            elif typ == "#jobs":
                parse_jobs = True
                self.export_parameters["jobs"] = []
                continue
            elif typ.startswith("job"):
                raise ValidationError("Job outside of #Jobs list")

            self.export_parameters[typ] = value

    def _parse_jobs_params(self, sheet_jobs_params):
        pass

    def _parse_target_peptides(self, sheet_target_peptides):
        """
        :type sheet_target_peptides: openpyxl.worksheet.Worksheet
        """
        cols = sheet_target_peptides.get_highest_column()
        headers = []
        values = map(lambda x: [], range(cols))
        for i in range(cols):
            headers.append(sheet_target_peptides.cell(row=0, column=i).value)

        # header omitted
        rows = sheet_target_peptides.get_highest_row() - 1
        for i in range(rows):
            for j in range(cols):
                values[i].append(sheet_target_peptides.cell(row=i+1, column=j).value)

    def _parse_decoy_peptides(self, sheet_decoy_peptides):
        """
        :type sheet_decoy_peptides: openpyxl.worksheet.Worksheet
        """
        cols = sheet_decoy_peptides.get_highest_column()
        headers = []
        values = map(lambda x: [], range(cols))
        for i in range(cols):
            headers.append(sheet_decoy_peptides.cell(row=0, column=i).value)

        # header omitted
        rows = sheet_decoy_peptides.get_highest_row() - 1
        for i in range(rows):
            for j in range(cols):
                values[i].append(sheet_decoy_peptides.cell(row=i+1, column=j).value)

    @commit_on_success
    def apply(self):
        self.detail.save()
        
        self.species.save(self.detail)
            
        if hasattr(self, "notify_progress"):
            outstr = "Assigning KEGG pathways"
            self.notify_progress(current = len(cds_map.values()), total = len(cds_map.values()), message = outstr)

