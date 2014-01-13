"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""
from collections import OrderedDict

from django.core.exceptions import ValidationError

import cyano.models as cmodels
from cyano.helpers import slugify
from bioparser import BioParser
from django.db.transaction import commit_on_success


class _NoOverwriteDict(dict):
    # overwriting still possible with update
    def __setitem__(self, key, value):
        if key in self:
            raise ValueError("Key {} already in use".format(key))
        super(_NoOverwriteDict, self).__setitem__(key, value)


class EasyProt(BioParser):
    def __init__(self, wid, user, reason):
        super(EasyProt, self).__init__(wid, user, reason)

        self.export_parameters = _NoOverwriteDict()
        self.job_params = _NoOverwriteDict()
        self.target_peptides = []
        self.decoy_peptides = []
        self.protein_details = []
        self.protein_summary = []
        self.quant_stats_info = OrderedDict()
        self.quant_stats_info_complex = OrderedDict()
        self.mascat_quant_details = []
        self.libra_quant_details = []
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

        self._parse_export_parameters(sheet_export_parameters)
        self._parse_jobs_params(sheet_jobs_params)

        sheet_target_peptides = workbook.get_sheet_by_name("Target Peptides")
        if sheet_target_peptides is not None:
            self._parse_target_peptides(sheet_target_peptides)

        sheet_decoy_peptides = workbook.get_sheet_by_name("Decoy Peptides")
        if sheet_decoy_peptides is not None:
            self._parse_decoy_peptides(sheet_decoy_peptides)

        sheet_protein_details = workbook.get_sheet_by_name("Protein Details")
        if sheet_protein_details is not None:
            self._parse_protein_details(sheet_protein_details)

        sheet_protein_summary = workbook.get_sheet_by_name("Protein Summary")
        if sheet_protein_summary is not None:
            self._parse_protein_summary(sheet_protein_summary)

        sheet_quant_stats_info = workbook.get_sheet_by_name("Quant Stats Info")
        if sheet_quant_stats_info is not None:
            self._parse_quant_stats_info(sheet_quant_stats_info)

        sheet_mascat_quant_details = workbook.get_sheet_by_name("Mascat Quant Details")
        if sheet_mascat_quant_details is not None:
            self._parse_mascat_quant_details(sheet_mascat_quant_details)

        sheet_libra_quant_details = workbook.get_sheet_by_name("Libra Quant Details")
        if sheet_libra_quant_details is not None:
            self._parse_libra_quant_details(sheet_libra_quant_details)

        sheet_protein_expression = workbook.get_sheet_by_name("Protein Expression 95%")
        if sheet_protein_details is not None:
            self._parse_protein_expression(sheet_protein_expression)

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
        parse_isotopic = False
        for i in range(row):
            typ = sheet_export_parameters.cell(row=i, column=0).value
            value = sheet_export_parameters.cell(row=i, column=1).value

            if typ is not None:
                typ = typ.lower()

            if parse_jobs:
                if typ.startswith("job"):
                    value2 = sheet_export_parameters.cell(row=i, column=2).value
                    self.export_parameters["jobs"].append([value, value2])
                    continue
                else:
                    parse_jobs = False
            elif parse_isotopic:
                if typ:
                    parse_isotopic = False
                else:
                    self.export_parameters["isotopic purity correction"].append(value)
                    continue

            if typ is None:
                continue

            if typ == "easyprot version":
                typ = "version"
            elif typ == "#jobs":
                parse_jobs = True
                self.export_parameters["jobs"] = []
                continue
            elif typ.startswith("job"):
                raise ValidationError("Job outside of #Jobs list")
            elif typ == "isotopic purity correction":
                parse_isotopic = True
                self.export_parameters["isotopic purity correction"] = []
                self.export_parameters["isotopic purity correction"].append(value)
                continue

            self.export_parameters[typ] = value

    def _parse_jobs_params(self, sheet_jobs_params):
        """
        :type sheet_jobs_params: openpyxl.worksheet.Worksheet
        """
        row = sheet_jobs_params.get_highest_row()
        cols = sheet_jobs_params.get_highest_column()
        parse_rounds = False
        parse_mods = False
        self.job_params["rounds"] = []
        for i in range(row):
            typ = sheet_jobs_params.cell(row=i, column=0).value
            value = sheet_jobs_params.cell(row=i, column=1).value
            #print i, typ, value, self.job_params

            if typ is not None:
                typ = typ.lower()

            if parse_rounds:
                first_mod = False
                if typ == "modifications":
                    parse_mods = True
                    first_mod = True
                    for x in self.job_params["rounds"]:
                        x[typ] = []

                for j, x in enumerate(self.job_params["rounds"]):
                    value = sheet_jobs_params.cell(row=i, column=1).value
                    if value is not None:
                        parse_rounds = False
                        break
                    value = sheet_jobs_params.cell(row=i, column=j+2).value
                    if parse_mods:
                        new_typ = sheet_jobs_params.cell(row=i, column=0).value
                        if not first_mod and new_typ is not None:
                            parse_mods = False
                            typ = new_typ.lower()
                        else:
                            x["modifications"].append(value)
                    if not parse_mods:
                        x[typ.lower()] = value

            if typ is None:
                # Could be start of rounds, check 3rd cell
                typ = sheet_jobs_params.cell(row=i, column=2).value
                if typ is not None and typ.startswith("Round"):
                    parse_rounds = True
                    for x in range(cols):
                        typ = sheet_jobs_params.cell(row=i, column=x+2).value
                        if typ is None or not typ.startswith("Round"):
                            break
                        self.job_params["rounds"].append(_NoOverwriteDict())
                continue

            if typ == "job id":
                if not any(value == job[1] for job in self.export_parameters["jobs"]):
                    raise ValidationError("Unknown Job ID {}".format(value))

            self.job_params[typ] = value

    def _parse_target_peptides(self, sheet_target_peptides):
        self._read_table(sheet_target_peptides, self.target_peptides)

    def _parse_decoy_peptides(self, sheet_decoy_peptides):
        self._read_table(sheet_decoy_peptides, self.decoy_peptides)

    def _parse_protein_details(self, sheet_protein_details):
        self._read_table(sheet_protein_details, self.protein_details)

    def _parse_protein_summary(self, sheet_protein_summary):
        self._read_table(sheet_protein_summary, self.protein_summary)

    def _parse_quant_stats_info(self, sheet_quant_stats_info):
        """
        :type sheet_quant_stats_info: openpyxl.worksheet.Worksheet
        """
        # Parse downwards until cell(1) contains text
        rows = sheet_quant_stats_info.get_highest_row()
        cols = sheet_quant_stats_info.get_highest_column()
        START = 0
        PARSE_FIRST = 1
        BETWEEN = 2
        PARSE_LAST_HEADER = 3
        PARSE_LAST_DATA = 4
        section_type = START

        for row in range(rows):
            #print self.quant_stats_info
            print self.quant_stats_info_complex
            if section_type == START:
                if sheet_quant_stats_info.cell(row=row, column=1).value is not None:
                    section_type = PARSE_FIRST
                    for col in range(cols):
                        value = sheet_quant_stats_info.cell(row=row, column=col+1).value
                        if value is not None:
                            self.quant_stats_info[value] = OrderedDict()
            elif section_type == PARSE_FIRST:
                if sheet_quant_stats_info.cell(row=row, column=0).value is None:
                    section_type = BETWEEN
                else:
                    typ = sheet_quant_stats_info.cell(row=row, column=0).value
                    for i, item in enumerate(self.quant_stats_info, start=1):
                        print item
                        self.quant_stats_info[item][typ] = sheet_quant_stats_info.cell(row=row, column=i).value
            elif section_type == BETWEEN:
                if sheet_quant_stats_info.cell(row=row, column=1).value is not None:
                    section_type = PARSE_LAST_HEADER
                    for col in range(cols)[::2]:
                        value = sheet_quant_stats_info.cell(row=row, column=col+1).value
                        if value is not None:
                            self.quant_stats_info_complex[value] = OrderedDict()
            elif section_type == PARSE_LAST_HEADER:
                section_type = PARSE_LAST_DATA
                for data, col in zip(self.quant_stats_info_complex, range(cols)[::2]):
                    self.quant_stats_info_complex[data] = OrderedDict()
                    self.quant_stats_info_complex[data][sheet_quant_stats_info.cell(row=row, column=col+1).value] = OrderedDict()
                    self.quant_stats_info_complex[data][sheet_quant_stats_info.cell(row=row, column=col+2).value] = OrderedDict()
            elif section_type == PARSE_LAST_DATA:
                if sheet_quant_stats_info.cell(row=row, column=0).value is None:
                    break
                else:
                    typ = sheet_quant_stats_info.cell(row=row, column=0).value
                    for i, item in zip(range(len(self.quant_stats_info_complex)*2)[::2], self.quant_stats_info_complex):
                        val1 = sheet_quant_stats_info.cell(row=row, column=i+1).value
                        val2 = sheet_quant_stats_info.cell(row=row, column=i+2).value
                        it = iter(self.quant_stats_info_complex[item])
                        self.quant_stats_info_complex[item][it.next()][typ] = val1
                        self.quant_stats_info_complex[item][it.next()][typ] = val2

        #self._read_table_with_annotation(sheet_quant_stats_info, self.quant_stats_info_header, self.quant_stats_info)

        # Strategy:
        # Parse downwards until cell(1) contains text
        #  - Create list with size of filled cells
        # Parse downwards and add to list until cell(0) is empty
        # Parse downwards until cell(1) contains text
        #  - Create list with size of filled cells
        #  - !!! Cells are connected, size 2 hardcoded
        #  - Create sublist (2 elements per connected cell)
        # Parse downwards and add to list until cell(0) is empty
        # RETURN

    def _parse_mascat_quant_details(self, sheet_mascat_quant_details):
        self._read_table_with_annotation(sheet_mascat_quant_details, self.mascat_quant_details)

    def _parse_libra_quant_details(self, sheet_libra_quant_details):
        self._read_table_with_annotation(sheet_libra_quant_details, self.libra_quant_details)

    def _parse_protein_expression(self, sheet_protein_expression):
        # Strategy:
        # Parse to the right until filled cell hit.
        #  - Remember pos + 1
        # Parse more to the right
        #  - Remember pos + 1
        # Read header in next lines
        # Normal table scan
        # When first cell empty and cell remembered from beg. filled:
        #  - Goto strategy from beginning
        # Do until col max reached -> RETURN
        pass

    def _read_table(self, sheet, data):
        """
        :type sheet: openpyxl.worksheet.Worksheet
        :type data: list
        """
        cols = sheet.get_highest_column()
        header = []

        for i in range(cols):
            val = sheet.cell(row=0, column=i).value
            if val is None:
                break
            header.append(sheet.cell(row=0, column=i).value)

        # header omitted
        rows = sheet.get_highest_row() - 1

        for i in range(rows):
            odict = OrderedDict()
            for j in range(len(header)):
                odict[header[j]] = sheet.cell(row=i+1, column=j).value
            data.append(odict)

    def _read_table_with_annotation(self, sheet, data):
        """
        :type sheet: openpyxl.worksheet.Worksheet
        :type data: list
        """
        cols = sheet.get_highest_column()
        header = []

        current_group = "main"
        header_group = dict()
        header_group[0] = current_group

        # Read header and group elements
        for i in range(cols):
            current_group = sheet.cell(row=0, column=i).value
            if current_group is not None:
                header_group[i+1] = current_group

            head = sheet.cell(row=1, column=i).value
            header.append(head)

        # header omitted, read values
        rows = sheet.get_highest_row() - 2
        for i in range(rows):
            odict = OrderedDict()

            for j in range(len(header)):
                if j in header_group:
                    current_group = header_group[j]
                    odict[current_group] = OrderedDict()
                odict[current_group][header[j]] = sheet.cell(row=i+2, column=j).value
            data.append(odict)

    @commit_on_success
    def apply(self):
        self.detail.save()

        self.species.save(self.detail)

        ms_job = cmodels.MassSpectrometryJob.for_species(self.species).for_wid(self.export_parameters["jobs"][0][0], create=True)
        ms_job.name = ms_job.wid
        ms_job.save(self.detail)
        ms_job.species.add(self.species)

        # Create types 'Target-Peptide' and 'Decoy-Peptide' if missing
        target_type = cmodels.Type.objects.for_wid("Target-Peptide", create=True)
        decoy_type = cmodels.Type.objects.for_wid("Decoy-Peptide", create=True)

        for item in self.target_peptides:
            peptide = cmodels.Peptide.for_species(self.species).for_wid("???", create=True)

            peptide.parent = ms_job
            peptide.sequence = item["Sequence"]
            peptide.proteotypic = item["Proteotypic"]
            peptide.charge = item["Charge"]
            peptide.mass = item["m/z"]
            peptide.zscore = item["zscore"]
            peptide.retention_time = item["RT"]

        if hasattr(self, "notify_progress"):
            outstr = "Importing EasyProt file"
            self.notify_progress(current=len(cds_map.values()), total=len(cds_map.values()), message=outstr)
