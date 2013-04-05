# This is an auto-generated Django model module.
# You'll have to do the following manually to clean this up:
#     * Rearrange models' order
#     * Make sure each model has one field with primary_key=True
# Feel free to rename the models, but don't rename db_table values or field names.
#
# Also note: You'll have to insert the output of 'django-admin.py sqlcustom [appname]'
# into your database.
from __future__ import unicode_literals

from django.db import models

class Archive(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    otherWid = models.BigIntegerField(db_column='OtherWID') # Field name made lowercase.
    format = models.CharField(max_length=10, db_column='Format') # Field name made lowercase.
    contents = models.TextField(db_column='Contents', blank=True) # Field name made lowercase. This field type is a guess.
    url = models.TextField(db_column='URL', blank=True) # Field name made lowercase.
    toolname = models.CharField(max_length=50, db_column='ToolName', blank=True) # Field name made lowercase.
    datasetWid = models.BigIntegerField(db_column='DataSetWID') # Field name made lowercase.
    class Meta:
        db_table = 'Archive'

class Arraydesign(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey('Dataset', db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    mageclass = models.CharField(max_length=100, db_column='MAGEClass') # Field name made lowercase.
    identifier = models.CharField(max_length=255, db_column='Identifier', blank=True) # Field name made lowercase.
    name = models.CharField(max_length=255, db_column='Name', blank=True) # Field name made lowercase.
    version = models.CharField(max_length=255, db_column='Version', blank=True) # Field name made lowercase.
    numberoffeatures = models.SmallIntegerField(null=True, db_column='NumberOfFeatures', blank=True) # Field name made lowercase.
    surfacetype = models.ForeignKey('Term', null=True, db_column='SurfaceType', blank=True, related_name = '+') # Field name made lowercase.
    compositegrp = models.ManyToManyField('Designelementgroup', related_name = 'arraydesign_compositegrp', through = 'ArraydesignWidcompositegrpWid')
    contact = models.ManyToManyField('Contact', related_name = 'arraydesign_contact', through = 'ArraydesignWidcontactWid')
    reportergroup = models.ManyToManyField('Designelementgroup', related_name = 'arraydesign_reportergroup', through = 'ArraydesignWidreportergroupWid')
    class Meta:
        db_table = 'ArrayDesign'

class ArraydesignWidcompositegrpWid(models.Model):
    arraydesignWid = models.ForeignKey(Arraydesign, db_column='ArrayDesignWID', related_name = '+') # Field name made lowercase.
    compositegroupWid = models.ForeignKey('Designelementgroup', db_column='CompositeGroupWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'ArrayDesignWIDCompositeGrpWID'

class ArraydesignWidcontactWid(models.Model):
    arraydesignWid = models.ForeignKey(Arraydesign, db_column='ArrayDesignWID', related_name = '+') # Field name made lowercase.
    contactWid = models.ForeignKey('Contact', db_column='ContactWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'ArrayDesignWIDContactWID'

class ArraydesignWidreportergroupWid(models.Model):
    arraydesignWid = models.ForeignKey(Arraydesign, db_column='ArrayDesignWID', related_name = '+') # Field name made lowercase.
    reportergroupWid = models.ForeignKey('Designelementgroup', db_column='ReporterGroupWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'ArrayDesignWIDReporterGroupWID'

class Arraygroup(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey('Dataset', db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    identifier = models.CharField(max_length=255, db_column='Identifier', blank=True) # Field name made lowercase.
    name = models.CharField(max_length=255, db_column='Name', blank=True) # Field name made lowercase.
    barcode = models.CharField(max_length=255, db_column='Barcode', blank=True) # Field name made lowercase.
    arrayspacingx = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='ArraySpacingX', blank=True) # Field name made lowercase.
    arrayspacingy = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='ArraySpacingY', blank=True) # Field name made lowercase.
    numarrays = models.SmallIntegerField(null=True, db_column='NumArrays', blank=True) # Field name made lowercase.
    orientationmark = models.CharField(max_length=255, db_column='OrientationMark', blank=True) # Field name made lowercase.
    orientationmarkposition = models.CharField(max_length=25, db_column='OrientationMarkPosition', blank=True) # Field name made lowercase.
    Width = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='Width', blank=True) # Field name made lowercase.
    length = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='Length', blank=True) # Field name made lowercase.
    arraygroup_substratetype = models.ForeignKey('Term', null=True, db_column='ArrayGroup_SubstrateType', blank=True, related_name = '+') # Field name made lowercase.
    arraygroup_distanceunit = models.ForeignKey('Unit', null=True, db_column='ArrayGroup_DistanceUnit', blank=True, related_name = '+') # Field name made lowercase.
    array = models.ManyToManyField('Array', related_name = 'arraygroup_array', through = 'ArraygroupWidarrayWid')

    class Meta:
        db_table = 'ArrayGroup'

class ArraygroupWidarrayWid(models.Model):
    arraygroupWid = models.ForeignKey(Arraygroup, db_column='ArrayGroupWID', related_name = '+') # Field name made lowercase.
    arrayWid = models.ForeignKey('Array', db_column='ArrayWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'ArrayGroupWIDArrayWID'

class Arraymanufacture(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey('Dataset', db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    identifier = models.CharField(max_length=255, db_column='Identifier', blank=True) # Field name made lowercase.
    name = models.CharField(max_length=255, db_column='Name', blank=True) # Field name made lowercase.
    manufacturingdate = models.CharField(max_length=255, db_column='ManufacturingDate', blank=True) # Field name made lowercase.
    tolerance = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='Tolerance', blank=True) # Field name made lowercase.
    array = models.ManyToManyField('Array', related_name = 'arraymanufacture_array', through = 'ArraymanufactureWidarrayWid')
    contact = models.ManyToManyField('Contact', related_name = 'arraymanufacture_contact', through = 'ArraymanufactureWidcontactWid')
    class Meta:
        db_table = 'ArrayManufacture'

class Arraymanufacturedeviation(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey('Dataset', db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    array_field = models.ForeignKey('Array', null=True, db_column='Array_', blank=True, related_name = '+') # Field name made lowercase. Field renamed because it ended with '_'.
    class Meta:
        db_table = 'ArrayManufactureDeviation'

class ArraymanufactureWidarrayWid(models.Model):
    arraymanufactureWid = models.ForeignKey(Arraymanufacture, db_column='ArrayManufactureWID', related_name = '+') # Field name made lowercase.
    arrayWid = models.ForeignKey('Array', db_column='ArrayWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'ArrayManufactureWIDArrayWID'

class ArraymanufactureWidcontactWid(models.Model):
    arraymanufactureWid = models.ForeignKey(Arraymanufacture, db_column='ArrayManufactureWID', related_name = '+') # Field name made lowercase.
    contactWid = models.ForeignKey('Contact', db_column='ContactWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'ArrayManufactureWIDContactWID'

class Array(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey('Dataset', db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    identifier = models.CharField(max_length=255, db_column='Identifier', blank=True) # Field name made lowercase.
    name = models.CharField(max_length=255, db_column='Name', blank=True) # Field name made lowercase.
    arrayidentifier = models.CharField(max_length=255, db_column='ArrayIdentifier', blank=True) # Field name made lowercase.
    arrayxorigin = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='ArrayXOrigin', blank=True) # Field name made lowercase.
    arrayyorigin = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='ArrayYOrigin', blank=True) # Field name made lowercase.
    originrelativeto = models.CharField(max_length=255, db_column='OriginRelativeTo', blank=True) # Field name made lowercase.
    arraydesign = models.ForeignKey(Arraydesign, null=True, db_column='ArrayDesign', blank=True, related_name = '+') # Field name made lowercase.
    information = models.ForeignKey(Arraymanufacture, null=True, db_column='Information', blank=True, related_name = '+') # Field name made lowercase.
    arraygroup = models.ForeignKey(Arraygroup, null=True, db_column='ArrayGroup', blank=True, related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'Array_'

class BassaymappingWidbassaymapWid(models.Model):
    bioassaymappingWid = models.ForeignKey('Bioassaymapping', db_column='BioAssayMappingWID', related_name = '+') # Field name made lowercase.
    bioassaymapWid = models.ForeignKey('Bioevent', db_column='BioAssayMapWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'BAssayMappingWIDBAssayMapWID'

class Bioassay(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey('Dataset', db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    mageclass = models.CharField(max_length=100, db_column='MAGEClass') # Field name made lowercase.
    identifier = models.CharField(max_length=255, db_column='Identifier', blank=True) # Field name made lowercase.
    name = models.CharField(max_length=255, db_column='Name', blank=True) # Field name made lowercase.
    derivedbioassay_type = models.ForeignKey('Term', null=True, db_column='DerivedBioAssay_Type', blank=True, related_name = '+') # Field name made lowercase.
    featureextraction = models.ForeignKey('Bioevent', null=True, db_column='FeatureExtraction', blank=True, related_name = '+') # Field name made lowercase.
    bioassaycreation = models.ForeignKey('Bioevent', null=True, db_column='BioAssayCreation', blank=True, related_name = '+') # Field name made lowercase.
    channel = models.ManyToManyField('Channel', related_name = 'bioassay_channel', through = 'BioassayWidchannelWid')
    factorvalue = models.ManyToManyField('Factorvalue', related_name = 'bioassay_factorvalue', through = 'BioassayWidfactorvalueWid')
    derivbioadata = models.ManyToManyField('Bioassaydata', related_name = 'derivbioa_derivbioadata', through = 'DerivbioaWidderivbioadataWid')
    bioassaymap = models.ManyToManyField('Bioevent', related_name = 'derivbioassay_bioassaymap', through = 'DerivbioassayWidbioassaymapWid')
    measbassaydata = models.ManyToManyField('Bioassaydata', related_name = 'measbassay_measbassaydata', through = 'MeasbassayWidmeasbassaydataWid')
    class Meta:
        db_table = 'BioAssay'

class Bioassaydata(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey('Dataset', db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    mageclass = models.CharField(max_length=100, db_column='MAGEClass') # Field name made lowercase.
    identifier = models.CharField(max_length=255, db_column='Identifier', blank=True) # Field name made lowercase.
    name = models.CharField(max_length=255, db_column='Name', blank=True) # Field name made lowercase.
    bioassaydimension = models.ForeignKey('Bioassaydimension', null=True, db_column='BioAssayDimension', blank=True, related_name = '+') # Field name made lowercase.
    designelementdimension = models.ForeignKey('Designelementdimension', null=True, db_column='DesignElementDimension', blank=True, related_name = '+') # Field name made lowercase.
    quantitationtypedimension = models.ForeignKey('Quantitationtypedimension', null=True, db_column='QuantitationTypeDimension', blank=True, related_name = '+') # Field name made lowercase.
    bioassaydata_biodatavalues = models.ForeignKey('Biodatavalues', null=True, db_column='BioAssayData_BioDataValues', blank=True, related_name = '+') # Field name made lowercase.
    producertransformation = models.ForeignKey('Bioevent', null=True, db_column='ProducerTransformation', blank=True, related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'BioAssayData'

class Bioassaydatacluster(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey('Dataset', db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    identifier = models.CharField(max_length=255, db_column='Identifier', blank=True) # Field name made lowercase.
    name = models.CharField(max_length=255, db_column='Name', blank=True) # Field name made lowercase.
    clusterbioassaydata = models.ForeignKey(Bioassaydata, null=True, db_column='ClusterBioAssayData', blank=True, related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'BioAssayDataCluster'

class BioassaydimensioWidbioassayWid(models.Model):
    bioassaydimensionWid = models.ForeignKey('Bioassaydimension', db_column='BioAssayDimensionWID', related_name = '+') # Field name made lowercase.
    bioassayWid = models.ForeignKey(Bioassay, db_column='BioAssayWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'BioAssayDimensioWIDBioAssayWID'

class Bioassaydimension(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey('Dataset', db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    identifier = models.CharField(max_length=255, db_column='Identifier', blank=True) # Field name made lowercase.
    name = models.CharField(max_length=255, db_column='Name', blank=True) # Field name made lowercase.
    bioassay = models.ManyToManyField('Bioassay', related_name = 'bioassaydimensio_bioassay', through = 'BioassaydimensioWidbioassayWid')
    class Meta:
        db_table = 'BioAssayDimension'

class BioassaymapWidbioassayWid(models.Model):
    bioassaymapWid = models.ForeignKey('Bioevent', db_column='BioAssayMapWID', related_name = '+') # Field name made lowercase.
    bioassayWid = models.ForeignKey(Bioassay, db_column='BioAssayWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'BioAssayMapWIDBioAssayWID'

class Bioassaymapping(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey('Dataset', db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    bassaymap = models.ManyToManyField('Bioevent', related_name = 'bassaymapping_bassaymap', through = 'BassaymappingWidbassaymapWid')
    class Meta:
        db_table = 'BioAssayMapping'

class Bioassaytuple(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey('Dataset', db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    bioassay = models.ForeignKey(Bioassay, null=True, db_column='BioAssay', blank=True, related_name = '+') # Field name made lowercase.
    biodatatuples_bioassaytuples = models.ForeignKey('Biodatavalues', null=True, db_column='BioDataTuples_BioAssayTuples', blank=True, related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'BioAssayTuple'

class BioassayWidchannelWid(models.Model):
    bioassayWid = models.ForeignKey(Bioassay, db_column='BioAssayWID', related_name = '+') # Field name made lowercase.
    channelWid = models.ForeignKey('Channel', db_column='ChannelWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'BioAssayWIDChannelWID'

class BioassayWidfactorvalueWid(models.Model):
    bioassayWid = models.ForeignKey(Bioassay, db_column='BioAssayWID', related_name = '+') # Field name made lowercase.
    factorvalueWid = models.ForeignKey('Factorvalue', db_column='FactorValueWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'BioAssayWIDFactorValueWID'

class Biodatavalues(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey('Dataset', db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    mageclass = models.CharField(max_length=100, db_column='MAGEClass') # Field name made lowercase.
    order_field = models.CharField(max_length=25, db_column='Order_', blank=True) # Field name made lowercase. Field renamed because it ended with '_'.
    biodatacube_datainternal = models.ForeignKey('Datainternal', null=True, db_column='BioDataCube_DataInternal', blank=True, related_name = '+') # Field name made lowercase.
    biodatacube_dataexternal = models.ForeignKey('Dataexternal', null=True, db_column='BioDataCube_DataExternal', blank=True, related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'BioDataValues'

class Bioevent(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey('Dataset', db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    mageclass = models.CharField(max_length=100, db_column='MAGEClass') # Field name made lowercase.
    identifier = models.CharField(max_length=255, db_column='Identifier', blank=True) # Field name made lowercase.
    name = models.CharField(max_length=255, db_column='Name', blank=True) # Field name made lowercase.
    compositesequence = models.ForeignKey('Designelement', null=True, db_column='CompositeSequence', blank=True, related_name = '+') # Field name made lowercase.
    reporter = models.ForeignKey('Designelement', null=True, db_column='Reporter', blank=True, related_name = '+') # Field name made lowercase.
    compositesequence2 = models.ForeignKey('Designelement', null=True, db_column='CompositeSequence2', blank=True, related_name = '+') # Field name made lowercase.
    bioassaymaptarget = models.ForeignKey(Bioassay, null=True, db_column='BioAssayMapTarget', blank=True, related_name = '+') # Field name made lowercase.
    targetquantitationtype = models.ForeignKey('Quantitationtype', null=True, db_column='TargetQuantitationType', blank=True, related_name = '+') # Field name made lowercase.
    derivedbioassaydatatarget = models.ForeignKey(Bioassaydata, null=True, db_column='DerivedBioAssayDataTarget', blank=True, related_name = '+') # Field name made lowercase.
    quantitationtypemapping = models.ForeignKey('Quantitationtypemapping', null=True, db_column='QuantitationTypeMapping', blank=True, related_name = '+') # Field name made lowercase.
    designelementmapping = models.ForeignKey('Designelementmapping', null=True, db_column='DesignElementMapping', blank=True, related_name = '+') # Field name made lowercase.
    transformation_bioassaymapping = models.ForeignKey(Bioassaymapping, null=True, db_column='Transformation_BioAssayMapping', blank=True, related_name = '+') # Field name made lowercase.
    biomaterial_treatments = models.ForeignKey('Biosource', null=True, db_column='BioMaterial_Treatments', blank=True, related_name = '+') # Field name made lowercase.
    order_field = models.SmallIntegerField(null=True, db_column='Order_', blank=True) # Field name made lowercase. Field renamed because it ended with '_'.
    treatment_action = models.ForeignKey('Term', null=True, db_column='Treatment_Action', blank=True, related_name = '+') # Field name made lowercase.
    treatment_actionmeasurement = models.ForeignKey('Measurement', null=True, db_column='Treatment_ActionMeasurement', blank=True, related_name = '+') # Field name made lowercase.
    array_field = models.ForeignKey(Array, null=True, db_column='Array_', blank=True, related_name = '+') # Field name made lowercase. Field renamed because it ended with '_'.
    physicalbioassaytarget = models.ForeignKey(Bioassay, null=True, db_column='PhysicalBioAssayTarget', blank=True, related_name = '+') # Field name made lowercase.
    physicalbioassay = models.ForeignKey(Bioassay, null=True, db_column='PhysicalBioAssay', blank=True, related_name = '+') # Field name made lowercase.
    target = models.ForeignKey(Bioassay, null=True, db_column='Target', blank=True, related_name = '+') # Field name made lowercase.
    physicalbioassaysource = models.ForeignKey(Bioassay, null=True, db_column='PhysicalBioAssaySource', blank=True, related_name = '+') # Field name made lowercase.
    measuredbioassaytarget = models.ForeignKey(Bioassay, null=True, db_column='MeasuredBioAssayTarget', blank=True, related_name = '+') # Field name made lowercase.
    physicalbioassay2 = models.ForeignKey(Bioassay, null=True, db_column='PhysicalBioAssay2', blank=True, related_name = '+') # Field name made lowercase.
    bioassay = models.ManyToManyField('Bioassay', related_name = 'bioassaymap_bioassay', through = 'BioassaymapWidbioassayWid')
    image = models.ManyToManyField('Image', related_name = 'imageacquisition_image', through = 'ImageacquisitionWidimageWid')
    quanttype = models.ManyToManyField('Quantitationtype', related_name = 'quanttypemap_quanttype', through = 'QuanttypemapWidquanttypeWid')
    bioassaydata = models.ManyToManyField('Bioassaydata', related_name = 'transform_bioassaydata', through = 'TransformWidbioassaydataWid')
    class Meta:
        db_table = 'BioEvent'

class Biomaterialmeasurement(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey('Dataset', db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    biomaterial = models.ForeignKey('Biosource', null=True, db_column='BioMaterial', blank=True, related_name = '+') # Field name made lowercase.
    measurement = models.ForeignKey('Measurement', null=True, db_column='Measurement', blank=True, related_name = '+') # Field name made lowercase.
    treatment = models.ForeignKey(Bioevent, null=True, db_column='Treatment', blank=True, related_name = '+') # Field name made lowercase.
    bioassaycreation = models.ForeignKey(Bioevent, null=True, db_column='BioAssayCreation', blank=True, related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'BioMaterialMeasurement'

class Biosource(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    mageclass = models.CharField(max_length=100, db_column='MAGEClass', blank=True) # Field name made lowercase.
    taxonWid = models.ForeignKey('Taxon', null=True, db_column='TaxonWID', blank=True, related_name = '+') # Field name made lowercase.
    name = models.CharField(max_length=200, db_column='Name', blank=True) # Field name made lowercase.
    strain = models.CharField(max_length=220, db_column='Strain', blank=True) # Field name made lowercase.
    organ = models.CharField(max_length=50, db_column='Organ', blank=True) # Field name made lowercase.
    organelle = models.CharField(max_length=50, db_column='Organelle', blank=True) # Field name made lowercase.
    tissue = models.CharField(max_length=100, db_column='Tissue', blank=True) # Field name made lowercase.
    celltype = models.CharField(max_length=50, db_column='CellType', blank=True) # Field name made lowercase.
    cellline = models.CharField(max_length=50, db_column='CellLine', blank=True) # Field name made lowercase.
    atccid = models.CharField(max_length=50, db_column='ATCCId', blank=True) # Field name made lowercase.
    diseased = models.CharField(max_length=1, db_column='Diseased', blank=True) # Field name made lowercase.
    disease = models.CharField(max_length=250, db_column='Disease', blank=True) # Field name made lowercase.
    developmentstage = models.CharField(max_length=50, db_column='DevelopmentStage', blank=True) # Field name made lowercase.
    sex = models.CharField(max_length=15, db_column='Sex', blank=True) # Field name made lowercase.
    datasetWid = models.ForeignKey('Dataset', db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    biosubtype = models.ManyToManyField('Biosubtype', related_name = 'biosource_biosubtype', through = 'BiosourceWidbiosubtypeWid')
    contact = models.ManyToManyField('Contact', related_name = 'biosource_contact', through = 'BiosourceWidcontactWid')
    gene = models.ManyToManyField('Gene', related_name = 'biosource_gene', through = 'BiosourceWidgeneWid')
    protein = models.ManyToManyField('Protein', related_name = 'biosource_protein', through = 'BiosourceWidproteinWid')
    compound = models.ManyToManyField('Chemical', related_name = 'labeledextract_compound', through = 'LabeledextractWidcompoundWid')

    class Meta:
        db_table = 'BioSource'

class BiosourceWidbiosubtypeWid(models.Model):
    biosourceWid = models.ForeignKey(Biosource, db_column='BioSourceWID', related_name = '+') # Field name made lowercase.
    biosubtypeWid = models.ForeignKey('Biosubtype', db_column='BioSubtypeWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'BioSourceWIDBioSubtypeWID'

class BiosourceWidcontactWid(models.Model):
    biosourceWid = models.ForeignKey(Biosource, db_column='BioSourceWID', related_name = '+') # Field name made lowercase.
    contactWid = models.ForeignKey('Contact', db_column='ContactWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'BioSourceWIDContactWID'

class BiosourceWidgeneWid(models.Model):
    biosourceWid = models.ForeignKey(Biosource, db_column='BioSourceWID', related_name = '+') # Field name made lowercase.
    geneWid = models.ForeignKey('Gene', db_column='GeneWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'BioSourceWIDGeneWID'

class BiosourceWidproteinWid(models.Model):
    biosourceWid = models.ForeignKey(Biosource, db_column='BioSourceWID', related_name = '+') # Field name made lowercase.
    proteinWid = models.ForeignKey('Protein', db_column='ProteinWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'BioSourceWIDProteinWID'

class Biosubtype(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    type = models.CharField(max_length=25, db_column='Type', blank=True) # Field name made lowercase.
    value = models.CharField(max_length=50, db_column='Value') # Field name made lowercase.
    datasetWid = models.ForeignKey('Dataset', db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'BioSubtype'

class Channel(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey('Dataset', db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    identifier = models.CharField(max_length=255, db_column='Identifier', blank=True) # Field name made lowercase.
    name = models.CharField(max_length=255, db_column='Name', blank=True) # Field name made lowercase.
    compound = models.ManyToManyField('Chemical', related_name = 'channel_compound', through = 'ChannelWidcompoundWid')
    class Meta:
        db_table = 'Channel'

class ChannelWidcompoundWid(models.Model):
    channelWid = models.ForeignKey(Channel, db_column='ChannelWID', related_name = '+') # Field name made lowercase.
    compoundWid = models.ForeignKey('Chemical', db_column='CompoundWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'ChannelWIDCompoundWID'

class Chemical(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    name = models.CharField(max_length=255, db_column='Name') # Field name made lowercase.
    class_field = models.CharField(max_length=1, db_column='Class', blank=True) # Field name made lowercase. Field renamed because it was a Python reserved word.
    beilsteinname = models.CharField(max_length=50, db_column='BeilsteinName', blank=True) # Field name made lowercase.
    systematicname = models.CharField(max_length=255, db_column='SystematicName', blank=True) # Field name made lowercase.
    cas = models.CharField(max_length=50, db_column='CAS', blank=True) # Field name made lowercase.
    charge = models.SmallIntegerField(null=True, db_column='Charge', blank=True) # Field name made lowercase.
    empiricalformula = models.CharField(max_length=50, db_column='EmpiricalFormula', blank=True) # Field name made lowercase.
    molecularweightcalc = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='MolecularWeightCalc', blank=True) # Field name made lowercase.
    molecularweightexp = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='MolecularWeightExp', blank=True) # Field name made lowercase.
    octh2opartitioncoeff = models.CharField(max_length=50, db_column='OctH2OPartitionCoeff', blank=True) # Field name made lowercase.
    pka1 = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='PKA1', blank=True) # Field name made lowercase.
    pka2 = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='PKA2', blank=True) # Field name made lowercase.
    pka3 = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='PKA3', blank=True) # Field name made lowercase.
    watersolubility = models.CharField(max_length=1, db_column='WaterSolubility', blank=True) # Field name made lowercase.
    smiles = models.CharField(max_length=255, db_column='Smiles', blank=True) # Field name made lowercase.
    datasetWid = models.ForeignKey('Dataset', db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'Chemical'

class Chemicalatom(models.Model):
    chemicalWid = models.ForeignKey(Chemical, db_column='ChemicalWID', related_name = '+') # Field name made lowercase.
    atomindex = models.SmallIntegerField(db_column='AtomIndex') # Field name made lowercase.
    atom = models.CharField(max_length=2, db_column='Atom') # Field name made lowercase.
    charge = models.SmallIntegerField(db_column='Charge') # Field name made lowercase.
    x = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='X', blank=True) # Field name made lowercase.
    y = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='Y', blank=True) # Field name made lowercase.
    z = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='Z', blank=True) # Field name made lowercase.
    stereoparity = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='StereoParity', blank=True) # Field name made lowercase.
    class Meta:
        db_table = 'ChemicalAtom'

class Chemicalbond(models.Model):
    chemicalWid = models.ForeignKey(Chemical, db_column='ChemicalWID', related_name = '+') # Field name made lowercase.
    atom1index = models.SmallIntegerField(db_column='Atom1Index') # Field name made lowercase.
    atom2index = models.SmallIntegerField(db_column='Atom2Index') # Field name made lowercase.
    bondtype = models.SmallIntegerField(db_column='BondType') # Field name made lowercase.
    bondstereo = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='BondStereo', blank=True) # Field name made lowercase.
    class Meta:
        db_table = 'ChemicalBond'

class Citation(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    citation = models.TextField(db_column='Citation', blank=True) # Field name made lowercase.
    pmid = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='PMID', blank=True) # Field name made lowercase.
    title = models.CharField(max_length=255, db_column='Title', blank=True) # Field name made lowercase.
    authors = models.CharField(max_length=255, db_column='Authors', blank=True) # Field name made lowercase.
    publication = models.CharField(max_length=255, db_column='Publication', blank=True) # Field name made lowercase.
    publisher = models.CharField(max_length=255, db_column='Publisher', blank=True) # Field name made lowercase.
    editor = models.CharField(max_length=255, db_column='Editor', blank=True) # Field name made lowercase.
    year = models.CharField(max_length=255, db_column='Year', blank=True) # Field name made lowercase.
    volume = models.CharField(max_length=255, db_column='Volume', blank=True) # Field name made lowercase.
    issue = models.CharField(max_length=255, db_column='Issue', blank=True) # Field name made lowercase.
    pages = models.CharField(max_length=255, db_column='Pages', blank=True) # Field name made lowercase.
    uri = models.CharField(max_length=255, db_column='URI', blank=True) # Field name made lowercase.
    datasetWid = models.ForeignKey('Dataset', db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'Citation'

class CitationWidotherWid(models.Model):
    otherWid = models.BigIntegerField(db_column='OtherWID') # Field name made lowercase.
    citationWid = models.ForeignKey(Citation, db_column='CitationWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'CitationWIDOtherWID'

class Commenttable(models.Model):
    otherWid = models.BigIntegerField(db_column='OtherWID') # Field name made lowercase.
    comm = models.TextField(db_column='Comm', blank=True) # Field name made lowercase.
    class Meta:
        db_table = 'CommentTable'

class ComposgrpWidcompossequenceWid(models.Model):
    compositegroupWid = models.ForeignKey('Designelementgroup', db_column='CompositeGroupWID', related_name = '+') # Field name made lowercase.
    compositesequenceWid = models.ForeignKey('Designelement', db_column='CompositeSequenceWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'ComposGrpWIDComposSequenceWID'

class ComposseqdimensWidcomposseqWid(models.Model):
    compositesequencedimensionWid = models.ForeignKey('Designelementdimension', db_column='CompositeSequenceDimensionWID', related_name = '+') # Field name made lowercase.
    compositesequenceWid = models.ForeignKey('Designelement', db_column='CompositeSequenceWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'ComposSeqDimensWIDComposSeqWID'

class ComposseqWidcomposcomposmapWid(models.Model):
    compositesequenceWid = models.ForeignKey('Designelement', db_column='CompositeSequenceWID', related_name = '+') # Field name made lowercase.
    compositecompositemapWid = models.ForeignKey(Bioevent, db_column='CompositeCompositeMapWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'ComposSeqWIDComposComposMapWID'

class ComposseqWidrepocomposmapWid(models.Model):
    compositesequenceWid = models.ForeignKey('Designelement', db_column='CompositeSequenceWID', related_name = '+') # Field name made lowercase.
    reportercompositemapWid = models.ForeignKey(Bioevent, db_column='ReporterCompositeMapWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'ComposSeqWIDRepoComposMapWID'

class CompositeseqWidbioseqWid(models.Model):
    compositesequenceWid = models.ForeignKey('Designelement', db_column='CompositeSequenceWID', related_name = '+') # Field name made lowercase.
    biosequenceWid = models.BigIntegerField(db_column='BioSequenceWID') # Field name made lowercase.
    class Meta:
        db_table = 'CompositeSeqWIDBioSeqWID'

class Compoundmeasurement(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey('Dataset', db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    compound_componentcompounds = models.ForeignKey(Chemical, null=True, db_column='Compound_ComponentCompounds', blank=True, related_name = '+') # Field name made lowercase.
    compound = models.ForeignKey(Chemical, null=True, db_column='Compound', blank=True, related_name = '+') # Field name made lowercase.
    measurement = models.ForeignKey('Measurement', null=True, db_column='Measurement', blank=True, related_name = '+') # Field name made lowercase.
    treatment_compoundmeasurements = models.ForeignKey(Bioevent, null=True, db_column='Treatment_CompoundMeasurements', blank=True, related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'CompoundMeasurement'

class Computation(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    name = models.CharField(max_length=50, db_column='Name') # Field name made lowercase.
    description = models.TextField(db_column='Description', blank=True) # Field name made lowercase.
    datasetWid = models.ForeignKey('Dataset', db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'Computation'

class Contact(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey('Dataset', db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    mageclass = models.CharField(max_length=100, db_column='MAGEClass') # Field name made lowercase.
    identifier = models.CharField(max_length=255, db_column='Identifier', blank=True) # Field name made lowercase.
    name = models.CharField(max_length=255, db_column='Name', blank=True) # Field name made lowercase.
    uri = models.CharField(max_length=255, db_column='URI', blank=True) # Field name made lowercase.
    address = models.CharField(max_length=255, db_column='Address', blank=True) # Field name made lowercase.
    phone = models.CharField(max_length=255, db_column='Phone', blank=True) # Field name made lowercase.
    tollfreephone = models.CharField(max_length=255, db_column='TollFreePhone', blank=True) # Field name made lowercase.
    email = models.CharField(max_length=255, db_column='Email', blank=True) # Field name made lowercase.
    fax = models.CharField(max_length=255, db_column='Fax', blank=True) # Field name made lowercase.
    parent = models.ForeignKey('self', null=True, db_column='Parent', blank=True, related_name = '+') # Field name made lowercase.
    lastname = models.CharField(max_length=255, db_column='LastName', blank=True) # Field name made lowercase.
    firstname = models.CharField(max_length=255, db_column='FirstName', blank=True) # Field name made lowercase.
    midinitials = models.CharField(max_length=255, db_column='MidInitials', blank=True) # Field name made lowercase.
    affiliation = models.ForeignKey('self', null=True, db_column='Affiliation', blank=True, related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'Contact'

class Crossreference(models.Model):
    otherWid = models.BigIntegerField(db_column='OtherWID') # Field name made lowercase.
    crossWid = models.BigIntegerField(null=True, db_column='CrossWID', blank=True) # Field name made lowercase.
    xid = models.CharField(max_length=50, db_column='XID', blank=True) # Field name made lowercase.
    type = models.CharField(max_length=20, db_column='Type', blank=True) # Field name made lowercase.
    version = models.CharField(max_length=10, db_column='Version', blank=True) # Field name made lowercase.
    relationship = models.CharField(max_length=50, db_column='Relationship', blank=True) # Field name made lowercase.
    databasename = models.CharField(max_length=255, db_column='DataBaseName', blank=True) # Field name made lowercase.
    class Meta:
        db_table = 'CrossReference'

class Dbid(models.Model):
    otherWid = models.BigIntegerField(primary_key = True, db_column='OtherWID') # Field name made lowercase.
    xid = models.CharField(max_length=150, db_column='XID') # Field name made lowercase.
    type = models.CharField(max_length=20, db_column='Type', blank=True) # Field name made lowercase.
    version = models.CharField(max_length=10, db_column='Version', blank=True) # Field name made lowercase.
    class Meta:
        db_table = 'DBID'

class Dataexternal(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey('Dataset', db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    dataformat = models.CharField(max_length=255, db_column='DataFormat', blank=True) # Field name made lowercase.
    dataformatinfouri = models.CharField(max_length=255, db_column='DataFormatInfoURI', blank=True) # Field name made lowercase.
    filenameuri = models.CharField(max_length=255, db_column='FilenameURI', blank=True) # Field name made lowercase.
    class Meta:
        db_table = 'DataExternal'

class Datainternal(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey('Dataset', db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'DataInternal'

class Dataset(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    name = models.CharField(max_length=255, db_column='Name') # Field name made lowercase.
    version = models.CharField(max_length=50, db_column='Version', blank=True) # Field name made lowercase.
    releasedate = models.DateTimeField(null=True, db_column='ReleaseDate', blank=True) # Field name made lowercase.
    loaddate = models.DateTimeField(db_column='LoadDate') # Field name made lowercase.
    changedate = models.DateTimeField(null=True, db_column='ChangeDate', blank=True) # Field name made lowercase.
    homeurl = models.CharField(max_length=255, db_column='HomeURL', blank=True) # Field name made lowercase.
    queryurl = models.CharField(max_length=255, db_column='QueryURL', blank=True) # Field name made lowercase.
    loadedby = models.CharField(max_length=255, db_column='LoadedBy', blank=True) # Field name made lowercase.
    application = models.CharField(max_length=255, db_column='Application', blank=True) # Field name made lowercase.
    applicationversion = models.CharField(max_length=255, db_column='ApplicationVersion', blank=True) # Field name made lowercase.
    class Meta:
        db_table = 'DataSet'

class Datasethierarchy(models.Model):
    superWid = models.ForeignKey(Dataset, db_column='SuperWID', related_name = '+') # Field name made lowercase.
    subWid = models.ForeignKey(Dataset, db_column='SubWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'DataSetHierarchy'

class DatabaseWidcontactWid(models.Model):
    databaseWid = models.ForeignKey('Database', db_column='DatabaseWID', related_name = '+') # Field name made lowercase.
    contactWid = models.ForeignKey(Contact, db_column='ContactWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'DatabaseWIDContactWID'

class Database(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    identifier = models.CharField(max_length=255, db_column='Identifier', blank=True) # Field name made lowercase.
    name = models.CharField(max_length=255, db_column='Name', blank=True) # Field name made lowercase.
    version = models.CharField(max_length=255, db_column='Version', blank=True) # Field name made lowercase.
    uri = models.CharField(max_length=255, db_column='URI', blank=True) # Field name made lowercase.
    contact = models.ManyToManyField('Contact', related_name = 'database_contact', through = 'DatabaseWidcontactWid')
    class Meta:
        db_table = 'Database_'

class Datum(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    value = models.CharField(max_length=255, db_column='Value', blank=True) # Field name made lowercase.
    class Meta:
        db_table = 'Datum'

class DerivbioaWidderivbioadataWid(models.Model):
    derivedbioassayWid = models.ForeignKey(Bioassay, db_column='DerivedBioAssayWID', related_name = '+') # Field name made lowercase.
    derivedbioassaydataWid = models.ForeignKey(Bioassaydata, db_column='DerivedBioAssayDataWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'DerivBioAWIDDerivBioADataWID'

class DerivbioassayWidbioassaymapWid(models.Model):
    derivedbioassayWid = models.ForeignKey(Bioassay, db_column='DerivedBioAssayWID', related_name = '+') # Field name made lowercase.
    bioassaymapWid = models.ForeignKey(Bioevent, db_column='BioAssayMapWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'DerivBioAssayWIDBioAssayMapWID'

class Description(models.Model):
    otherWid = models.BigIntegerField(db_column='OtherWID') # Field name made lowercase.
    tablename = models.CharField(max_length=30, db_column='TableName') # Field name made lowercase.
    comm = models.TextField(db_column='Comm', blank=True) # Field name made lowercase.
    class Meta:
        db_table = 'Description'

class Designelement(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    mageclass = models.CharField(max_length=100, db_column='MAGEClass') # Field name made lowercase.
    identifier = models.CharField(max_length=255, db_column='Identifier', blank=True) # Field name made lowercase.
    name = models.CharField(max_length=255, db_column='Name', blank=True) # Field name made lowercase.
    featuregroup_features = models.ForeignKey('Designelementgroup', null=True, db_column='FeatureGroup_Features', blank=True, related_name = '+') # Field name made lowercase.
    designelement_controltype = models.ForeignKey('Term', null=True, db_column='DesignElement_ControlType', blank=True, related_name = '+') # Field name made lowercase.
    feature_position = models.ForeignKey('Position', null=True, db_column='Feature_Position', blank=True, related_name = '+') # Field name made lowercase.
    zone = models.ForeignKey('Zone', null=True, db_column='Zone', blank=True, related_name = '+') # Field name made lowercase.
    feature_featurelocation = models.ForeignKey('Featurelocation', null=True, db_column='Feature_FeatureLocation', blank=True, related_name = '+') # Field name made lowercase.
    featuregroup = models.ForeignKey('Designelementgroup', null=True, db_column='FeatureGroup', blank=True, related_name = '+') # Field name made lowercase.
    reporter_warningtype = models.ForeignKey('Term', null=True, db_column='Reporter_WarningType', blank=True, related_name = '+') # Field name made lowercase.
    composcomposmap = models.ManyToManyField('Bioevent', related_name = 'composseq_composcomposmap', through = 'ComposseqWidcomposcomposmapWid')
    repocomposmap = models.ManyToManyField('Bioevent', related_name = 'composseq_repocomposmap', through = 'ComposseqWidrepocomposmapWid')
    feature = models.ManyToManyField('Designelement', related_name = 'feature_feature', through = 'FeatureWidfeatureWid')
    feature2 = models.ManyToManyField('Designelement', related_name = 'feature_feature2', through = 'FeatureWidfeatureWid2')
    featurerepormap = models.ManyToManyField('Bioevent', related_name = 'reporter_featurerepormap', through = 'ReporterWidfeaturerepormapWid')
    class Meta:
        db_table = 'DesignElement'

class Designelementdimension(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    mageclass = models.CharField(max_length=100, db_column='MAGEClass') # Field name made lowercase.
    identifier = models.CharField(max_length=255, db_column='Identifier', blank=True) # Field name made lowercase.
    name = models.CharField(max_length=255, db_column='Name', blank=True) # Field name made lowercase.
    composseq = models.ManyToManyField('Designelement', related_name = 'composseqdimens_composseq', through = 'ComposseqdimensWidcomposseqWid')
    feature = models.ManyToManyField('Designelement', related_name = 'featuredimension_feature', through = 'FeaturedimensionWidfeatureWid')
    reporter = models.ManyToManyField('Designelement', related_name = 'reporterdimens_reporter', through = 'ReporterdimensWidreporterWid')
    class Meta:
        db_table = 'DesignElementDimension'

class Designelementgroup(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    mageclass = models.CharField(max_length=100, db_column='MAGEClass') # Field name made lowercase.
    identifier = models.CharField(max_length=255, db_column='Identifier', blank=True) # Field name made lowercase.
    name = models.CharField(max_length=255, db_column='Name', blank=True) # Field name made lowercase.
    arraydesign_featuregroups = models.ForeignKey(Arraydesign, null=True, db_column='ArrayDesign_FeatureGroups', blank=True, related_name = '+') # Field name made lowercase.
    designelementgroup_species = models.ForeignKey('Term', null=True, db_column='DesignElementGroup_Species', blank=True, related_name = '+') # Field name made lowercase.
    featureWidth = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='FeatureWidth', blank=True) # Field name made lowercase.
    featurelength = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='FeatureLength', blank=True) # Field name made lowercase.
    featureheight = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='FeatureHeight', blank=True) # Field name made lowercase.
    featuregroup_technologytype = models.ForeignKey('Term', null=True, db_column='FeatureGroup_TechnologyType', blank=True, related_name = '+') # Field name made lowercase.
    featuregroup_featureshape = models.ForeignKey('Term', null=True, db_column='FeatureGroup_FeatureShape', blank=True, related_name = '+') # Field name made lowercase.
    featuregroup_distanceunit = models.ForeignKey('Unit', null=True, db_column='FeatureGroup_DistanceUnit', blank=True, related_name = '+') # Field name made lowercase.
    compossequence = models.ManyToManyField('Designelement', related_name = 'composgrp_compossequence', through = 'ComposgrpWidcompossequenceWid')
    reporter = models.ManyToManyField('Designelement', related_name = 'reportergroup_reporter', through = 'ReportergroupWidreporterWid')
    class Meta:
        db_table = 'DesignElementGroup'

class Designelementmapping(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    desnelmap = models.ManyToManyField('Bioevent', related_name = 'desnelmapping_desnelmap', through = 'DesnelmappingWiddesnelmapWid')
    class Meta:
        db_table = 'DesignElementMapping'

class Designelementtuple(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    bioassaytuple = models.ForeignKey(Bioassaytuple, null=True, db_column='BioAssayTuple', blank=True, related_name = '+') # Field name made lowercase.
    designelement = models.ForeignKey(Designelement, null=True, db_column='DesignElement', blank=True, related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'DesignElementTuple'

class DesnelmappingWiddesnelmapWid(models.Model):
    designelementmappingWid = models.ForeignKey(Designelementmapping, db_column='DesignElementMappingWID', related_name = '+') # Field name made lowercase.
    designelementmapWid = models.ForeignKey(Bioevent, db_column='DesignElementMapWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'DesnElMappingWIDDesnElMapWID'

class Division(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    code = models.CharField(max_length=10, db_column='Code', blank=True) # Field name made lowercase.
    name = models.CharField(max_length=100, db_column='Name', blank=True) # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'Division'

class Element(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    name = models.CharField(max_length=15, db_column='Name') # Field name made lowercase.
    elementsymbol = models.CharField(max_length=2, db_column='ElementSymbol') # Field name made lowercase.
    atomicweight = models.DecimalField(decimal_places=65535, max_digits=65535, db_column='AtomicWeight') # Field name made lowercase.
    atomicnumber = models.SmallIntegerField(db_column='AtomicNumber') # Field name made lowercase.
    class Meta:
        db_table = 'Element'

class Entry(models.Model):
    otherWid = models.BigIntegerField(primary_key=True, db_column='OtherWID') # Field name made lowercase.
    insertdate = models.DateTimeField(db_column='InsertDate') # Field name made lowercase.
    creationdate = models.DateTimeField(null=True, db_column='CreationDate', blank=True) # Field name made lowercase.
    modifieddate = models.DateTimeField(null=True, db_column='ModifiedDate', blank=True) # Field name made lowercase.
    loaderror = models.CharField(max_length=1, db_column='LoadError') # Field name made lowercase.
    linenumber = models.BigIntegerField(null=True, db_column='LineNumber', blank=True) # Field name made lowercase.
    errormessage = models.TextField(db_column='ErrorMessage', blank=True) # Field name made lowercase.
    datasetWid = models.BigIntegerField(db_column='DataSetWID') # Field name made lowercase.
    class Meta:
        db_table = 'Entry'

class Enumeration(models.Model):
    tablename = models.CharField(max_length=50, db_column='TableName') # Field name made lowercase.
    columnname = models.CharField(max_length=50, db_column='ColumnName') # Field name made lowercase.
    value = models.CharField(max_length=50, db_column='Value') # Field name made lowercase.
    meaning = models.TextField(db_column='Meaning', blank=True) # Field name made lowercase.
    class Meta:
        db_table = 'Enumeration'

class Enzreactionaltcompound(models.Model):
    enzymaticreactionWid = models.ForeignKey('Enzymaticreaction', db_column='EnzymaticReactionWID', related_name = '+') # Field name made lowercase.
    primaryWid = models.ForeignKey(Chemical, db_column='PrimaryWID', related_name = '+') # Field name made lowercase.
    alternativeWid = models.ForeignKey(Chemical, db_column='AlternativeWID', related_name = '+') # Field name made lowercase.
    cofactor = models.CharField(max_length=1, db_column='Cofactor', blank=True) # Field name made lowercase.
    class Meta:
        db_table = 'EnzReactionAltCompound'

class Enzreactioncofactor(models.Model):
    enzymaticreactionWid = models.ForeignKey('Enzymaticreaction', db_column='EnzymaticReactionWID', related_name = '+') # Field name made lowercase.
    chemicalWid = models.ForeignKey(Chemical, db_column='ChemicalWID', related_name = '+') # Field name made lowercase.
    prosthetic = models.CharField(max_length=1, db_column='Prosthetic', blank=True) # Field name made lowercase.
    class Meta:
        db_table = 'EnzReactionCofactor'

class Enzreactioninhibitoractivator(models.Model):
    enzymaticreactionWid = models.ForeignKey('Enzymaticreaction', db_column='EnzymaticReactionWID', related_name = '+') # Field name made lowercase.
    compoundWid = models.BigIntegerField(db_column='CompoundWID') # Field name made lowercase.
    inhibitoractivate = models.CharField(max_length=1, db_column='InhibitOrActivate', blank=True) # Field name made lowercase.
    mechanism = models.CharField(max_length=1, db_column='Mechanism', blank=True) # Field name made lowercase.
    physiorelevant = models.CharField(max_length=1, db_column='PhysioRelevant', blank=True) # Field name made lowercase.
    class Meta:
        db_table = 'EnzReactionInhibitorActivator'

class Enzymaticreaction(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    reactionWid = models.ForeignKey('Reaction', db_column='ReactionWID', related_name = '+') # Field name made lowercase.
    proteinWid = models.ForeignKey('Protein', db_column='ProteinWID', related_name = '+') # Field name made lowercase.
    complexWid = models.ForeignKey('Protein', null=True, db_column='ComplexWID', blank=True, related_name = '+') # Field name made lowercase.
    reactiondirection = models.CharField(max_length=30, db_column='ReactionDirection', blank=True) # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'EnzymaticReaction'

class ExperimWidbioassaydataclustWid(models.Model):
    experimentWid = models.ForeignKey('Experiment', db_column='ExperimentWID', related_name = '+') # Field name made lowercase.
    bioassaydataclusterWid = models.ForeignKey(Bioassaydatacluster, db_column='BioAssayDataClusterWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'ExperimWIDBioAssayDataClustWID'

class Experiment(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    type = models.CharField(max_length=50, db_column='Type') # Field name made lowercase.
    contactWid = models.ForeignKey(Contact, null=True, db_column='ContactWID', blank=True, related_name = '+') # Field name made lowercase.
    archiveWid = models.ForeignKey(Archive, null=True, db_column='ArchiveWID', blank=True, related_name = '+') # Field name made lowercase.
    startdate = models.DateTimeField(null=True, db_column='StartDate', blank=True) # Field name made lowercase.
    enddate = models.DateTimeField(null=True, db_column='EndDate', blank=True) # Field name made lowercase.
    description = models.TextField(db_column='Description', blank=True) # Field name made lowercase.
    groupWid = models.ForeignKey('self', null=True, db_column='GroupWID', blank=True, related_name = '+') # Field name made lowercase.
    grouptype = models.CharField(max_length=50, db_column='GroupType', blank=True) # Field name made lowercase.
    groupsize = models.IntegerField(db_column='GroupSize') # Field name made lowercase.
    groupindex = models.BigIntegerField(null=True, db_column='GroupIndex', blank=True) # Field name made lowercase.
    timepoint = models.BigIntegerField(null=True, db_column='TimePoint', blank=True) # Field name made lowercase.
    timeunit = models.CharField(max_length=20, db_column='TimeUnit', blank=True) # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    biosourceWid = models.ForeignKey(Biosource, null=True, db_column='BioSourceWID', blank=True, related_name = '+') # Field name made lowercase.
    bioassaydataclust = models.ManyToManyField('Bioassaydatacluster', related_name = 'experim_bioassaydataclust', through = 'ExperimWidbioassaydataclustWid')
    bioassaydata = models.ManyToManyField('Bioassaydata', related_name = 'experiment_bioassaydata', through = 'ExperimentWidbioassaydataWid')
    bioassay = models.ManyToManyField('Bioassay', related_name = 'experiment_bioassay', through = 'ExperimentWidbioassayWid')
    contact = models.ManyToManyField('Contact', related_name = 'experiment_contact', through = 'ExperimentWidcontactWid')
    class Meta:
        db_table = 'Experiment'

class Experimentdata(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    experimentWid = models.ForeignKey(Experiment, db_column='ExperimentWID', related_name = '+') # Field name made lowercase.
    data = models.TextField(db_column='Data', blank=True) # Field name made lowercase.
    magedata = models.ForeignKey('Parametervalue', null=True, db_column='MageData', blank=True, related_name = '+') # Field name made lowercase.
    role = models.CharField(max_length=50, db_column='Role') # Field name made lowercase.
    kind = models.CharField(max_length=1, db_column='Kind') # Field name made lowercase.
    dateproduced = models.DateTimeField(null=True, db_column='DateProduced', blank=True) # Field name made lowercase.
    otherWid = models.BigIntegerField(null=True, db_column='OtherWID', blank=True) # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    ##bioassay = models.ManyToManyField('Bioassay', related_name = 'experimentdesign_bioassay', through = 'ExperimentdesignWidbioassayWid')
    class Meta:
        db_table = 'ExperimentData'

class Experimentdesign(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    experiment_experimentdesigns = models.ForeignKey(Experiment, null=True, db_column='Experiment_ExperimentDesigns', blank=True, related_name = '+') # Field name made lowercase.
    qualitycontroldescription = models.BigIntegerField(null=True, db_column='QualityControlDescription', blank=True) # Field name made lowercase.
    normalizationdescription = models.BigIntegerField(null=True, db_column='NormalizationDescription', blank=True) # Field name made lowercase.
    replicatedescription = models.BigIntegerField(null=True, db_column='ReplicateDescription', blank=True) # Field name made lowercase.
    class Meta:
        db_table = 'ExperimentDesign'

class ExperimentdesignWidbioassayWid(models.Model):
    experimentdesignWid = models.ForeignKey(Experimentdesign, db_column='ExperimentDesignWID', related_name = '+') # Field name made lowercase.
    bioassayWid = models.ForeignKey(Bioassay, db_column='BioAssayWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'ExperimentDesignWIDBioAssayWID'

class Experimentrelationship(models.Model):
    experimentWid = models.ForeignKey(Experiment, db_column='ExperimentWID', related_name = '+') # Field name made lowercase.
    relatedexperimentWid = models.ForeignKey(Experiment, db_column='RelatedExperimentWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'ExperimentRelationship'

class ExperimentWidbioassaydataWid(models.Model):
    experimentWid = models.ForeignKey(Experiment, db_column='ExperimentWID', related_name = '+') # Field name made lowercase.
    bioassaydataWid = models.ForeignKey(Bioassaydata, db_column='BioAssayDataWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'ExperimentWIDBioAssayDataWID'

class ExperimentWidbioassayWid(models.Model):
    experimentWid = models.ForeignKey(Experiment, db_column='ExperimentWID', related_name = '+') # Field name made lowercase.
    bioassayWid = models.ForeignKey(Bioassay, db_column='BioAssayWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'ExperimentWIDBioAssayWID'

class ExperimentWidcontactWid(models.Model):
    experimentWid = models.ForeignKey(Experiment, db_column='ExperimentWID', related_name = '+') # Field name made lowercase.
    contactWid = models.ForeignKey(Contact, db_column='ContactWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'ExperimentWIDContactWID'

class Experimentalfactor(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    identifier = models.CharField(max_length=255, db_column='Identifier', blank=True) # Field name made lowercase.
    name = models.CharField(max_length=255, db_column='Name', blank=True) # Field name made lowercase.
    experimentdesign = models.ForeignKey(Experimentdesign, null=True, db_column='ExperimentDesign', blank=True, related_name = '+') # Field name made lowercase.
    experimentalfactor_category = models.ForeignKey('Term', null=True, db_column='ExperimentalFactor_Category', blank=True, related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'ExperimentalFactor'

class Factorvalue(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    identifier = models.CharField(max_length=255, db_column='Identifier', blank=True) # Field name made lowercase.
    name = models.CharField(max_length=255, db_column='Name', blank=True) # Field name made lowercase.
    experimentalfactor = models.ForeignKey(Experimentalfactor, null=True, db_column='ExperimentalFactor', blank=True, related_name = '+') # Field name made lowercase.
    experimentalfactor2 = models.ForeignKey(Experimentalfactor, null=True, db_column='ExperimentalFactor2', blank=True, related_name = '+') # Field name made lowercase.
    factorvalue_measurement = models.ForeignKey('Measurement', null=True, db_column='FactorValue_Measurement', blank=True, related_name = '+') # Field name made lowercase.
    factorvalue_value = models.ForeignKey('Term', null=True, db_column='FactorValue_Value', blank=True, related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'FactorValue'

class Feature(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    description = models.CharField(max_length=1300, db_column='Description', blank=True) # Field name made lowercase.
    type = models.CharField(max_length=50, db_column='Type', blank=True) # Field name made lowercase.
    class_field = models.CharField(max_length=50, db_column='Class', blank=True) # Field name made lowercase. Field renamed because it was a Python reserved word.
    sequencetype = models.CharField(max_length=1, db_column='SequenceType') # Field name made lowercase.
    sequenceWid = models.BigIntegerField(null=True, db_column='SequenceWID', blank=True) # Field name made lowercase.
    variant = models.TextField(db_column='Variant', blank=True) # Field name made lowercase.
    regionorpoint = models.CharField(max_length=10, db_column='RegionOrPoint', blank=True) # Field name made lowercase.
    pointtype = models.CharField(max_length=10, db_column='PointType', blank=True) # Field name made lowercase.
    startposition = models.BigIntegerField(null=True, db_column='StartPosition', blank=True) # Field name made lowercase.
    endposition = models.BigIntegerField(null=True, db_column='EndPosition', blank=True) # Field name made lowercase.
    startpositionapproximate = models.CharField(max_length=10, db_column='StartPositionApproximate', blank=True) # Field name made lowercase.
    endpositionapproximate = models.CharField(max_length=10, db_column='EndPositionApproximate', blank=True) # Field name made lowercase.
    experimentalsupport = models.CharField(max_length=1, db_column='ExperimentalSupport', blank=True) # Field name made lowercase.
    computationalsupport = models.CharField(max_length=1, db_column='ComputationalSupport', blank=True) # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'Feature'

class Featuredefect(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    arraymanufacturedeviation = models.ForeignKey(Arraymanufacturedeviation, null=True, db_column='ArrayManufactureDeviation', blank=True, related_name = '+') # Field name made lowercase.
    featuredefect_defecttype = models.ForeignKey('Term', null=True, db_column='FeatureDefect_DefectType', blank=True, related_name = '+') # Field name made lowercase.
    featuredefect_positiondelta = models.ForeignKey('Positiondelta', null=True, db_column='FeatureDefect_PositionDelta', blank=True, related_name = '+') # Field name made lowercase.
    feature = models.ForeignKey(Designelement, null=True, db_column='Feature', blank=True, related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'FeatureDefect'

class FeaturedimensionWidfeatureWid(models.Model):
    featuredimensionWid = models.ForeignKey(Designelementdimension, db_column='FeatureDimensionWID', related_name = '+') # Field name made lowercase.
    featureWid = models.ForeignKey(Designelement, db_column='FeatureWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'FeatureDimensionWIDFeatureWID'

class Featureinformation(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    feature = models.ForeignKey(Designelement, null=True, db_column='Feature', blank=True, related_name = '+') # Field name made lowercase.
    featurereportermap = models.ForeignKey(Bioevent, null=True, db_column='FeatureReporterMap', blank=True, related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'FeatureInformation'

class Featurelocation(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    row_field = models.SmallIntegerField(null=True, db_column='Row_', blank=True) # Field name made lowercase. Field renamed because it ended with '_'.
    column_field = models.SmallIntegerField(null=True, db_column='Column_', blank=True) # Field name made lowercase. Field renamed because it ended with '_'.
    class Meta:
        db_table = 'FeatureLocation'

class FeatureWidfeatureWid(models.Model):
    featureWid1 = models.ForeignKey(Designelement, db_column='FeatureWID1', related_name = '+') # Field name made lowercase.
    featureWid2 = models.ForeignKey(Designelement, db_column='FeatureWID2', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'FeatureWIDFeatureWID'

class FeatureWidfeatureWid2(models.Model):
    featureWid1 = models.ForeignKey(Designelement, db_column='FeatureWID1', related_name = '+') # Field name made lowercase.
    featureWid2 = models.ForeignKey(Designelement, db_column='FeatureWID2', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'FeatureWIDFeatureWID2'

class Fiducial(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    arraygroup_fiducials = models.ForeignKey(Arraygroup, null=True, db_column='ArrayGroup_Fiducials', blank=True, related_name = '+') # Field name made lowercase.
    fiducial_fiducialtype = models.ForeignKey('Term', null=True, db_column='Fiducial_FiducialType', blank=True, related_name = '+') # Field name made lowercase.
    fiducial_distanceunit = models.ForeignKey('Unit', null=True, db_column='Fiducial_DistanceUnit', blank=True, related_name = '+') # Field name made lowercase.
    fiducial_position = models.ForeignKey('Position', null=True, db_column='Fiducial_Position', blank=True, related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'Fiducial'

class Flowcytometryprobe(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    type = models.CharField(max_length=100, db_column='Type') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'FlowCytometryProbe'

class Flowcytometrysample(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    biosourceWid = models.ForeignKey(Biosource, null=True, db_column='BioSourceWID', blank=True, related_name = '+') # Field name made lowercase.
    flowcytometryprobeWid = models.ForeignKey(Flowcytometryprobe, null=True, db_column='FlowCytometryProbeWID', blank=True, related_name = '+') # Field name made lowercase.
    measurementWid = models.ForeignKey('Measurement', null=True, db_column='MeasurementWID', blank=True, related_name = '+') # Field name made lowercase.
    manufacturerWid = models.ForeignKey(Contact, null=True, db_column='ManufacturerWID', blank=True, related_name = '+') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'FlowCytometrySample'

class Function(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    name = models.CharField(max_length=255, db_column='Name', blank=True) # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'Function'

class Gellocation(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    spotWid = models.ForeignKey('Spot', db_column='SpotWID', related_name = '+') # Field name made lowercase.
    xcoord = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='Xcoord', blank=True) # Field name made lowercase.
    ycoord = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='Ycoord', blank=True) # Field name made lowercase.
    refgel = models.CharField(max_length=1, db_column='refGel', blank=True) # Field name made lowercase.
    experimentWid = models.ForeignKey(Experiment, db_column='ExperimentWID', related_name = '+') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DatasetWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'GelLocation'

class Gene(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    name = models.CharField(max_length=255, db_column='Name', blank=True) # Field name made lowercase.
    nucleicacidWid = models.ForeignKey('Nucleicacid', null=True, db_column='NucleicAcidWID', blank=True, related_name = '+') # Field name made lowercase.
    subsequenceWid = models.BigIntegerField(null=True, db_column='SubsequenceWID', blank=True) # Field name made lowercase.
    type = models.CharField(max_length=100, db_column='Type', blank=True) # Field name made lowercase.
    genomeid = models.CharField(max_length=35, db_column='GenomeID', blank=True) # Field name made lowercase.
    codingregionstart = models.BigIntegerField(null=True, db_column='CodingRegionStart', blank=True) # Field name made lowercase.
    codingregionend = models.BigIntegerField(null=True, db_column='CodingRegionEnd', blank=True) # Field name made lowercase.
    codingregionstartapproximate = models.CharField(max_length=10, db_column='CodingRegionStartApproximate', blank=True) # Field name made lowercase.
    codingregionendapproximate = models.CharField(max_length=10, db_column='CodingRegionEndApproximate', blank=True) # Field name made lowercase.
    direction = models.CharField(max_length=25, db_column='Direction', blank=True) # Field name made lowercase.
    interrupted = models.CharField(max_length=1, db_column='Interrupted', blank=True) # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    nucleicacid = models.ManyToManyField('Nucleicacid', related_name = 'gene_nucleicacid', through = 'GeneWidnucleicacidWid')
    protein = models.ManyToManyField('Protein', related_name = 'gene_protein', through = 'GeneWidproteinWid')
    class Meta:
        db_table = 'Gene'

class Geneexpressiondata(models.Model):
    b = models.SmallIntegerField(db_column='B') # Field name made lowercase.
    d = models.SmallIntegerField(db_column='D') # Field name made lowercase.
    q = models.SmallIntegerField(db_column='Q') # Field name made lowercase.
    value = models.CharField(max_length=100, db_column='Value') # Field name made lowercase.
    bioassayvaluesWid = models.ForeignKey(Biodatavalues, db_column='BioAssayValuesWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'GeneExpressionData'

class GeneWidnucleicacidWid(models.Model):
    geneWid = models.ForeignKey(Gene, db_column='GeneWID', related_name = '+') # Field name made lowercase.
    nucleicacidWid = models.ForeignKey('Nucleicacid', db_column='NucleicAcidWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'GeneWIDNucleicAcidWID'

class GeneWidproteinWid(models.Model):
    geneWid = models.ForeignKey(Gene, db_column='GeneWID', related_name = '+') # Field name made lowercase.
    proteinWid = models.ForeignKey('Protein', db_column='ProteinWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'GeneWIDProteinWID'

class Geneticcode(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    ncbiid = models.CharField(max_length=2, db_column='NCBIID', blank=True) # Field name made lowercase.
    name = models.CharField(max_length=100, db_column='Name', blank=True) # Field name made lowercase.
    translationtable = models.CharField(max_length=64, db_column='TranslationTable', blank=True) # Field name made lowercase.
    startcodon = models.CharField(max_length=64, db_column='StartCodon', blank=True) # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'GeneticCode'

class HardwareWidcontactWid(models.Model):
    hardwareWid = models.ForeignKey('Parameterizable', db_column='HardwareWID', related_name = '+') # Field name made lowercase.
    contactWid = models.ForeignKey(Contact, db_column='ContactWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'HardwareWIDContactWID'

class HardwareWidsoftwareWid(models.Model):
    hardwareWid = models.ForeignKey('Parameterizable', db_column='HardwareWID', related_name = '+') # Field name made lowercase.
    softwareWid = models.ForeignKey('Parameterizable', db_column='SoftwareWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'HardwareWIDSoftwareWID'

class Image(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    identifier = models.CharField(max_length=255, db_column='Identifier', blank=True) # Field name made lowercase.
    name = models.CharField(max_length=255, db_column='Name', blank=True) # Field name made lowercase.
    uri = models.CharField(max_length=255, db_column='URI', blank=True) # Field name made lowercase.
    image_format = models.ForeignKey('Term', null=True, db_column='Image_Format', blank=True, related_name = '+') # Field name made lowercase.
    physicalbioassay = models.ForeignKey(Bioassay, null=True, db_column='PhysicalBioAssay', blank=True, related_name = '+') # Field name made lowercase.
    channel = models.ManyToManyField('Channel', related_name = 'image_channel', through = 'ImageWidchannelWid')
    class Meta:
        db_table = 'Image'

class ImageacquisitionWidimageWid(models.Model):
    imageacquisitionWid = models.ForeignKey(Bioevent, db_column='ImageAcquisitionWID', related_name = '+') # Field name made lowercase.
    imageWid = models.ForeignKey(Image, db_column='ImageWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'ImageAcquisitionWIDImageWID'

class ImageWidchannelWid(models.Model):
    imageWid = models.ForeignKey(Image, db_column='ImageWID', related_name = '+') # Field name made lowercase.
    channelWid = models.ForeignKey(Channel, db_column='ChannelWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'ImageWIDChannelWID'

class Interaction(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    type = models.CharField(max_length=100, db_column='Type', blank=True) # Field name made lowercase.
    name = models.CharField(max_length=120, db_column='Name', blank=True) # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'Interaction'

class Interactionparticipant(models.Model):
    interactionWid = models.ForeignKey(Interaction, db_column='InteractionWID', related_name = '+') # Field name made lowercase.
    otherWid = models.BigIntegerField(db_column='OtherWID') # Field name made lowercase.
    coefficient = models.SmallIntegerField(null=True, db_column='Coefficient', blank=True) # Field name made lowercase.
    class Meta:
        db_table = 'InteractionParticipant'

class LabeledextractWidcompoundWid(models.Model):
    labeledextractWid = models.ForeignKey(Biosource, db_column='LabeledExtractWID', related_name = '+') # Field name made lowercase.
    compoundWid = models.ForeignKey(Chemical, db_column='CompoundWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'LabeledExtractWIDCompoundWID'

class Lightsource(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    wavelength = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='Wavelength', blank=True) # Field name made lowercase.
    type = models.CharField(max_length=100, db_column='Type', blank=True) # Field name made lowercase.
    instrumentWid = models.BigIntegerField(null=True, db_column='InstrumentWID', blank=True) # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'LightSource'

class Location(models.Model):
    proteinWid = models.ForeignKey('Protein', db_column='ProteinWID', related_name = '+') # Field name made lowercase.
    location = models.CharField(max_length=100, db_column='Location') # Field name made lowercase.
    class Meta:
        db_table = 'Location'

class Manufacturelims(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    mageclass = models.CharField(max_length=100, db_column='MAGEClass') # Field name made lowercase.
    arraymanufacture_featurelimss = models.ForeignKey(Arraymanufacture, null=True, db_column='ArrayManufacture_FeatureLIMSs', blank=True, related_name = '+') # Field name made lowercase.
    quality = models.CharField(max_length=255, db_column='Quality', blank=True) # Field name made lowercase.
    feature = models.ForeignKey(Designelement, null=True, db_column='Feature', blank=True, related_name = '+') # Field name made lowercase.
    biomaterial = models.ForeignKey(Biosource, null=True, db_column='BioMaterial', blank=True, related_name = '+') # Field name made lowercase.
    manufacturelims_identifierlims = models.BigIntegerField(null=True, db_column='ManufactureLIMS_IdentifierLIMS', blank=True) # Field name made lowercase.
    biomaterialplateidentifier = models.CharField(max_length=255, db_column='BioMaterialPlateIdentifier', blank=True) # Field name made lowercase.
    biomaterialplaterow = models.CharField(max_length=255, db_column='BioMaterialPlateRow', blank=True) # Field name made lowercase.
    biomaterialplatecol = models.CharField(max_length=255, db_column='BioMaterialPlateCol', blank=True) # Field name made lowercase.
    class Meta:
        db_table = 'ManufactureLIMS'

class MeasbassayWidmeasbassaydataWid(models.Model):
    measuredbioassayWid = models.ForeignKey(Bioassay, db_column='MeasuredBioAssayWID', related_name = '+') # Field name made lowercase.
    measuredbioassaydataWid = models.ForeignKey(Bioassaydata, db_column='MeasuredBioAssayDataWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'MeasBAssayWIDMeasBAssayDataWID'

class Measurement(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    type_field = models.CharField(max_length=25, db_column='Type_', blank=True) # Field name made lowercase. Field renamed because it ended with '_'.
    value = models.CharField(max_length=255, db_column='Value', blank=True) # Field name made lowercase.
    kindcv = models.CharField(max_length=25, db_column='KindCV', blank=True) # Field name made lowercase.
    otherkind = models.CharField(max_length=255, db_column='OtherKind', blank=True) # Field name made lowercase.
    measurement_unit = models.ForeignKey('Unit', null=True, db_column='Measurement_Unit', blank=True, related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'Measurement'

class Mismatchinformation(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    compositeposition = models.ForeignKey('Sequenceposition', null=True, db_column='CompositePosition', blank=True, related_name = '+') # Field name made lowercase.
    featureinformation = models.ForeignKey(Featureinformation, null=True, db_column='FeatureInformation', blank=True, related_name = '+') # Field name made lowercase.
    startcoord = models.SmallIntegerField(null=True, db_column='StartCoord', blank=True) # Field name made lowercase.
    newsequence = models.CharField(max_length=255, db_column='NewSequence', blank=True) # Field name made lowercase.
    replacedlength = models.SmallIntegerField(null=True, db_column='ReplacedLength', blank=True) # Field name made lowercase.
    reporterposition = models.ForeignKey('Sequenceposition', null=True, db_column='ReporterPosition', blank=True, related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'MismatchInformation'

class Namevaluetype(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    name = models.CharField(max_length=255, db_column='Name', blank=True) # Field name made lowercase.
    value = models.CharField(max_length=255, db_column='Value', blank=True) # Field name made lowercase.
    type_field = models.CharField(max_length=255, db_column='Type_', blank=True) # Field name made lowercase. Field renamed because it ended with '_'.
    namevaluetype_propertysets = models.ForeignKey('self', null=True, db_column='NameValueType_PropertySets', blank=True, related_name = '+') # Field name made lowercase.
    otherWid = models.BigIntegerField(null=True, db_column='OtherWID', blank=True) # Field name made lowercase.
    class Meta:
        db_table = 'NameValueType'

class Node(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    bioassaydatacluster_nodes = models.ForeignKey(Bioassaydatacluster, null=True, db_column='BioAssayDataCluster_Nodes', blank=True, related_name = '+') # Field name made lowercase.
    node_nodes = models.ForeignKey('self', null=True, db_column='Node_Nodes', blank=True, related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'Node'

class Nodecontents(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    node_nodecontents = models.ForeignKey(Node, null=True, db_column='Node_NodeContents', blank=True, related_name = '+') # Field name made lowercase.
    bioassaydimension = models.ForeignKey(Bioassaydimension, null=True, db_column='BioAssayDimension', blank=True, related_name = '+') # Field name made lowercase.
    designelementdimension = models.ForeignKey(Designelementdimension, null=True, db_column='DesignElementDimension', blank=True, related_name = '+') # Field name made lowercase.
    quantitationdimension = models.ForeignKey('Quantitationtypedimension', null=True, db_column='QuantitationDimension', blank=True, related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'NodeContents'

class Nodevalue(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    node_nodevalue = models.ForeignKey(Node, null=True, db_column='Node_NodeValue', blank=True, related_name = '+') # Field name made lowercase.
    name = models.CharField(max_length=255, db_column='Name', blank=True) # Field name made lowercase.
    value = models.CharField(max_length=255, db_column='Value', blank=True) # Field name made lowercase.
    nodevalue_type = models.ForeignKey('Term', null=True, db_column='NodeValue_Type', blank=True, related_name = '+') # Field name made lowercase.
    nodevalue_scale = models.ForeignKey('Term', null=True, db_column='NodeValue_Scale', blank=True, related_name = '+') # Field name made lowercase.
    nodevalue_datatype = models.ForeignKey('Term', null=True, db_column='NodeValue_DataType', blank=True, related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'NodeValue'

class Nucleicacid(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    name = models.CharField(max_length=200, db_column='Name', blank=True) # Field name made lowercase.
    type = models.CharField(max_length=30, db_column='Type') # Field name made lowercase.
    class_field = models.CharField(max_length=30, db_column='Class', blank=True) # Field name made lowercase. Field renamed because it was a Python reserved word.
    topology = models.CharField(max_length=30, db_column='Topology', blank=True) # Field name made lowercase.
    strandedness = models.CharField(max_length=30, db_column='Strandedness', blank=True) # Field name made lowercase.
    sequencederivation = models.CharField(max_length=30, db_column='SequenceDerivation', blank=True) # Field name made lowercase.
    fragment = models.CharField(max_length=1, db_column='Fragment', blank=True) # Field name made lowercase.
    fullysequenced = models.CharField(max_length=1, db_column='FullySequenced', blank=True) # Field name made lowercase.
    moleculelength = models.BigIntegerField(null=True, db_column='MoleculeLength', blank=True) # Field name made lowercase.
    moleculelengthapproximate = models.CharField(max_length=10, db_column='MoleculeLengthApproximate', blank=True) # Field name made lowercase.
    cumulativelength = models.BigIntegerField(null=True, db_column='CumulativeLength', blank=True) # Field name made lowercase.
    cumulativelengthapproximate = models.CharField(max_length=10, db_column='CumulativeLengthApproximate', blank=True) # Field name made lowercase.
    geneticcodeWid = models.ForeignKey(Geneticcode, null=True, db_column='GeneticCodeWID', blank=True, related_name = '+') # Field name made lowercase.
    biosourceWid = models.ForeignKey(Biosource, null=True, db_column='BioSourceWID', blank=True, related_name = '+') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'NucleicAcid'

class Parameter(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    identifier = models.CharField(max_length=255, db_column='Identifier', blank=True) # Field name made lowercase.
    name = models.CharField(max_length=255, db_column='Name', blank=True) # Field name made lowercase.
    parameter_defaultvalue = models.ForeignKey(Measurement, null=True, db_column='Parameter_DefaultValue', blank=True, related_name = '+') # Field name made lowercase.
    parameter_datatype = models.ForeignKey('Term', null=True, db_column='Parameter_DataType', blank=True, related_name = '+') # Field name made lowercase.
    parameterizable_parametertypes = models.ForeignKey('Parameterizable', null=True, db_column='Parameterizable_ParameterTypes', blank=True, related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'Parameter'

class Parametervalue(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    value = models.CharField(max_length=255, db_column='Value', blank=True) # Field name made lowercase.
    parametertype = models.ForeignKey(Parameter, null=True, db_column='ParameterType', blank=True, related_name = '+') # Field name made lowercase.
    parameterizableapplication = models.ForeignKey('Parameterizableapplication', null=True, db_column='ParameterizableApplication', blank=True, related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'ParameterValue'

class Parameterizable(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    mageclass = models.CharField(max_length=100, db_column='MAGEClass') # Field name made lowercase.
    identifier = models.CharField(max_length=255, db_column='Identifier', blank=True) # Field name made lowercase.
    name = models.CharField(max_length=255, db_column='Name', blank=True) # Field name made lowercase.
    uri = models.CharField(max_length=255, db_column='URI', blank=True) # Field name made lowercase.
    model = models.CharField(max_length=255, db_column='Model', blank=True) # Field name made lowercase.
    make = models.CharField(max_length=255, db_column='Make', blank=True) # Field name made lowercase.
    hardware_type = models.ForeignKey('Term', null=True, db_column='Hardware_Type', blank=True, related_name = '+') # Field name made lowercase.
    text = models.CharField(max_length=1000, db_column='Text', blank=True) # Field name made lowercase.
    title = models.CharField(max_length=255, db_column='Title', blank=True) # Field name made lowercase.
    protocol_type = models.ForeignKey('Term', null=True, db_column='Protocol_Type', blank=True, related_name = '+') # Field name made lowercase.
    software_type = models.ForeignKey('Term', null=True, db_column='Software_Type', blank=True, related_name = '+') # Field name made lowercase.
    hardware = models.ForeignKey('self', null=True, db_column='Hardware', blank=True, related_name = '+') # Field name made lowercase.
    contact = models.ManyToManyField('Contact', related_name = 'hardware_contact', through = 'HardwareWidcontactWid')
    hardware_software = models.ManyToManyField('Parameterizable', related_name = 'hardware_software_key', through = 'HardwareWidsoftwareWid')
    protocol_hardware = models.ManyToManyField('Parameterizable', related_name = 'protocol_hardware_key', through = 'ProtocolWidhardwareWid')
    protocol_software = models.ManyToManyField('Parameterizable', related_name = 'protocol_software_key', through = 'ProtocolWidsoftwareWid')
    software_contact = models.ManyToManyField('Contact', related_name = 'software_contact_key', through = 'SoftwareWidcontactWid')
    software_software = models.ManyToManyField('Parameterizable', related_name = 'software_software_key', through = 'SoftwareWidsoftwareWid')
    class Meta:
        db_table = 'Parameterizable'

class Parameterizableapplication(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    mageclass = models.CharField(max_length=100, db_column='MAGEClass') # Field name made lowercase.
    arraydesign = models.ForeignKey(Arraydesign, null=True, db_column='ArrayDesign', blank=True, related_name = '+') # Field name made lowercase.
    arraymanufacture = models.ForeignKey(Arraymanufacture, null=True, db_column='ArrayManufacture', blank=True, related_name = '+') # Field name made lowercase.
    bioevent_protocolapplications = models.ForeignKey(Bioevent, null=True, db_column='BioEvent_ProtocolApplications', blank=True, related_name = '+') # Field name made lowercase.
    serialnumber = models.CharField(max_length=255, db_column='SerialNumber', blank=True) # Field name made lowercase.
    hardware = models.ForeignKey(Parameterizable, null=True, db_column='Hardware', blank=True, related_name = '+') # Field name made lowercase.
    activitydate = models.CharField(max_length=255, db_column='ActivityDate', blank=True) # Field name made lowercase.
    protocolapplication = models.ForeignKey('self', null=True, db_column='ProtocolApplication', blank=True, related_name = '+') # Field name made lowercase.
    protocolapplication2 = models.ForeignKey('self', null=True, db_column='ProtocolApplication2', blank=True, related_name = '+') # Field name made lowercase.
    protocol = models.ForeignKey(Parameterizable, null=True, db_column='Protocol', blank=True, related_name = '+') # Field name made lowercase.
    version = models.CharField(max_length=255, db_column='Version', blank=True) # Field name made lowercase.
    releasedate = models.DateTimeField(null=True, db_column='ReleaseDate', blank=True) # Field name made lowercase.
    software = models.ForeignKey(Parameterizable, null=True, db_column='Software', blank=True, related_name = '+') # Field name made lowercase.
    person = models.ManyToManyField('Contact', related_name = 'protocolappl_person', through = 'ProtocolapplWidpersonWid')
    class Meta:
        db_table = 'ParameterizableApplication'

class Pathway(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    name = models.CharField(max_length=255, db_column='Name') # Field name made lowercase.
    type = models.CharField(max_length=1, db_column='Type') # Field name made lowercase.
    biosourceWid = models.ForeignKey(Biosource, null=True, db_column='BioSourceWID', blank=True, related_name = '+') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'Pathway'

class Pathwaylink(models.Model):
    pathway1Wid = models.ForeignKey(Pathway, db_column='Pathway1WID', related_name = '+') # Field name made lowercase.
    pathway2Wid = models.ForeignKey(Pathway, db_column='Pathway2WID', related_name = '+') # Field name made lowercase.
    chemicalWid = models.ForeignKey(Chemical, db_column='ChemicalWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'PathwayLink'

class Pathwayreaction(models.Model):
    pathwayWid = models.ForeignKey(Pathway, db_column='PathwayWID', related_name = '+') # Field name made lowercase.
    reactionWid = models.BigIntegerField(db_column='ReactionWID') # Field name made lowercase.
    priorreactionWid = models.ForeignKey('Reaction', null=True, db_column='PriorReactionWID', blank=True, related_name = '+') # Field name made lowercase.
    hypothetical = models.CharField(max_length=1, db_column='Hypothetical') # Field name made lowercase.
    class Meta:
        db_table = 'PathwayReaction'

class Positiondelta(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    deltax = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='DeltaX', blank=True) # Field name made lowercase.
    deltay = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='DeltaY', blank=True) # Field name made lowercase.
    positiondelta_distanceunit = models.ForeignKey('Unit', null=True, db_column='PositionDelta_DistanceUnit', blank=True, related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'PositionDelta'

class Position(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    x = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='X', blank=True) # Field name made lowercase.
    y = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='Y', blank=True) # Field name made lowercase.
    position_distanceunit = models.ForeignKey('Unit', null=True, db_column='Position_DistanceUnit', blank=True, related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'Position_'

class Product(models.Model):
    reactionWid = models.ForeignKey('Reaction', db_column='ReactionWID', related_name = '+') # Field name made lowercase.
    otherWid = models.BigIntegerField(db_column='OtherWID') # Field name made lowercase.
    coefficient = models.SmallIntegerField(db_column='Coefficient') # Field name made lowercase.
    class Meta:
        db_table = 'Product'

class Protein(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    name = models.TextField(db_column='Name', blank=True) # Field name made lowercase.
    aasequence = models.TextField(db_column='AASequence', blank=True) # Field name made lowercase.
    length = models.BigIntegerField(null=True, db_column='Length', blank=True) # Field name made lowercase.
    lengthapproximate = models.CharField(max_length=10, db_column='LengthApproximate', blank=True) # Field name made lowercase.
    charge = models.SmallIntegerField(null=True, db_column='Charge', blank=True) # Field name made lowercase.
    fragment = models.CharField(max_length=1, db_column='Fragment', blank=True) # Field name made lowercase.
    molecularweightcalc = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='MolecularWeightCalc', blank=True) # Field name made lowercase.
    molecularweightexp = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='MolecularWeightExp', blank=True) # Field name made lowercase.
    picalc = models.CharField(max_length=50, db_column='PICalc', blank=True) # Field name made lowercase.
    piexp = models.CharField(max_length=50, db_column='PIExp', blank=True) # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    function = models.ManyToManyField('Function', related_name = 'protein_function', through = 'ProteinWidfunctionWid')
    spot = models.ManyToManyField('Spot', related_name = 'protein_spot', through = 'ProteinWidspotWid')
    class Meta:
        db_table = 'Protein'

class ProteinWidfunctionWid(models.Model):
    proteinWid = models.ForeignKey(Protein, db_column='ProteinWID', related_name = '+') # Field name made lowercase.
    functionWid = models.ForeignKey(Function, db_column='FunctionWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'ProteinWIDFunctionWID'

class ProteinWidspotWid(models.Model):
    proteinWid = models.ForeignKey(Protein, db_column='ProteinWID', related_name = '+') # Field name made lowercase.
    spotWid = models.ForeignKey('Spot', db_column='SpotWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'ProteinWIDSpotWID'

class ProtocolapplWidpersonWid(models.Model):
    protocolapplicationWid = models.ForeignKey(Parameterizableapplication, db_column='ProtocolApplicationWID', related_name = '+') # Field name made lowercase.
    personWid = models.ForeignKey(Contact, db_column='PersonWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'ProtocolApplWIDPersonWID'

class ProtocolWidhardwareWid(models.Model):
    protocolWid = models.ForeignKey(Parameterizable, db_column='ProtocolWID', related_name = '+') # Field name made lowercase.
    hardwareWid = models.ForeignKey(Parameterizable, db_column='HardwareWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'ProtocolWIDHardwareWID'

class ProtocolWidsoftwareWid(models.Model):
    protocolWid = models.ForeignKey(Parameterizable, db_column='ProtocolWID', related_name = '+') # Field name made lowercase.
    softwareWid = models.ForeignKey(Parameterizable, db_column='SoftwareWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'ProtocolWIDSoftwareWID'

class QuanttymapWidquanttymapwi(models.Model):
    quantitationtypemappingWid = models.ForeignKey('Quantitationtypemapping', db_column='QuantitationTypeMappingWID', related_name = '+') # Field name made lowercase.
    quantitationtypemapWid = models.ForeignKey(Bioevent, db_column='QuantitationTypeMapWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'QuantTyMapWIDQuantTyMapWI'

class QuanttypedimensWidquanttypeWid(models.Model):
    quantitationtypedimensionWid = models.ForeignKey('Quantitationtypedimension', db_column='QuantitationTypeDimensionWID', related_name = '+') # Field name made lowercase.
    quantitationtypeWid = models.ForeignKey('Quantitationtype', db_column='QuantitationTypeWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'QuantTypeDimensWIDQuantTypeWID'

class QuanttypemapWidquanttypeWid(models.Model):
    quantitationtypemapWid = models.ForeignKey(Bioevent, db_column='QuantitationTypeMapWID', related_name = '+') # Field name made lowercase.
    quantitationtypeWid = models.ForeignKey('Quantitationtype', db_column='QuantitationTypeWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'QuantTypeMapWIDQuantTypeWID'

class QuanttypeWidconfidenceindWid(models.Model):
    quantitationtypeWid = models.ForeignKey('Quantitationtype', db_column='QuantitationTypeWID', related_name = '+') # Field name made lowercase.
    confidenceindicatorWid = models.ForeignKey('Quantitationtype', db_column='ConfidenceIndicatorWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'QuantTypeWIDConfidenceIndWID'

class QuanttypeWidquanttypemapWid(models.Model):
    quantitationtypeWid = models.ForeignKey('Quantitationtype', db_column='QuantitationTypeWID', related_name = '+') # Field name made lowercase.
    quantitationtypemapWid = models.ForeignKey(Bioevent, db_column='QuantitationTypeMapWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'QuantTypeWIDQuantTypeMapWID'

class Quantitationtype(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    mageclass = models.CharField(max_length=100, db_column='MAGEClass') # Field name made lowercase.
    identifier = models.CharField(max_length=255, db_column='Identifier', blank=True) # Field name made lowercase.
    name = models.CharField(max_length=255, db_column='Name', blank=True) # Field name made lowercase.
    isbackground = models.CharField(max_length=1, db_column='IsBackground', blank=True) # Field name made lowercase.
    channel = models.ForeignKey(Channel, null=True, db_column='Channel', blank=True, related_name = '+') # Field name made lowercase.
    quantitationtype_scale = models.ForeignKey('Term', null=True, db_column='QuantitationType_Scale', blank=True, related_name = '+') # Field name made lowercase.
    quantitationtype_datatype = models.ForeignKey('Term', null=True, db_column='QuantitationType_DataType', blank=True, related_name = '+') # Field name made lowercase.
    targetquantitationtype = models.ForeignKey('self', null=True, db_column='TargetQuantitationType', blank=True, related_name = '+') # Field name made lowercase.
    confidenceind = models.ManyToManyField('Quantitationtype', related_name = 'quanttype_confidenceind', through = 'QuanttypeWidconfidenceindWid')
    quanttypemap = models.ManyToManyField('Bioevent', related_name = 'quanttype_quanttypemap', through = 'QuanttypeWidquanttypemapWid')
    class Meta:
        db_table = 'QuantitationType'

class Quantitationtypedimension(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    identifier = models.CharField(max_length=255, db_column='Identifier', blank=True) # Field name made lowercase.
    name = models.CharField(max_length=255, db_column='Name', blank=True) # Field name made lowercase.
    quanttype = models.ManyToManyField('Quantitationtype', related_name = 'quanttypedimens_quanttype', through = 'QuanttypedimensWidquanttypeWid')
    class Meta:
        db_table = 'QuantitationTypeDimension'

class Quantitationtypemapping(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'QuantitationTypeMapping'

class Quantitationtypetuple(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    designelementtuple = models.ForeignKey(Designelementtuple, null=True, db_column='DesignElementTuple', blank=True, related_name = '+') # Field name made lowercase.
    quantitationtype = models.ForeignKey(Quantitationtype, null=True, db_column='QuantitationType', blank=True, related_name = '+') # Field name made lowercase.
    quantitationtypetuple_datum = models.ForeignKey(Datum, null=True, db_column='QuantitationTypeTuple_Datum', blank=True, related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'QuantitationTypeTuple'

class Reactant(models.Model):
    reactionWid = models.ForeignKey('Reaction', db_column='ReactionWID', related_name = '+') # Field name made lowercase.
    otherWid = models.BigIntegerField(db_column='OtherWID') # Field name made lowercase.
    coefficient = models.SmallIntegerField(db_column='Coefficient') # Field name made lowercase.
    class Meta:
        db_table = 'Reactant'

class Reaction(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    name = models.CharField(max_length=250, db_column='Name', blank=True) # Field name made lowercase.
    deltag = models.CharField(max_length=50, db_column='DeltaG', blank=True) # Field name made lowercase.
    ecnumber = models.CharField(max_length=50, db_column='ECNumber', blank=True) # Field name made lowercase.
    ecnumberproposed = models.CharField(max_length=50, db_column='ECNumberProposed', blank=True) # Field name made lowercase.
    spontaneous = models.CharField(max_length=1, db_column='Spontaneous', blank=True) # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'Reaction'

class Relatedterm(models.Model):
    termWid = models.ForeignKey('Term', db_column='TermWID', related_name = '+') # Field name made lowercase.
    otherWid = models.BigIntegerField(db_column='OtherWID') # Field name made lowercase.
    relationship = models.CharField(max_length=50, db_column='Relationship', blank=True) # Field name made lowercase.
    class Meta:
        db_table = 'RelatedTerm'

class ReporterdimensWidreporterWid(models.Model):
    reporterdimensionWid = models.ForeignKey(Designelementdimension, db_column='ReporterDimensionWID', related_name = '+') # Field name made lowercase.
    reporterWid = models.ForeignKey(Designelement, db_column='ReporterWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'ReporterDimensWIDReporterWID'

class ReportergroupWidreporterWid(models.Model):
    reportergroupWid = models.ForeignKey(Designelementgroup, db_column='ReporterGroupWID', related_name = '+') # Field name made lowercase.
    reporterWid = models.ForeignKey(Designelement, db_column='ReporterWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'ReporterGroupWIDReporterWID'

class ReporterWidbiosequenceWid(models.Model):
    reporterWid = models.ForeignKey(Designelement, db_column='ReporterWID', related_name = '+') # Field name made lowercase.
    biosequenceWid = models.BigIntegerField(db_column='BioSequenceWID') # Field name made lowercase.
    class Meta:
        db_table = 'ReporterWIDBioSequenceWID'

class ReporterWidfeaturerepormapWid(models.Model):
    reporterWid = models.ForeignKey(Designelement, db_column='ReporterWID', related_name = '+') # Field name made lowercase.
    featurereportermapWid = models.ForeignKey(Bioevent, db_column='FeatureReporterMapWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'ReporterWIDFeatureReporMapWID'

class Seqfeaturelocation(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    seqfeature_regions = models.ForeignKey(Feature, null=True, db_column='SeqFeature_Regions', blank=True, related_name = '+') # Field name made lowercase.
    strandtype = models.CharField(max_length=255, db_column='StrandType', blank=True) # Field name made lowercase.
    seqfeaturelocation_subregions = models.ForeignKey('self', null=True, db_column='SeqFeatureLocation_Subregions', blank=True, related_name = '+') # Field name made lowercase.
    seqfeaturelocation_coordinate = models.ForeignKey('Sequenceposition', null=True, db_column='SeqFeatureLocation_Coordinate', blank=True, related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'SeqFeatureLocation'

class Sequencematch(models.Model):
    queryWid = models.BigIntegerField(db_column='QueryWID') # Field name made lowercase.
    matchWid = models.BigIntegerField(db_column='MatchWID') # Field name made lowercase.
    computationWid = models.ForeignKey(Computation, db_column='ComputationWID', related_name = '+') # Field name made lowercase.
    evalue = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='EValue', blank=True) # Field name made lowercase.
    pvalue = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='PValue', blank=True) # Field name made lowercase.
    percentidentical = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='PercentIdentical', blank=True) # Field name made lowercase.
    percentsimilar = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='PercentSimilar', blank=True) # Field name made lowercase.
    rank = models.SmallIntegerField(null=True, db_column='Rank', blank=True) # Field name made lowercase.
    length = models.BigIntegerField(null=True, db_column='Length', blank=True) # Field name made lowercase.
    querystart = models.BigIntegerField(null=True, db_column='QueryStart', blank=True) # Field name made lowercase.
    queryend = models.BigIntegerField(null=True, db_column='QueryEnd', blank=True) # Field name made lowercase.
    matchstart = models.BigIntegerField(null=True, db_column='MatchStart', blank=True) # Field name made lowercase.
    matchend = models.IntegerField(null=True, db_column='MatchEnd', blank=True) # Field name made lowercase.
    class Meta:
        db_table = 'SequenceMatch'

class Sequenceposition(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    mageclass = models.CharField(max_length=100, db_column='MAGEClass') # Field name made lowercase.
    start_field = models.SmallIntegerField(null=True, db_column='Start_', blank=True) # Field name made lowercase. Field renamed because it ended with '_'.
    end = models.SmallIntegerField(null=True, db_column='End', blank=True) # Field name made lowercase.
    compositecompositemap = models.ForeignKey(Bioevent, null=True, db_column='CompositeCompositeMap', blank=True, related_name = '+') # Field name made lowercase.
    composite = models.ForeignKey(Designelement, null=True, db_column='Composite', blank=True, related_name = '+') # Field name made lowercase.
    reportercompositemap = models.ForeignKey(Bioevent, null=True, db_column='ReporterCompositeMap', blank=True, related_name = '+') # Field name made lowercase.
    reporter = models.ForeignKey(Designelement, null=True, db_column='Reporter', blank=True, related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'SequencePosition'

class SoftwareWidcontactWid(models.Model):
    softwareWid = models.ForeignKey(Parameterizable, db_column='SoftwareWID', related_name = '+') # Field name made lowercase.
    contactWid = models.ForeignKey(Contact, db_column='ContactWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'SoftwareWIDContactWID'

class SoftwareWidsoftwareWid(models.Model):
    softwareWid1 = models.ForeignKey(Parameterizable, db_column='SoftwareWID1', related_name = '+') # Field name made lowercase.
    softwareWid2 = models.ForeignKey(Parameterizable, db_column='SoftwareWID2', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'SoftwareWIDSoftwareWID'

class SpecialWidtable(models.Model):
    previousWid = models.BigIntegerField(primary_key=True, db_column='PreviousWID') # Field name made lowercase.
    class Meta:
        db_table = 'SpecialWIDTable'

class Spot(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    spotid = models.CharField(max_length=25, db_column='SpotId', blank=True) # Field name made lowercase.
    molecularweightest = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='MolecularWeightEst', blank=True) # Field name made lowercase.
    piest = models.CharField(max_length=50, db_column='PIEst', blank=True) # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DatasetWID', related_name = '+') # Field name made lowercase.
    spotidmethod = models.ManyToManyField('Spotidmethod', related_name = 'spot_spotidmethod', through = 'SpotWidspotidmethodWid')
    class Meta:
        db_table = 'Spot'

class Spotidmethod(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    methodname = models.CharField(max_length=100, db_column='MethodName') # Field name made lowercase.
    methoddesc = models.CharField(max_length=500, db_column='MethodDesc', blank=True) # Field name made lowercase.
    methodabbrev = models.CharField(max_length=10, db_column='MethodAbbrev', blank=True) # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DatasetWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'SpotIdMethod'

class SpotWidspotidmethodWid(models.Model):
    spotWid = models.ForeignKey(Spot, db_column='SpotWID', related_name = '+') # Field name made lowercase.
    spotidmethodWid = models.ForeignKey(Spotidmethod, db_column='SpotIdMethodWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'SpotWIDSpotIdMethodWID'

class Subsequence(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    nucleicacidWid = models.ForeignKey(Nucleicacid, db_column='NucleicAcidWID', related_name = '+') # Field name made lowercase.
    fullsequence = models.CharField(max_length=1, db_column='FullSequence', blank=True) # Field name made lowercase.
    sequence = models.TextField(db_column='Sequence', blank=True) # Field name made lowercase.
    length = models.BigIntegerField(null=True, db_column='Length', blank=True) # Field name made lowercase.
    lengthapproximate = models.CharField(max_length=10, db_column='LengthApproximate', blank=True) # Field name made lowercase.
    percentgc = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='PercentGC', blank=True) # Field name made lowercase.
    version = models.CharField(max_length=30, db_column='Version', blank=True) # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'Subsequence'

class Subunit(models.Model):
    complexWid = models.ForeignKey(Protein, db_column='ComplexWID', related_name = '+') # Field name made lowercase.
    subunitWid = models.ForeignKey(Protein, db_column='SubunitWID', related_name = '+') # Field name made lowercase.
    coefficient = models.SmallIntegerField(null=True, db_column='Coefficient', blank=True) # Field name made lowercase.
    class Meta:
        db_table = 'Subunit'

class Superpathway(models.Model):
    subpathwayWid = models.ForeignKey(Pathway, db_column='SubPathwayWID', related_name = '+') # Field name made lowercase.
    superpathwayWid = models.ForeignKey(Pathway, db_column='SuperPathwayWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'SuperPathway'

class Support(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    otherWid = models.BigIntegerField(db_column='OtherWID') # Field name made lowercase.
    type = models.CharField(max_length=100, db_column='Type', blank=True) # Field name made lowercase.
    evidencetype = models.CharField(max_length=100, db_column='EvidenceType', blank=True) # Field name made lowercase.
    confidence = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='Confidence', blank=True) # Field name made lowercase.
    datasetWid = models.BigIntegerField(db_column='DataSetWID') # Field name made lowercase.
    class Meta:
        db_table = 'Support'

class Synonymtable(models.Model):
    otherWid = models.BigIntegerField(db_column='OtherWID') # Field name made lowercase.
    syn = models.CharField(max_length=255, db_column='Syn') # Field name made lowercase.
    class Meta:
        db_table = 'SynonymTable'

class Taxon(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    parentWid = models.BigIntegerField(null=True, db_column='ParentWID', blank=True) # Field name made lowercase.
    name = models.CharField(max_length=100, db_column='Name', blank=True) # Field name made lowercase.
    rank = models.CharField(max_length=100, db_column='Rank', blank=True) # Field name made lowercase.
    divisionWid = models.ForeignKey(Division, null=True, db_column='DivisionWID', blank=True, related_name = '+') # Field name made lowercase.
    inheriteddivision = models.CharField(max_length=1, db_column='InheritedDivision', blank=True) # Field name made lowercase.
    gencodeWid = models.ForeignKey(Geneticcode, null=True, db_column='GencodeWID', blank=True, related_name = '+') # Field name made lowercase.
    inheritedgencode = models.CharField(max_length=1, db_column='InheritedGencode', blank=True) # Field name made lowercase.
    mcgencodeWid = models.BigIntegerField(null=True, db_column='MCGencodeWID', blank=True) # Field name made lowercase.
    inheritedmcgencode = models.CharField(max_length=1, db_column='InheritedMCGencode', blank=True) # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'Taxon'

class Term(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    name = models.CharField(max_length=255, db_column='Name') # Field name made lowercase.
    definition = models.TextField(db_column='Definition', blank=True) # Field name made lowercase.
    hierarchical = models.CharField(max_length=1, db_column='Hierarchical', blank=True) # Field name made lowercase.
    root = models.CharField(max_length=1, db_column='Root', blank=True) # Field name made lowercase.
    obsolete = models.CharField(max_length=1, db_column='Obsolete', blank=True) # Field name made lowercase.
    datasetWid = models.BigIntegerField(db_column='DataSetWID') # Field name made lowercase.
    class Meta:
        db_table = 'Term'

class Termrelationship(models.Model):
    termWid = models.ForeignKey(Term, db_column='TermWID', related_name = '+') # Field name made lowercase.
    relatedtermWid = models.ForeignKey(Term, db_column='RelatedTermWID', related_name = '+') # Field name made lowercase.
    relationship = models.CharField(max_length=10, db_column='Relationship') # Field name made lowercase.
    class Meta:
        db_table = 'TermRelationship'

class Tooladvice(models.Model):
    otherWid = models.BigIntegerField(db_column='OtherWID') # Field name made lowercase.
    toolname = models.CharField(max_length=50, db_column='ToolName') # Field name made lowercase.
    advice = models.TextField(db_column='Advice', blank=True) # Field name made lowercase.
    class Meta:
        db_table = 'ToolAdvice'

class Transcriptionunit(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    name = models.CharField(max_length=255, db_column='Name') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'TranscriptionUnit'

class Transcriptionunitcomponent(models.Model):
    type = models.CharField(max_length=100, db_column='Type') # Field name made lowercase.
    transcriptionunitWid = models.ForeignKey(Transcriptionunit, primary_key = True, db_column='TranscriptionUnitWID', related_name = '+') # Field name made lowercase.
    otherWid = models.BigIntegerField(db_column='OtherWID') # Field name made lowercase.
    class Meta:
        db_table = 'TranscriptionUnitComponent'

class TransformWidbioassaydataWid(models.Model):
    transformationWid = models.ForeignKey(Bioevent, db_column='TransformationWID', related_name = '+') # Field name made lowercase.
    bioassaydataWid = models.ForeignKey(Bioassaydata, db_column='BioAssayDataWID', related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'TransformWIDBioAssayDataWID'

class Unit(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    mageclass = models.CharField(max_length=100, db_column='MAGEClass') # Field name made lowercase.
    unitname = models.CharField(max_length=255, db_column='UnitName', blank=True) # Field name made lowercase.
    unitnamecv = models.CharField(max_length=25, db_column='UnitNameCV', blank=True) # Field name made lowercase.
    unitnamecv2 = models.CharField(max_length=25, db_column='UnitNameCV2', blank=True) # Field name made lowercase.
    unitnamecv3 = models.CharField(max_length=25, db_column='UnitNameCV3', blank=True) # Field name made lowercase.
    unitnamecv4 = models.CharField(max_length=25, db_column='UnitNameCV4', blank=True) # Field name made lowercase.
    unitnamecv5 = models.CharField(max_length=25, db_column='UnitNameCV5', blank=True) # Field name made lowercase.
    unitnamecv6 = models.CharField(max_length=25, db_column='UnitNameCV6', blank=True) # Field name made lowercase.
    unitnamecv7 = models.CharField(max_length=25, db_column='UnitNameCV7', blank=True) # Field name made lowercase.
    class Meta:
        db_table = 'Unit'

class Valence(models.Model):
    otherWid = models.ForeignKey(Element, db_column='OtherWID', related_name = '+') # Field name made lowercase.
    valence = models.SmallIntegerField(db_column='Valence') # Field name made lowercase.
    class Meta:
        db_table = 'Valence'

class Widtable(models.Model):
    previousWid = models.BigIntegerField(primary_key=True, db_column='PreviousWID') # Field name made lowercase.
    class Meta:
        db_table = 'WIDTable'

class Warehouse(models.Model):
    version = models.DecimalField(decimal_places=3, max_digits=6, db_column='Version') # Field name made lowercase.
    loaddate = models.DateTimeField(db_column='LoadDate') # Field name made lowercase.
    maxspecialWid = models.BigIntegerField(db_column='MaxSpecialWID') # Field name made lowercase.
    maxreservedWid = models.BigIntegerField(db_column='MaxReservedWID') # Field name made lowercase.
    description = models.TextField(db_column='Description', blank=True) # Field name made lowercase.
    class Meta:
        db_table = 'Warehouse'

class Zone(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    identifier = models.CharField(max_length=255, db_column='Identifier', blank=True) # Field name made lowercase.
    name = models.CharField(max_length=255, db_column='Name', blank=True) # Field name made lowercase.
    row_field = models.SmallIntegerField(null=True, db_column='Row_', blank=True) # Field name made lowercase. Field renamed because it ended with '_'.
    column_field = models.SmallIntegerField(null=True, db_column='Column_', blank=True) # Field name made lowercase. Field renamed because it ended with '_'.
    upperleftx = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='UpperLeftX', blank=True) # Field name made lowercase.
    upperlefty = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='UpperLeftY', blank=True) # Field name made lowercase.
    lowerrightx = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='LowerRightX', blank=True) # Field name made lowercase.
    lowerrighty = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='LowerRightY', blank=True) # Field name made lowercase.
    zone_distanceunit = models.ForeignKey(Unit, null=True, db_column='Zone_DistanceUnit', blank=True, related_name = '+') # Field name made lowercase.
    zonegroup_zonelocations = models.ForeignKey('Zonegroup', null=True, db_column='ZoneGroup_ZoneLocations', blank=True, related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'Zone'

class Zonedefect(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    arraymanufacturedeviation = models.ForeignKey(Arraymanufacturedeviation, null=True, db_column='ArrayManufactureDeviation', blank=True, related_name = '+') # Field name made lowercase.
    zonedefect_defecttype = models.ForeignKey(Term, null=True, db_column='ZoneDefect_DefectType', blank=True, related_name = '+') # Field name made lowercase.
    zonedefect_positiondelta = models.ForeignKey(Positiondelta, null=True, db_column='ZoneDefect_PositionDelta', blank=True, related_name = '+') # Field name made lowercase.
    zone = models.ForeignKey(Zone, null=True, db_column='Zone', blank=True, related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'ZoneDefect'

class Zonegroup(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    physicalarraydesign_zonegroups = models.ForeignKey(Arraydesign, null=True, db_column='PhysicalArrayDesign_ZoneGroups', blank=True, related_name = '+') # Field name made lowercase.
    spacingsbetweenzonesx = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='SpacingsBetweenZonesX', blank=True) # Field name made lowercase.
    spacingsbetweenzonesy = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='SpacingsBetweenZonesY', blank=True) # Field name made lowercase.
    zonesperx = models.SmallIntegerField(null=True, db_column='ZonesPerX', blank=True) # Field name made lowercase.
    zonespery = models.SmallIntegerField(null=True, db_column='ZonesPerY', blank=True) # Field name made lowercase.
    zonegroup_distanceunit = models.ForeignKey(Unit, null=True, db_column='ZoneGroup_DistanceUnit', blank=True, related_name = '+') # Field name made lowercase.
    zonegroup_zonelayout = models.ForeignKey('Zonelayout', null=True, db_column='ZoneGroup_ZoneLayout', blank=True, related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'ZoneGroup'

class Zonelayout(models.Model):
    Wid = models.BigIntegerField(primary_key=True, db_column='WID') # Field name made lowercase.
    datasetWid = models.ForeignKey(Dataset, db_column='DataSetWID', related_name = '+') # Field name made lowercase.
    numfeaturesperrow = models.SmallIntegerField(null=True, db_column='NumFeaturesPerRow', blank=True) # Field name made lowercase.
    numfeaturespercol = models.SmallIntegerField(null=True, db_column='NumFeaturesPerCol', blank=True) # Field name made lowercase.
    spacingbetweenrows = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='SpacingBetweenRows', blank=True) # Field name made lowercase.
    spacingbetweencols = models.DecimalField(decimal_places=65535, null=True, max_digits=65535, db_column='SpacingBetweenCols', blank=True) # Field name made lowercase.
    zonelayout_distanceunit = models.ForeignKey(Unit, null=True, db_column='ZoneLayout_DistanceUnit', blank=True, related_name = '+') # Field name made lowercase.
    class Meta:
        db_table = 'ZoneLayout'

