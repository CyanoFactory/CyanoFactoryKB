from django.core.urlresolvers import reverse

def get_database_url_from_organism(organism):
    """Takes an db_xref identifier and converts it into an URL.
    A db_xref identifier consists of two strings delimited by a **:**.

    :param organism: db_xref identifier (**Must** contain a **:**)
    :type organism: str

    :returns: str -- Target URL for the given db_xref identifier.
    """
    stri = organism.split(':', 1)
    return reverse("db_xref.views.dbxref", args=[stri[0], stri[1]])

def get_database_url(database, organism):
    """Takes an db_xref identifier (already split in database and organism part)
    and converts it into an URL.

    :param database: Part on left side of **:**
    :type database: str
    :param organism: Part on right side of **:**
    :type organism: str

    :returns: str -- Target URL for the given db_xref identifier.
    """  
    database = database.lower()
    
    return {
        "asap": lambda : "https://asap.ahabs.wisc.edu/annotation/php/feature_info.php?FeatureID=" + organism,
        "ec": lambda : "http://enzyme.expasy.org/EC/" + organism,
        "ecocyc": lambda : "http://biocyc.org/ECOLI/new-image?type=GENE&object=" + organism,
        "ecogene": lambda : "http://ecogene.org/geneInfo.php?eg_id=" + organism,
        "geneid": lambda : "http://www.ncbi.nlm.nih.gov/sites/entrez?db=gene&cmd=Retrieve&dopt=full_report&list_uids=" + organism,
        "gi": lambda: "http://www.ncbi.nlm.nih.gov/nuccore/" + organism,
        "project": lambda: "http://www.ncbi.nlm.nih.gov/bioproject/" + organism,
        "bioproject": lambda: "http://www.ncbi.nlm.nih.gov/bioproject?term=" + organism,
        "pubmed": lambda: "http://www.ncbi.nlm.nih.gov/pubmed/" + organism,
        "refseq": lambda: "http://www.ncbi.nlm.nih.gov/nuccore/" + organism,
        "uniprotkb/swiss-prot" : lambda : "http://www.uniprot.org/uniprot/" + organism,        
            }.get(database, lambda : "")
