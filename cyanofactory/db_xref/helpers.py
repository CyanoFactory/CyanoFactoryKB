
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
    
    result = {
        "asap": lambda: "https://asap.ahabs.wisc.edu/annotation/php/feature_info.php?FeatureID=%s",
        "ec": lambda: "http://enzyme.expasy.org/EC/%s",
        "ecocyc": lambda: "http://biocyc.org/ECOLI/new-image?type=GENE&object=%s",
        "ecogene": lambda: "http://ecogene.org/geneInfo.php?eg_id=%s",
        "geneid": lambda: "https://www.ncbi.nlm.nih.gov/sites/entrez?db=gene&cmd=Retrieve&dopt=full_report&list_uids=%s",
        "gi": lambda: "https://www.ncbi.nlm.nih.gov/nuccore/%s",
        "project": lambda: "https://www.ncbi.nlm.nih.gov/bioproject/%s",
        "bioproject": lambda: "https://www.ncbi.nlm.nih.gov/bioproject?term=%s",
        "pubmed": lambda: "https://www.ncbi.nlm.nih.gov/pubmed/%s",
        "rebase": lambda: "http://rebase.neb.com/rebase/enz/%s.html",
        "refseq": lambda: "https://www.ncbi.nlm.nih.gov/nuccore/%s",
        "taxon": lambda: "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=%s",
        "uniprotkb/swiss-prot": lambda: "http://www.uniprot.org/uniprot/%s",
             }.get(database, lambda: None)()

    if result is not None:
        return result % organism

    return ""
