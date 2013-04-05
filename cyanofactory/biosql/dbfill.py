# Utility class to load a Genebank File into the database
from Bio import Entrez
from Bio import SeqIO
from BioSQL import BioSeqDatabase
server = BioSeqDatabase.open_database(driver="psycopg2", user="cyano",
                     passwd = "cyano", host = "localhost", db="cyano")
db = server.new_database("cyanofactory", authority="HSMW", description="cyano-Organismen")
server.commit()
db = server["cyanofactory"]
Entrez.email = "gk@trash-mail.com"
handle = Entrez.efetch(db="nuccore", id="NC_017277.1", rettype="gb", retmode="text")
count = db.load(SeqIO.parse(handle, "genbank"))
print "Loaded %i records" % count
server.commit()
