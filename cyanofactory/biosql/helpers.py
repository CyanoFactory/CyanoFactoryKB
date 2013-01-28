import settings
from BioSQL import BioSeqDatabase
from cyano.cache import Cache
import Bio.Seq

database_name = "cyanofactory"

class Res:
	def __init__(self, dbrecord):
		dbrecord.features
		for x in dbrecord.features:
			for y in x.qualifiers.iterkeys():
				y
				x.qualifiers[y]
		
		self.features = dbrecord.features
		self.seq = Bio.Seq.Seq(str(dbrecord.seq))
		
def get_database_handle():
	dbdata = settings.DATABASES['default']
	server = BioSeqDatabase.open_database(driver="psycopg2", user=dbdata['USER'],
                     passwd = dbdata['PASSWORD'], host = dbdata['HOST'], db=dbdata['NAME'])
	return server

def get_database_item(organism):
	handle = get_database_handle()
	db = handle[database_name]
	dbrecord = db.lookup(name = organism.name)

	record = Cache.try_get(organism.get_absolute_url(),
						   lambda : __retrieve_item(dbrecord))
	return record

def __retrieve_item(dbrecord):
	#dbrecord.features
	#Bio.Seq.Seq(str(dbrecord.seq))
	#for x in dbrecord.features:
#		for y in x.qualifiers.iterkeys():
#			y
#			x.qualifiers[y]
	return Res(dbrecord)
