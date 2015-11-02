# via https://stackoverflow.com/questions/1960516
import json

class DecimalEncoder(json.JSONEncoder):
    def default(self, o):
        import decimal
        if isinstance(o, decimal.Decimal):
            return float(o)
        return super(DecimalEncoder, self).default(o)