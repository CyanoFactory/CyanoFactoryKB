from cyano.helpers import render_queryset_to_response

class Importer(object):
    def __init__(self):
        self.messages = []
    
    def load(self, filename):
        pass
    
    def preview(self, request, species):        
        return render_queryset_to_response(
            species = species,
            request = request,
            template = 'cyano/preview.html', 
            data = {
                })
    
    def submit(self, request, species):
        pass

class ParserError:
    def __init__(self, expr, line, msg):
        self.expr = expr
        self.msg = msg
        self.line = line
        
    def __str__(self):
        return "L{}: {} ({})".format(self.line, self.msg, self.expr)
