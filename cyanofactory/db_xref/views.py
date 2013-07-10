from django.http import HttpResponseBadRequest
from django.shortcuts import render_to_response
from django.shortcuts import redirect
from django.template import Context, loader
from db_xref.helpers import get_database_url

def index(request):
    """Displayed when no arguments are passed. Renders the intro page.

    :return: Website -- ``db_xref/index.html``
    """
    data = {"db" : "ECOCYC",
            "value": "G7954",
            "item" : "ECOCYC:G7954",
            "formats" : ["txt", "html", "xml", "json", "redirect"]
            }
    
    return render_to_response("db_xref/index.html", data)

def dbxref(request, source, xid):
    """Converts the db_xref identifier to a database url.

    :URL:
        ``(?P<source>[a-zA-Z0-9_-]+)\s*:\s*(?P<xid>[a-zA-Z0-9_.-]+)``

    :Example:
        ``ECOCYC:G7954``
    
    :GET:
        ``format:``
            * ``redirect``: Default value. Redirects to the database.
            * ``txt``: Returns a text file containing the url.
            * ``json``: Returns a json file.
            * ``xml``: Returns a XML file.
            * ``html``: Returns a html file.
    
    :raises:
        ``HttpResponseBadRequest:`` -- Unsupported database or format.
      
    :return:
       Website -- ``db_xref/output.format``
    """
    _format = request.GET.get('format', 'redirect')
    
    database = source
    organism = xid
    
    data = {
        "database" : database,
        "item" : organism,
        "url" : get_database_url(database, organism)()
    }
    print data["url"]
    
    template = "db_xref/output." + _format

    do_error = False

    if _format == "txt":
        mimetype = "text/plain"
    elif _format == "json":
        mimetype = "application/json"
    elif _format == "xml":
        mimetype = "application/xml"
    elif _format == "html":
        mimetype = "application/xhtml+xml"
    elif _format == "redirect":
        if len(data["url"]) == 0:
            do_error = True
            data["error"] = "Unsupported Database " + database
        else:
            return redirect(data["url"])
    else:
        do_error = True
        data["error"] = "Unknown format " + _format
        
    if do_error:
        t = loader.get_template('db_xref/error.html')
        c = Context(data)
        return HttpResponseBadRequest(
            t.render(c),
            mimetype = 'text/html; charset=UTF-8',
            content_type = 'text/html; charset=UTF-8')

    return render_to_response(template, data, mimetype = mimetype)
