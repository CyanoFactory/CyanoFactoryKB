"""
Copyright (c) 2014 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

import StringIO
from collections import OrderedDict
import os
import tempfile
from crispy_forms.utils import render_crispy_form
from django.core.exceptions import ObjectDoesNotExist
from django.http.response import HttpResponse, HttpResponseBadRequest
from django.template.context import RequestContext
from django.views.decorators.csrf import ensure_csrf_cookie
from jsonview.decorators import json_view
from cyano.decorators import ajax_required, permission_required
from PyNetMet2.metabolism import Metabolism
from cyano.helpers import render_queryset_to_response, render_queryset_to_response_error
from cyanodesign.forms import UploadModelForm
from cyanodesign.helpers import model_from_string, apply_commandlist, calc_reactions, get_selected_reaction
from .models import DesignModel, Revision
import networkx as nx
from networkx.readwrite import json_graph
import pygraphviz
import json


@permission_required("access_cyanodesign")
def index(request):
    upload_form = UploadModelForm(None)

    models = DesignModel.objects.filter(user=request.user.profile)

    return render_queryset_to_response(
        request,
        template="cyanodesign/list.html",
        queryset=models,
        data={'upload_form': upload_form}
    )


@ensure_csrf_cookie
@permission_required("access_cyanodesign")
def design(request, pk):
    try:
        item = DesignModel.objects.get(user=request.user.profile, pk=pk)
        current = item.get_latest_revision()
    except ObjectDoesNotExist:
        return render_queryset_to_response_error(request, error=404, msg="Model not found")

    try:
        revision = request.GET["revision"]
        try:
            revision = Revision.objects.get(model=item, pk=revision)
        except ObjectDoesNotExist:
            return render_queryset_to_response_error(request, error=404, msg="Revision not found")
    except KeyError:
        revision = current

    data = {"pk": pk, "name": item.name, "revision": None if current.pk == revision.pk else revision}

    return render_queryset_to_response(request, template="cyanodesign/design.html", data=data)


@ajax_required
@permission_required("access_cyanodesign")
def get_reactions(request, pk):
    try:
        item = DesignModel.objects.get(user=request.user.profile, pk=pk)
    except ObjectDoesNotExist:
        return HttpResponseBadRequest("Bad Model")

    try:
        revision = request.GET["revision"]
        try:
            revision = Revision.objects.get(model=item, pk=revision)
        except ObjectDoesNotExist:
            return HttpResponseBadRequest("Bad Revision")
    except KeyError:
        revision = item.get_latest_revision()

    org = model_from_string(revision.content)

    if org is None:
        return HttpResponseBadRequest("Bad Model")

    ret = []

    for enzyme in org.enzymes:
        # Filter out auto-created external transport funcs
        if not enzyme.name.endswith("_transp"):
            ret.append({
                "name": enzyme.name,
                "stoichiometric": enzyme.stoic,
                "substrates": enzyme.substrates,
                "products": enzyme.products,
                "reversible": enzyme.reversible,
                "constraints": enzyme.constraint,
                "pathway": enzyme.pathway,
                "disabled": enzyme.disabled
            })

    outgraph = ""

    return HttpResponse(json.dumps({
        "external": filter(lambda x: not x.startswith("#"), org.external),
        "enzymes": ret,
        "objective": org.obj[0] if org.obj else None,
        "graph": outgraph}), content_type="application/json")


@permission_required("access_cyanodesign")
def simulate(request, pk):
    if not all(x in request.POST for x in ["changes", "objectives", "display", "auto_flux"]):
        return HttpResponseBadRequest("Request incomplete")

    output_format = request.POST.get("format", "json")
    dry_run = request.POST.get("dry_run", False)

    try:
        model = DesignModel.objects.get(user=request.user.profile, pk=pk)
    except ObjectDoesNotExist:
        return HttpResponseBadRequest("Bad Model")

    try:
        revision = request.POST["revision"]
        try:
            revision = Revision.objects.get(model=model, pk=revision)
        except ObjectDoesNotExist:
            return HttpResponseBadRequest("Bad Revision")
    except KeyError:
        revision = model.get_latest_revision()

    org = model_from_string(revision.content)

    try:
        changes = json.loads(request.POST["changes"])
        objectives = json.loads(request.POST["objectives"])
        display = json.loads(request.POST["display"])
        auto_flux = json.loads(request.POST["auto_flux"])
    except ValueError:
        return HttpResponseBadRequest("Invalid JSON data")

    try:
        auto_flux = bool(auto_flux)
    except ValueError:
        return HttpResponseBadRequest("Invalid data type")

    try:
        org = apply_commandlist(org, changes)

        if not objectives:
            raise ValueError("No objective specified")
        org.obj = []
        for obj in objectives:
            obj_reac = org.get_reaction(obj["name"])
            if obj_reac is None:
                raise ValueError("Objective not in model: " + obj["name"])
            org.obj.append([obj_reac.name, "1" if obj["maximize"] else -1])

            if obj_reac.disabled:
                return HttpResponseBadRequest("Objective disabled: " + obj_reac.name)
    except ValueError as e:
        return HttpResponseBadRequest("Model error: " + e.message)

    for reaction in org.enzymes:
        reaction.disabled = False

    display = filter(lambda x: not len(x) == 0, display)

    for item in display:
        if not org.has_reaction(item):
            return HttpResponseBadRequest("Unknown reaction in display list: " + item)

    try:
        fba = org.fba()
    except ValueError as e:
        return HttpResponseBadRequest("FBA error: " + e.message)

    if dry_run:
        return HttpResponse(
            json.dumps({"solution": fba.get_status()}),
            content_type="application/json")

    full_g, nodeIDs = calc_reactions(org, fba)

    display = json.loads(request.POST["display"])
    display = filter(lambda x: len(x) > 0, display)

    dflux = {}
    for reac, flux in zip(map(lambda x: x.name, fba.reacs), fba.flux):
        dflux[reac] = flux

    # Auto filter by FBA
    if auto_flux:
        if len(full_g.edges()) <= 30:
            full_eg = full_g
        else:
            full_eg = nx.ego_graph(full_g.reverse(), nodeIDs[org.obj[0][0]], radius=3, center=False, undirected=False)

        full_g.remove_edges_from(full_g.in_edges(nodeIDs[org.obj[0][0]]) + full_g.out_edges(nodeIDs[org.obj[0][0]]))
        all_edges = map(lambda x: full_eg.get_edge_data(*x)["object"].name, full_eg.edges())
        # Get fluxes of "edges"
        flux = []
        for reac in all_edges:
            flux.append([reac, dflux[reac]])
        flux = sorted(flux, key=lambda x: -x[1])
        display = map(lambda x: x[0], flux[:30])
    else:
        full_g.remove_edges_from(full_g.in_edges(nodeIDs[org.obj[0][0]]) + full_g.out_edges(nodeIDs[org.obj[0][0]]))

    display.append(org.obj[0][0])

    g = get_selected_reaction(json_graph.node_link_data(full_g), nodeIDs, display, org)

    graph = nx.to_agraph(g)
    graph.graph_attr.update(splines=True, overlap=False, rankdir="LR")
    graph.node_attr.update(style="filled", colorscheme="pastel19")
    outgraph = str(graph.to_string())
    outgraph = pygraphviz.AGraph(outgraph)
    outgraph = outgraph.draw(format="svg", prog="dot")

    if output_format == "json":
        return HttpResponse(json.dumps(
            {"graph": outgraph,
            "solution": fba.get_status(),
            "flux": dflux
            }),
            content_type="application/json"
        )
    elif output_format == "png":
        import wand.image
        with wand.image.Image(blob=outgraph, format="svg") as image:
            png_image = image.make_blob("png")

        r = HttpResponse(png_image, content_type="image/png")
        r['Content-Disposition'] = 'attachment; filename={}.png'.format(model.filename)
        return r
    elif output_format == "svg":
        r = HttpResponse(outgraph, content_type="image/svg+xml")
        r['Content-Disposition'] = 'attachment; filename={}.svg'.format(model.filename)
        return r
    elif output_format == "csv":
        s = StringIO.StringIO()
        for reac, flux in zip(fba.reacs, fba.flux):
            s.write(reac.name)
            s.write("\t")
            s.write(flux)
            s.write("\r\n")
        r = HttpResponse(s.getvalue(), content_type="text/csv")
        r['Content-Disposition'] = 'attachment; filename={}.csv'.format(model.filename)
        return r
    else:
        return HttpResponseBadRequest("Unknown format")


@permission_required("access_cyanodesign")
def export(request, pk):
    try:
        model = DesignModel.objects.get(user=request.user.profile, pk=pk)
        content = model.get_latest_revision().content
    except ObjectDoesNotExist:
        return HttpResponseBadRequest("Bad Model")

    response = HttpResponse(
        content,
        mimetype="application/x-bioopt; charset=UTF-8",
        content_type="application/x-bioopt; charset=UTF-8"
    )
    response['Content-Disposition'] = "attachment; filename=" + model.filename

    return response


@ajax_required
@permission_required("access_cyanodesign")
def save(request, pk):
    if not all(x in request.POST for x in ["changes", "objectives", "favourites"]):
        return HttpResponseBadRequest("Request incomplete")

    try:
        model = DesignModel.objects.get(user=request.user.profile, pk=pk)
    except ObjectDoesNotExist:
        return HttpResponseBadRequest("Bad Model")

    try:
        revision = request.POST["revision"]
        try:
            revision = Revision.objects.get(model=model, pk=revision)
        except ObjectDoesNotExist:
            return HttpResponseBadRequest("Bad Revision")
    except KeyError:
        revision = model.get_latest_revision()

    org = model_from_string(revision.content)

    try:
        changes = json.loads(request.POST["changes"])
        objectives = json.loads(request.POST["objectives"])
        favourites = json.loads(request.POST["favourites"])
        summary = json.loads(request.POST.get("summary", ""))
    except ValueError:
        return HttpResponseBadRequest("Invalid JSON data")

    try:
        org = apply_commandlist(org, changes)
        org.obj = []
        for obj in objectives:
            if obj["name"]:
                obj_reac = org.get_reaction(obj["name"])
                if obj_reac is None:
                    raise ValueError("Objective not in model: " + obj["name"])
                org.obj.append([obj_reac.name, "1" if obj["maximize"] else -1])
        else:
            org.obj = None
    except ValueError as e:
        return HttpResponseBadRequest("Model error: " + e.message)

    with tempfile.NamedTemporaryFile(delete=False) as fid:
        path = fid.name
        org.dump(path)

    with open(path) as f:
        model.content = f.read()
        model.save()

    os.remove(path)

    Revision(
        model=model,
        content=model.content,
        changes=dict(changes=changes,objectives=objectives,favourite=favourites),
        reason=summary
    ).save()

    return HttpResponse("OK")


@permission_required("access_cyanodesign")
@json_view
def upload(request):
    data = {}

    if request.method == 'POST':
        form = UploadModelForm(request.POST, request.FILES)

        if form.is_valid():
            name = form.cleaned_data.get('name')

            #save to temporary file
            filename = request.FILES['file'].name

            ss = StringIO.StringIO()

            with tempfile.NamedTemporaryFile(delete=False) as fid:
                path = fid.name

                for chunk in request.FILES['file'].chunks():
                    ss.write(chunk)
                    fid.write(chunk)

            try:
                Metabolism(path)
            except:
                return HttpResponseBadRequest("Bad Model")
            finally:
                os.remove(path)

            DesignModel.objects.create(
                user=request.user.profile,
                name=name,
                filename=filename,
                content=ss.getvalue()
            )

            return {'success': True}
        else:
            form_html = render_crispy_form(form, context=RequestContext(request))
            return {'success': False, 'form_html': form_html}

    return HttpResponseBadRequest()


@permission_required("access_cyanodesign")
@ajax_required
def delete(request):
    pk = request.POST.get("id", 0)

    try:
        DesignModel.objects.get(user=request.user.profile, pk=pk).delete()
    except ObjectDoesNotExist:
        return HttpResponseBadRequest("Bad Model")

    return HttpResponse("ok")


@permission_required("access_cyanodesign")
def history(request, pk):
    try:
        model = DesignModel.objects.get(user=request.user.profile, pk=pk)
        revisions = model.revisions.all()
    except ObjectDoesNotExist:
        return render_queryset_to_response_error(request, error=404, msg="Model not found")

    entries = []

    from itertools import groupby
    for k,v in groupby(revisions, key=lambda x: x.date.date()):
        entries.append([k, list(v)[::-1]])

    data = {"model": model, "revisions": entries}

    return render_queryset_to_response(request, template="cyanodesign/history.html", data=data)
