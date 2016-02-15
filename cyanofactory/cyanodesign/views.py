"""
Copyright (c) 2014 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from io import StringIO
import os
import tempfile
from crispy_forms.utils import render_crispy_form
from django.contrib.auth.decorators import login_required
from django.core.exceptions import ObjectDoesNotExist
from django.core.urlresolvers import reverse
from django.db.transaction import atomic
from django.http.response import HttpResponse, HttpResponseBadRequest
from django.template.context import RequestContext
from django.views.decorators.csrf import ensure_csrf_cookie
from jsonview.decorators import json_view
from cyano.decorators import ajax_required, permission_required, global_permission_required
from PyNetMet2.metabolism import Metabolism
from cyano.helpers import render_queryset_to_response, render_queryset_to_response_error
from cyano.models import UserProfile
from cyanodesign.forms import UploadModelForm, ModelFromTemplateForm, SaveModelAsForm, SaveModelForm
from cyanodesign.helpers import model_from_string, apply_commandlist, calc_reactions, get_selected_reaction, \
    compress_command_list
from cyanodesign.json_model import JsonModel
from .models import DesignModel, Revision, DesignTemplate
import networkx as nx
from networkx.readwrite import json_graph
import pygraphviz
import json
from jsonview.exceptions import BadRequest

@global_permission_required("access_cyanodesign")
def index(request):
    upload_form = UploadModelForm(None)

    templates = DesignTemplate.objects.values_list("pk", "name")
    template_form = ModelFromTemplateForm(choices=templates)

    models = DesignModel.objects.filter(user=UserProfile.get_profile(request.user))

    return render_queryset_to_response(
        request,
        template="cyanodesign/list.html",
        queryset=models,
        data={'upload_form': upload_form, "template_form": template_form}
    )


@ensure_csrf_cookie
@global_permission_required("access_cyanodesign")
def design(request, pk):
    try:
        item = DesignModel.objects.get(user=UserProfile.get_profile(request.user), pk=pk)
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

    data = {
        "pk": pk,
        "name": item.name,
        "revision": None if current.pk == revision.pk else revision,
        "save_form": SaveModelForm(None),
        "saveas_form": SaveModelAsForm(None)
    }

    return render_queryset_to_response(request, template="cyanodesign/design.html", data=data)

@ajax_required
@global_permission_required("access_cyanodesign")
@json_view
def get_reactions(request, pk):
    try:
        item = DesignModel.objects.get(user=UserProfile.get_profile(request.user), pk=pk)
    except ObjectDoesNotExist:
        return BadRequest("Bad Model")

    try:
        revision = request.GET["revision"]
        try:
            revision = Revision.objects.get(model=item, pk=revision)
        except ObjectDoesNotExist:
            return BadRequest("Bad Revision")
    except KeyError:
        revision = item.get_latest_revision()

    return revision.content


@global_permission_required("access_cyanodesign")
def simulate(request, pk):
    if not all(x in request.POST for x in ["changes", "objectives", "display", "auto_flux"]):
        return HttpResponseBadRequest("Request incomplete")

    output_format = request.POST.get("format", "json")
    dry_run = request.POST.get("dry_run", False)

    try:
        model = DesignModel.objects.get(user=UserProfile.get_profile(request.user), pk=pk)
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
        apply_commandlist(org, changes)
        org = model_from_string(revision.content)
        changes = compress_command_list(changes)
        apply_commandlist(org, changes)

        if not objectives:
            raise ValueError("No objective specified")
        org.objectives = []
        for obj in objectives:
            if len(obj["name"]) > 0:
                obj_reac = org.get_reaction(obj["name"])
                if obj_reac is None:
                    raise ValueError("Objective not in model: " + obj["name"])

                if obj_reac.disabled:
                    return HttpResponseBadRequest("Objective disabled: " + obj_reac.name)

                org.objectives.append(JsonModel.Objective(**obj))
            else:
                raise ValueError("No objective specified")
    except ValueError as e:
        return HttpResponseBadRequest("Model error: " + str(e))

    display = filter(lambda x: not len(x) == 0, display)

    for item in display:
        if not org.has_reaction(item):
            return HttpResponseBadRequest("Unknown reaction in display list: " + item)

    org = org.to_model()

    try:
       fba = org.fba()
    except ValueError as e:
        return HttpResponseBadRequest("FBA error: " + str(e))

    if dry_run:
        return HttpResponse(
            json.dumps({"solution": fba.get_status()}),
            content_type="application/json")

    full_g, nodeIDs = calc_reactions(org, fba)

    display = json.loads(request.POST["display"])
    display = list(filter(lambda x: len(x) > 0 and org.has_reaction(x), display))

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
        display = list(map(lambda x: x[0], flux[:30]))
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
            {"graph": outgraph.decode("utf-8"),
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
        r = HttpResponse(outgraph.decode("utf-8"), content_type="image/svg+xml")
        r['Content-Disposition'] = 'attachment; filename={}.svg'.format(model.filename)
        return r
    elif output_format == "csv":
        s = StringIO()
        for reac, flux in zip(fba.reacs, fba.flux):
            s.write(reac.name)
            s.write("\t")
            s.write(str(flux))
            s.write("\r\n")
        r = HttpResponse(s.getvalue(), content_type="text/csv")
        r['Content-Disposition'] = 'attachment; filename={}.csv'.format(model.filename)
        return r
    else:
        return HttpResponseBadRequest("Unknown format")


@global_permission_required("access_cyanodesign")
def export(request, pk):
    form = request.GET.get("format", "bioopt")

    if form not in ["bioopt", "sbml"]:
        return HttpResponseBadRequest("Bad format")

    try:
        model = DesignModel.objects.get(user=UserProfile.get_profile(request.user), pk=pk)
        content = model.get_latest_revision().content
    except ObjectDoesNotExist:
        return HttpResponseBadRequest("Bad Model")

    org = model_from_string(content)

    ss = StringIO()

    org.to_model().dump(fileout=ss, filetype="opt" if form == "bioopt" else "sbml")

    response = HttpResponse(
        ss.getvalue(),
        content_type="application/x-bioopt" if form=="bioopt" else "application/sbml+xml"
    )

    ss.close()

    response['Content-Disposition'] = "attachment; filename=" + model.filename

    return response


@login_required
@ajax_required
@global_permission_required("access_cyanodesign")
@json_view
def save(request, pk):
    if not all(x in request.POST for x in ["changes", "objectives"]):
        return BadRequest("Request incomplete")

    try:
        model = DesignModel.objects.get(user=UserProfile.get_profile(request.user), pk=pk)
    except ObjectDoesNotExist:
        return BadRequest("Bad Model")

    try:
        revision = request.POST["revision"]
        try:
            revision = Revision.objects.get(model=model, pk=revision)
        except ObjectDoesNotExist:
            return BadRequest("Bad Revision")
    except KeyError:
        revision = model.get_latest_revision()

    org = model_from_string(revision.content)

    try:
        changes = json.loads(request.POST["changes"])
        objectives = json.loads(request.POST["objectives"])
    except ValueError:
        return BadRequest("Invalid JSON data")

    summary = request.POST.get("summary")

    try:
        apply_commandlist(org, changes)
        org = model_from_string(revision.content)
        changes = compress_command_list(changes)
        apply_commandlist(org, changes)

        if objectives:
            org.objectives = []
            for obj in objectives:
                if len(obj["name"]) > 0:
                    obj_reac = org.get_reaction(obj["name"])
                    if obj_reac is None:
                        raise ValueError("Objective not in model: " + obj["name"])

                    org.objectives.append(JsonModel.Objective(**obj))
    except ValueError as e:
        return BadRequest("Model error: " + e.message)

    Revision(
        model=model,
        content=org.to_json(),
        changes=dict(changes=changes,objectives=objectives),
        reason=summary
    ).save()

    return {}


@login_required
@global_permission_required("access_cyanodesign")
@json_view
def save_as(request, pk):
    if request.method == 'POST':
        form = SaveModelAsForm(request.POST, request.FILES)

        if form.is_valid():
            name = form.cleaned_data.get('saveas_name')
            summary = form.cleaned_data.get('saveas_summary')

            request.POST = request.POST.copy()
            request.POST["summary"] = summary
            save(request, pk)
    
            model = DesignModel.objects.get(user=UserProfile.get_profile(request.user), pk=pk)

            dm = DesignModel.objects.create(
                user=request.user.profile,
                name=name,
                filename=model.filename,
                content=model.content
            )

            # Repoint to new model
            rev = model.get_latest_revision()
            rev.model = dm
            rev.save()

            return {'success': True, 'url': reverse("cyano-design-design", kwargs={"pk":dm.pk})}
        else:
            form_html = render_crispy_form(form, context=RequestContext(request))
            return {'success': False, 'form_html': form_html}

    return BadRequest()

@login_required
@global_permission_required("access_cyanodesign")
@json_view
@atomic
def upload(request, pk):
    data = {}

    if request.method == 'POST':
        if pk == "1":
            # Uploaded model
            form = UploadModelForm(request.POST, request.FILES)

            if form.is_valid():
                name = form.cleaned_data.get('name')

                #save to temporary file
                freq = request.FILES['file']
                filename = freq.name

                ss = StringIO()
                for chunk in freq.chunks():
                    ss.write(chunk.decode("utf-8"))

                try:
                    model = Metabolism(ss)
                except:
                    return HttpResponseBadRequest("Bad Model")

                dm = DesignModel.objects.create(
                    user=request.user.profile,
                    name=name,
                    filename=filename,
                    content=ss.getvalue()
                )

                try:
                    jm = JsonModel.from_model(model).to_json()
                except ValueError:
                    return BadRequest(str(ValueError))

                Revision(
                    model=dm,
                    content=jm,
                    reason="Initial version"
                ).save()

                return {'success': True}
            else:
                form_html = render_crispy_form(form, context=request)
                return {'success': False, 'form_html': form_html}
        if pk == "2":
            # from template
            templates = DesignTemplate.objects.values_list("pk", "name")
            form = ModelFromTemplateForm(templates, request.POST, request.FILES)

            if form.is_valid():
                name = form.cleaned_data.get('name')
                choice = form.cleaned_data.get('choice')

                template = DesignTemplate.objects.get(pk=choice)

                dm = DesignModel.objects.create(
                    user=UserProfile.get_profile(request.user),
                    name=name,
                    filename=template.filename,
                    content=""
                )

                Revision(
                    model=dm,
                    content=template.content,
                    reason="Initial version"
                ).save()

                return {'success': True}
            else:
                form_html = render_crispy_form(form, context=request)
                return {'success': False, 'form_html': form_html}

    return HttpResponseBadRequest()


@login_required
@global_permission_required("access_cyanodesign")
@ajax_required
@json_view
def delete(request):
    pk = request.POST.get("id", 0)

    try:
        DesignModel.objects.get(user=UserProfile.get_profile(request.user), pk=pk).delete()
    except ObjectDoesNotExist:
        return BadRequest("Bad Model")

    return {}


@global_permission_required("access_cyanodesign")
def history(request, pk):
    try:
        model = DesignModel.objects.get(user=UserProfile.get_profile(request.user), pk=pk)
        revisions = model.revisions.all()
    except ObjectDoesNotExist:
        return render_queryset_to_response_error(request, error=404, msg="Model not found")

    entries = []

    from itertools import groupby
    for k,v in groupby(revisions, key=lambda x: x.date.date()):
        v = list(v)[::-1]
        for vv in v:
            try:
                vv.changes = json.dumps(vv.changes["changes"])
            except KeyError:
                vv.changes = json.dumps({})

        entries.append([k, list(v)[::-1]])

    data = {"model": model, "revisions": entries}

    return render_queryset_to_response(request, template="cyanodesign/history.html", data=data)
