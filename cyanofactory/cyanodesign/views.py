"""
Copyright (c) 2014 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from io import StringIO
import os
import tempfile
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
from cyano.helpers import render_queryset_to_response, render_queryset_to_response_error, render_crispy_form
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
    if not all(x in request.POST for x in ["changes", "objectives", "design_objectives", "target_reactions", "display", "auto_flux", "type"]):
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

    org_json = model_from_string(revision.content)

    try:
        changes = json.loads(request.POST["changes"])
        objectives = json.loads(request.POST["objectives"])
        design_objectives = json.loads(request.POST["design_objectives"])
        target_reactions = json.loads(request.POST["target_reactions"])
        display = json.loads(request.POST["display"])
        auto_flux = json.loads(request.POST["auto_flux"])
        simtype = json.loads(request.POST["type"])
    except ValueError:
        return HttpResponseBadRequest("Invalid JSON data")

    try:
        auto_flux = bool(auto_flux)
    except ValueError:
        return HttpResponseBadRequest("Invalid data type")

    if str(simtype) not in ["fba", "mba", "sa"]:
        return HttpResponseBadRequest("Unsupported simulation type: {}".format(simtype))

    try:
        apply_commandlist(org_json, changes)
        org_json = model_from_string(revision.content)
        changes = compress_command_list(changes)
        apply_commandlist(org_json, changes)

        if not objectives:
            raise ValueError("No objective specified")
        org_json.objectives = []
        for obj in objectives:
            obj_reac = org_json.get_reaction(obj["name"])
            if obj_reac is None:
                raise ValueError("Objective not in model: " + obj["name"])

            if obj_reac.disabled:
                return HttpResponseBadRequest("Objective disabled: " + obj_reac.name)

            org_json.objectives.append(JsonModel.Objective(**obj))

        if simtype == "mba" or simtype == "sa":
            org_json.design_objectives = []
            for obj in design_objectives:
                obj_reac = org_json.get_reaction(obj["name"])
                if obj_reac is None:
                    raise ValueError("Design objective not in model: " + obj["name"])

                if obj_reac.disabled:
                    return HttpResponseBadRequest("Design objective disabled: " + obj_reac.name)

                org_json.design_objectives.append(JsonModel.Objective(**obj))

            org_json.target_reactions = []
            if simtype == "sa":
                for obj in target_reactions:
                    obj_reac = org_json.get_reaction(obj["name"])
                    if obj_reac is None:
                        raise ValueError("Target reaction not in model: " + obj["name"])

                    if obj_reac.disabled:
                        return HttpResponseBadRequest("Target reaction disabled: " + obj_reac.name)

                    org_json.target_reactions.append(JsonModel.Objective(**obj))

    except ValueError as e:
        return HttpResponseBadRequest("Model error: " + str(e))

    org = org_json.to_model()

    try:
        fba = org.fba()
    except ValueError as e:
        return HttpResponseBadRequest("FBA error: " + str(e))

    if simtype == "fba":
        if dry_run:
            return HttpResponse(
                json.dumps({"solution": fba.get_status()}),
                content_type="application/json")

        full_g, nodeIDs = calc_reactions(org, fba)

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

        # Work around a bug in nx 1.11 breaking nx.to_agraph function
        graph = nx.drawing.nx_agraph.to_agraph(g)
        graph.graph_attr.update(splines=True, overlap=False, rankdir="LR")
        graph.node_attr.update(style="filled", colorscheme="pastel19")
        outgraph = str(graph.to_string())
        outgraph = pygraphviz.AGraph(outgraph)
        #outgraph = outgraph.draw(format="svg", prog="dot")

        if output_format == "json":
            return HttpResponse(json.dumps(
                {"graph": None,
                 "graphstr": str(graph.to_string()),
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
    elif simtype == "mba":
        # Determine optimal steps
        obj = org.obj[:]
        dobj = org.design_obj[:]

        # Absolute
        design_result = [["x"], [dobj[0][0]]]
        design_result[0] += [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
        sim_result = fba.design_fba(design_result[0][1:])
        dflux = [list(str(x[fba.reac_names.index(obj[0][0])]) for x in sim_result),
                 list(x[fba.reac_names.index(dobj[0][0])] for x in sim_result)]

        design_result[1] += dflux[1]

        if output_format == "json":
            return HttpResponse(json.dumps(
                {
                "solution": fba.get_status(),
                "flux": design_result
                }),
                content_type="application/json"
            )
        else:
            return HttpResponseBadRequest("Unknown format")
    elif simtype == "sa":
        # Percentage
        # 1. WT conditions
        obj = org.obj[:]
        dobj = org.design_obj[:]
        tobj = org_json.target_reactions[0].name
        obj_r = org.get_reaction(obj[0][0])
        target = org.get_reaction(tobj)
        target_flux = fba.flux[fba.reac_names.index(target.name)]

        design_result = [["x"], [obj[0][0]], [dobj[0][0]], ["Yield"]]

        # 2. Limit target reaction
        for limit in [1.0, 0.8, 0.6, 0.4, 0.2, 0.0]:
            # Limit to a % of the target flux
            target.constraint = (target.constraint[0], target_flux * limit)

            # Optimize limited growth
            org.obj = obj
            growth = org.fba().Z

            # Optimize production
            obj_r.constraint = (growth, growth)
            org.obj = dobj
            production = org.fba().Z

            # Reset production constraint
            obj_r.constraint = (0, None)

            design_result[0].append(str(int(limit * 100)) + "%")
            design_result[1].append(round(growth, 4))
            design_result[2].append(round(production, 4))
            design_result[3].append(round(growth * production, 4))

        if output_format == "json":
            return HttpResponse(json.dumps(
                {
                "solution": fba.get_status(),
                "flux": design_result
                }),
                content_type="application/json"
            )
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
        raise BadRequest("Request incomplete")

    try:
        model = DesignModel.objects.get(user=UserProfile.get_profile(request.user), pk=pk)
    except ObjectDoesNotExist:
        raise BadRequest("Bad Model")

    try:
        revision = request.POST["revision"]
        try:
            revision = Revision.objects.get(model=model, pk=revision)
        except ObjectDoesNotExist:
            raise BadRequest("Bad Revision")
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
                        raise BadRequest("Objective not in model: " + obj["name"])

                    org.objectives.append(JsonModel.Objective(**obj))
    except ValueError as e:
        raise BadRequest("Model error: " + str(e))

    if len(changes) == 0:
        raise BadRequest("Model not saved: No changes found")

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

            return {'success': True, 'url': reverse("cyanodesign:design", kwargs={"pk":dm.pk})}
        else:
            form_html = render_crispy_form(form, context=request)
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
                    try:
                        ss.write(chunk.decode("utf-8"))
                    except UnicodeDecodeError:
                        form.add_error("file", "File does not have UTF-8 encoding")
                        form_html = render_crispy_form(form, context=request)
                        return {'success': False, 'form_html': form_html}

                try:
                    ss.seek(0)
                    model = Metabolism(ss)

                    if len(model.enzymes) == 0 or len(model.mets) == 0:
                        raise ValueError("Model is empty")
                except:
                    form.add_error("file", "Not a valid model")
                    form_html = render_crispy_form(form, context=request)
                    return {'success': False, 'form_html': form_html}

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

    return BadRequest()


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
