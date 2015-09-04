"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""
from django.utils.functional import wraps
from django.contrib.auth.models import User
from django.core.exceptions import ObjectDoesNotExist, MultipleObjectsReturned
from cyano.models import Entry

import cyano.models as models
from cyano.helpers import render_queryset_to_response_error, getModel

def __assign_if_not_false(kw, key, obj):
    if obj:
        kw[key] = obj


def resolve_to_objects(function):
    @wraps(function)
    def wrapper(request, *args, **kw):
        # False instead of None, None is default value for passed but not set
        species_wid = kw.pop("species_wid", False)
        model_type = kw.pop("model_type", False)
        wid = kw.pop("wid", False)

        species = False
        model = False
        item = False

        if species_wid:
            species = models.Species.objects.for_wid(species_wid, get=False)
            try:
                species = species.get()
            except ObjectDoesNotExist as e:
                return render_queryset_to_response_error(
                    request,
                    error=404,
                    msg="The requested species \"{}\" was not found.".format(species_wid),
                    msg_debug=repr(e))
            except MultipleObjectsReturned as e:
                return render_queryset_to_response_error(
                    request,
                    error=500,
                    msg="Database error for species \"{}\". Please report this to an administrator!.".format(
                        species_wid),
                    msg_debug=repr(e))
        if model_type:
            model = getModel(model_type)
            # issubclass test needed to allow /add/Species/ but not /species_wid/Species/
            if model is None or (species_wid and not issubclass(model, models.SpeciesComponent)):
                return render_queryset_to_response_error(
                    request,
                    species=species,
                    error=404,
                    msg="The requested species \"{}\" has no model \"{}\".".format(species_wid, model_type))
        if model and wid:
            item = model.objects.for_species(species).for_wid(wid, get=False)
            try:
                item = item.get()
            except ObjectDoesNotExist as e:
                return render_queryset_to_response_error(
                    request,
                    error=404,
                    species=species,
                    model=model,
                    msg="The requested species \"{}\" has no item \"{}\" of type \"{}\".".format(species_wid, wid,
                                                                                                 model_type),
                    msg_debug=repr(e))
            except MultipleObjectsReturned as e:
                return render_queryset_to_response_error(
                    request,
                    error=500,
                    species=species,
                    model=model,
                    msg="Database error for species \"{}\" accessing item \"{}\" of type \"{}\"."
                        " Please report this to an administrator!".format(
                        species_wid, wid, model_type),
                    msg_debug=repr(e))

        # Prepare keyword arguments
        __assign_if_not_false(kw, "species", species)
        __assign_if_not_false(kw, "model", model)
        __assign_if_not_false(kw, "item", item)

        return function(request, *args, **kw)

    return wrapper


def permission_required(permission):
    """Checks whether a user has a permission for the specified entry
    """

    def decorator(function):
        @wraps(function)
        def wrapper(request, *args, **kw):
            from cyano.views import login

            species = kw.get("species", False)
            item = kw.get("item", False)

            user = request.user

            if not user.is_authenticated():
                # Special handling for guests
                user = User.objects.get(username="guest")
            else:
                pass

            # Admins always have full access
            if user.is_superuser:
                return function(request, *args, **kw)

            if species:
                species_old_cls = species.__class__
                species.__class__ = Entry
                allow_item = None

                if item:
                    item_old_cls = item.__class__
                    item = item.__class__ = Entry
                    allow_item = user.profile.has_perm(permission, item)
                    item.__class__ = item_old_cls

                allow_species = user.profile.has_perm(permission, species)
                species._class__ = species_old_cls

                if allow_species or allow_item:
                    # Has permission
                    return function(request, *args, **kw)
            else:
                # Global perm check
                if user.profile.has_perm(permission):
                    return function(request, *args, **kw)

            # No perm fallthrough
            if species:
                msg = "You need permission \"{}\" to access ".format(permission)

                if item:
                    msg += "item {} of ".format(item.wid)

                msg += "species {}".format(species.wid)
            else:
                msg = "You need global permission \"{}\"".format(permission)

            return login(
                request,
                species=species,
                error=403,
                message=msg,
                force_next=request.build_absolute_uri(),
                required_perm=permission)

        return wrapper

    return decorator


def ajax_required(f):
    """
    AJAX request required decorator
    use it in your views:

    @ajax_required
    def my_view(request):
    ....

    via http://djangosnippets.org/snippets/771/
    """
    def wrap(request, *args, **kwargs):
        from django.http import HttpResponseBadRequest
        if not request.is_ajax():
            return HttpResponseBadRequest("Ajax only")
        return f(request, *args, **kwargs)
    wrap.__doc__ = f.__doc__
    wrap.__name__ = f.__name__
    return wrap
