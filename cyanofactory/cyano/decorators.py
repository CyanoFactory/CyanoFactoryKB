import cyano.models as model
from cyano.helpers import get_queryset_object_or_404, render_queryset_to_response_error
from django.utils.functional import wraps
from django.contrib.auth.models import AnonymousUser

def read_permission_required(function):
    """Checks whether a user has a read_permission for the specified entry
    """
    @wraps(function)
    def wrapper(request, species_wid = None, wid = None, *args, **kw):
        if species_wid != None:
            species = get_queryset_object_or_404(model.Species.objects.filter(wid = species_wid))
        
        user = request.user
        
        if user.is_authenticated():
            # Do something for logged-in users.
            profile = user.profile
            if not profile.has_permission(species, model.PermissionEnum.READ_NORMAL):
                return render_queryset_to_response_error(
                    request,
                    species = species,
                    error = 403,
                    msg = "You don't have permission to access this item")
            else:
                return function(request, species_wid, *args, **kw)
        else:
            # Do something for anonymous users.
            # TODO: Maybe add a guest user?
            #main_group = model.GroupProfile.objects.get(group__name = "Guest")
            return render_queryset_to_response_error(
                request,
                species = species,
                error = 403,
                msg = "Please login")

    return wrapper
