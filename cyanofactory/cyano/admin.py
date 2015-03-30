'''
Whole-cell knowledge base admin interface

Author: Jonathan Karr, jkarr@stanford.edu
Affiliation: Covert Lab, Department of Bioengineering, Stanford University
Last updated: 2012-07-17

Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
'''

from django.contrib import admin
from django.contrib.auth.admin import UserAdmin, GroupAdmin
from django.contrib.auth.models import User, Group
from cyano import models

''' User profile admin '''
class UserProfileInline(admin.StackedInline):
    model = models.UserProfile


class UserProfileAdmin(UserAdmin):
    inlines = [ UserProfileInline, ]


class GroupProfileInline(admin.StackedInline):
    import cyano.models as models
    model = models.GroupProfile


class GroupProfileAdmin(GroupAdmin):
    inlines = [ GroupProfileInline, ]


''' Register admins '''
admin.site.unregister(User)
admin.site.register(User, UserProfileAdmin)

admin.site.unregister(Group)
admin.site.register(Group, GroupProfileAdmin)

