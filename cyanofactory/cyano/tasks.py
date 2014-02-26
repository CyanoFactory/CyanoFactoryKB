"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from celery import current_task
from celery.task.base import task

from cyano.models import UserProfile
from celery.app.task import Task

class ReportingTask(Task):
    abstract = True

    def after_return(self, status, retval, task_id, args, kwargs, einfo):
        current_task.update_state(state=status,
                                  meta={'result': retval, 'user': kwargs["user"], 'message': kwargs["reason"]})

class NotifyProgress:
    def __init__(self, user):
        if isinstance(user, UserProfile):
            self.user = user.pk
        else:            
            try:
                self.user = int(user)
            except ValueError:
                raise
    
    def __call__(self, current, total, message):
        current_task.update_state(state='PROGRESS',
                        meta={'current': current, 'total': total, 'user': self.user, 'message': message})

def handle_task(obj, filename):
    obj.notify_progress = NotifyProgress(obj.user)
    with open(filename, "r") as f:
        obj.parse(f)
    obj.apply()

@task(base = ReportingTask)
def genbank(filename, wid, user, reason, chromosome, name):
    from bioparser.genbank import Genbank

    g = Genbank(wid = wid, user = user, reason = reason,
                chromosome = chromosome, name = name)

    handle_task(g, filename)

    return True

@task(base = ReportingTask)
def proopdb(filename, wid, user, reason):
    from bioparser.proopdb import ProOpDB
    
    p = ProOpDB(wid = wid, user = user, reason = reason)

    handle_task(p, filename)

    return True

@task(base = ReportingTask)
def sbml(filename, wid, user, reason):
    from bioparser.sbml import SBML
    
    s = SBML(wid = wid, user = user, reason = reason)
    
    handle_task(s, filename)

    return True

@task(base = ReportingTask)
def interproscan(filename, wid, user, reason):
    from bioparser.interproscan import InterProScan

    i = InterProScan(wid = wid, user = user, reason = reason)

    handle_task(i, filename)

    return True

@task(base = ReportingTask)
def fastafeature(filename, wid, user, reason, chromosome, feature_type):
    from bioparser.fastafeature import FastaFeature

    i = FastaFeature(wid=wid, user=user, reason=reason, chromosome=chromosome, feature_type=feature_type)

    handle_task(i, filename)

    return True
