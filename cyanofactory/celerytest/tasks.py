from celery import current_task
from time import sleep
from celery.task.base import task, Task
import datetime
from django.db.transaction import commit_on_success

class ReportingTask(Task):
    abstract = True

    def after_return(self, status, retval, task_id, args, kwargs, einfo):
        print("after return",retval,task_id,args,kwargs,einfo)
        current_task.update_state(state=status,
                                  meta={'result': retval, 'user': kwargs["user"], 'description': kwargs["description"]})

@task(base = ReportingTask)
def add(x, y, user = None, description = None):
    current_task.update_state(state='PROGRESS',
                        meta={'current': 1, 'total': 1, 'user': user, 'description': description})
    
    #self.user = user.pk
    
    return x + y

@task(base = ReportingTask)
@commit_on_success
def progress(x, user = None, description = None):
    for i in range(x):
        sleep(5)
        current_task.update_state(state='PROGRESS',
                        meta={'current': i, 'total': x, 'user': user, 'description': description})
    return x
