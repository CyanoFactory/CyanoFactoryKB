from celery import current_task
from time import sleep
from celery.task.base import task, Task

from celery import current_app as celery
from celery.signals import task_sent

class ReportingTask(Task):
    abstract = True

    def after_return(self, status, retval, task_id, args, kwargs, einfo):
        print("after return",retval,task_id,args,kwargs)
        current_task.update_state(state=status,
                                  meta={'result': retval, 'user': kwargs["user"], 'description': kwargs["description"]})

@task(base = ReportingTask)
def add(x, y, user = None, description = None):
    current_task.update_state(state='PROGRESS',
                        meta={'current': 1, 'total': 1, 'user': user, 'description': description})
    
    #self.user = user.pk
    
    return x + y

@task(base = ReportingTask)
def progress(x, user = None, description = None):
    for i in range(x):
        sleep(5)
        current_task.update_state(state='PROGRESS',
                        meta={'current': i, 'total': x, 'user': user, 'description': description})
    return x

# via http://stackoverflow.com/questions/9824172
def update_sent_state(sender=None, id=None, **kwargs):
    # the task may not exist if sent using `send_task` which
    # sends tasks by name, so fall back to the default result backend
    # if that is the case.
    task = celery.tasks.get(sender)
    backend = task.backend if task else celery.backend

    backend.store_result(id, None, "SENT")
task_sent.connect(update_sent_state)
