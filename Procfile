web: gunicorn formalign.wsgi --log-file -
worker: celery multi start w1 -A formalign -B --scheduler=djcelery.schedulers.DatabaseScheduler