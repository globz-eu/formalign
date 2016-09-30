web: gunicorn formalign.wsgi --log-file -
worker: celery -A formalign beat -S djcelery.schedulers.DatabaseScheduler -l info