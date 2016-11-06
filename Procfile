web: gunicorn formalign.wsgi --log-file -
worker: celery -A formalign worker -l info
worker: celery -A formalign beat -S djcelery.schedulers.DatabaseScheduler -l info