web: gunicorn formalign.wsgi --log-file -
worker: python manage.py celeryd -c 3 --beat --scheduler=djcelery.schedulers.DatabaseScheduler