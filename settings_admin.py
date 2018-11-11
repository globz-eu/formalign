import os
from formalign.settings import *  # noqa

DATABASES['default']['USER'] = os.environ['DB_USER']  # noqa
DATABASES['default']['PASSWORD'] = os.environ['DB_PASSWORD']  # noqa
