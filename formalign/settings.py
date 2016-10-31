"""
=====================================================================
Formalign.eu format and display multiple sequence alignments
Copyright (C) 2016 Stefan Dieterle
e-mail: golgoths@yahoo.fr

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
=====================================================================

Django settings for formalign project.

Generated by 'django-admin startproject' using Django 1.9.

For more information on this file, see
https://docs.djangoproject.com/en/1.9/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/1.9/ref/settings/
"""

import os
import json
import dj_database_url
from unittest import TestCase
from django.contrib.staticfiles.testing import StaticLiveServerTestCase
from django.core.exceptions import ImproperlyConfigured
from datetime import timedelta


# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Load settings from json settings file
with open(os.path.join(BASE_DIR, 'settings.json')) as secrets_file:
    secrets = json.load(secrets_file)


def json_setting(setting, secrets=secrets):
    """Get secret setting or fail with ImproperlyConfigured"""
    try:
        return secrets[setting]
    except KeyError:
        raise ImproperlyConfigured("Set the %s setting" % setting)

# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/1.9/howto/deployment/checklist/

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = json_setting('SECRET_KEY')

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = json_setting('DEBUG')

# Application definition

INSTALLED_APPS = [
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'formalign',
    'base',
    'debug_toolbar',
    'django_nose',
    'djcelery',
    'randomslugfield',
    'behave_django',
]

TEST_RUNNER = 'django_nose.NoseTestSuiteRunner'

MIDDLEWARE_CLASSES = [
    'django.middleware.security.SecurityMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.auth.middleware.SessionAuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
]

ROOT_URLCONF = 'formalign.urls'

TEMPLATES = [
    {'BACKEND': 'django.template.backends.jinja2.Jinja2',
     'DIRS': [os.path.join(BASE_DIR, 'templates/jinja2')],
     'APP_DIRS': True,
     'OPTIONS': {
         'environment': 'formalign.jinja2.environment',
     }
     },
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [os.path.join(BASE_DIR, 'templates/base')],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.contrib.auth.context_processors.auth',
                'django.contrib.messages.context_processors.messages',
            ],
        },
    },

]

WSGI_APPLICATION = 'formalign.wsgi.application'

# Database
# https://docs.djangoproject.com/en/1.9/ref/settings/#databases
HEROKU = json_setting('HEROKU')
if HEROKU:
    db_from_env = dj_database_url.config(conn_max_age=500)
    DATABASES = {'default': db_from_env}
else:
    DATABASES = {
        'default': {
            'ENGINE': json_setting('DB_ENGINE'),
            'NAME': json_setting('DB_NAME'),
            'USER': json_setting('DB_USER'),
            "PASSWORD": json_setting('DB_PASSWORD'),
            'HOST': json_setting('DB_HOST'),
            'TEST': {
                'NAME': json_setting('TEST_DB_NAME'),
            }
        }
    }

# Password validation
# https://docs.djangoproject.com/en/1.9/ref/settings/#auth-password-validators

AUTH_PASSWORD_VALIDATORS = [
    {
        'NAME': 'django.contrib.auth.password_validation.UserAttributeSimilarityValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.MinimumLengthValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.CommonPasswordValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.NumericPasswordValidator',
    },
]

# Internationalization
# https://docs.djangoproject.com/en/1.9/topics/i18n/

LANGUAGE_CODE = 'en-us'

TIME_ZONE = 'UTC'

USE_I18N = True

USE_L10N = True

USE_TZ = True

# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/1.9/howto/static-files/

STATIC_URL = '/static/'
STATIC_ROOT = os.path.join(BASE_DIR, 'static')
STATICFILES_DIRS = [os.path.join(BASE_DIR, 'base/static'), ]

# Access control and SSL
ALLOWED_HOSTS = json_setting('ALLOWED_HOSTS')
SECURE_SSL_REDIRECT = json_setting('SECURE_SSL_REDIRECT')
SECURE_PROXY_SSL_HEADER = json_setting('SECURE_PROXY_SSL_HEADER')

# Celery configuration
BROKER_URL = json_setting('BROKER_URL')
CELERY_RESULT_BACKEND = json_setting('CELERY_RESULT_BACKEND')
BROKER_TRANSPORT_OPTIONS = {
    'fanout_prefix': True,
    'fanout_patterns': True,
    'visibility_timeout': 3600
}
CELERY_ACCEPT_CONTENT = ['json']
CELERY_TASK_SERIALIZER = 'json'
CELERY_RESULT_SERIALIZER = 'json'

# Celery beat configuration
CELERYBEAT_SCHEDULE = {
    'remove_more_than_week_old_alignments': {
        'task': 'base.tasks.clean_alignments',
        'schedule': timedelta(hours=1),
    },
}

# Selenium configuration
CHROME_DRIVER = json_setting('CHROME_DRIVER')
FIREFOX_BINARY = json_setting('FIREFOX_BINARY')
SERVER_URL = json_setting('SERVER_URL')
if SERVER_URL == 'liveserver':
    TEST_CASE = StaticLiveServerTestCase
else:
    TEST_CASE = TestCase

# Heroku configuration
if HEROKU:
    MIDDLEWARE_CLASSES.insert(1, 'whitenoise.middleware.WhiteNoiseMiddleware')
    BROKER_URL = os.environ['REDIS_URL']
    CELERY_RESULT_BACKEND = os.environ['REDIS_URL']
