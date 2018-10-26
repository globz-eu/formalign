import requests
from behave import use_step_matcher


use_step_matcher('re')


def before_scenario(context, scenario):
    context.client = requests.Session()


def after_scenario(context, scenario):
    context.client.close()
