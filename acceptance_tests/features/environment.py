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
"""

from behave import use_step_matcher
from selenium import webdriver
from formalign.settings import CHROME_DRIVER

__author__ = 'Stefan Dieterle'


use_step_matcher('re')


def before_scenario(context, scenario):
    context.browser = webdriver.Chrome(CHROME_DRIVER)


def after_scenario(context, scenario):
    context.browser.quit()
