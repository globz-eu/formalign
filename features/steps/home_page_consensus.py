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

from behave import then

__author__ = 'Stefan Dieterle'


@then(r'there is a "(?P<dropdown>[^"]*)" dropdown menu')
def dropdown_present(context, dropdown):
    dropdown_menu = {'consensus': context.display.cssselect('input[id="id_consensus_select"]')[0].label.text_content()}
    assert 'consensus' == dropdown_menu[dropdown], 'Got %s' % dropdown_menu[dropdown]


@then(r'the "(?P<dropdown>[^"]*)" dropdown menu contains (?P<consensus_choice>.*)')
def check_dropdown_content(context, dropdown, consensus_choice):
    dropdown_menu = {'consensus': context.display.cssselect('input[id="id_consensus_select"]')[0]}
    assert consensus_choice in dropdown_menu[dropdown], 'Got %s' % dropdown_menu[dropdown]


@then(r'the "(?P<dropdown>[^"]*)" dropdown menu contents (?P<consensus_choice>.*) are hidden')
def check_dropdown_hidden_content(context, dropdown, consensus_choice):
    dropdown_menu = {'consensus': context.display.cssselect('input[id="id_consensus_select"]')[0]}
    assert consensus_choice in dropdown_menu[dropdown], 'Got %s' % dropdown_menu[dropdown]
