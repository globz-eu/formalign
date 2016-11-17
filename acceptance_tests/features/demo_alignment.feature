# =====================================================================
# Formalign.eu format and display multiple sequence alignments
# Copyright (C) 2016 Stefan Dieterle
# e-mail: golgoths@yahoo.fr
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# =====================================================================

Feature: Demo alignment
  Tests submission of the demo alignment

  Scenario Outline: User visits the Formalign home page and submits the demo alignment
    Given a user visits the URL "/" with "<browser>"
    When the user clicks the "Demo" submit button
    Then the user is redirected to the "sequence display" page
    And the current URL is the "sequence display" URL
    And there are demo sequences displayed
    And the sequences are displayed in lines of 80 characters
    And the correct demo sequences are displayed
    And the correct demo sequence metadata are displayed
    And there is a "Render" button with "get" method and "align-display" action
    And the action URL of the "Render" button contains a 16 character slug

  Examples: browsers
     | browser |
     | Chrome  |
     | Firefox |

  Scenario Outline: User submits the demo alignment and renders it
    Given a user visits the URL "/" with "<browser>"
    When the user clicks the "Demo" submit button
    Then the user is redirected to the "sequence display" page
    When the user clicks the "Render" submit button
    Then the user is redirected to the "alignment display" page

  Examples: browsers
     | browser |
     | Chrome  |
     | Firefox |
