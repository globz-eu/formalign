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

Feature: Formalign alignment display page
  Checks that the alignment display page displays the correct sequences
  with the correct alignment formatting

  Scenario Outline: User submits a custom alignment and renders it
    Given a user visits the URL "/" with "<browser>"
    When the user clicks the "Protein" radio button
    And the user pastes a protein alignment: "spa protein alignment" in the form text area
    And the user clicks the "Submit" submit button
    Then the user is redirected to the "sequence display" page
    When the user clicks the "Render" submit button
    Then the user is redirected to the "alignment display" page
    And the alignment is displayed with 80 characters per line in blocks of 10 with sequence IDs
    And the expected alignments are displayed
    And the expected consensus sequence is displayed
    And the sequence elements have the expected color classes
    And there is a "Formalign.eu" button with "/" href

  Examples: browsers
     | browser |
     | Chrome  |
     | Firefox |
