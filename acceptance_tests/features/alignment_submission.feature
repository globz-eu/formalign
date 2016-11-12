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

@pending
@current
Feature: Protein alignment submission
  Checks the submission of a protein alignment

  Scenario: User submits a protein alignment
    Given a user visits the URL "/"
    When the user clicks the "Protein" radio button
    And the user pastes a protein alignment in the form text area
    And the user clicks the "Submit" submit button
    Then the user is redirected to the "sequence display" page
    And there are protein sequences displayed
    And the sequences are displayed in lines of 80 characters
    And the correct protein sequences are displayed
    And the correct protein sequence metadata are displayed
    And there is a consensus sequence displayed
    And the correct consensus sequence is displayed
    And the correct consensus sequence metadata is displayed

  Scenario: User submits a protein alignment and renders it
    Given a user visits the URL "/"
    When the user clicks the "Protein" radio button
    And the user pastes a protein alignment in the form text area
    And the user clicks the "Submit" submit button
    And the user clicks the "Render" submit button
    Then the user is redirected to the "alignment display" page
