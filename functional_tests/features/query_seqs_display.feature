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

Feature: Formalign query sequences display page
  Checks that the query sequences display page displays the sequences with the
  correct formatting

  Scenario: User visits the Formalign home page and submits a protein alignment
    Given a user visits the URL "/"
    When a protein alignment: "spa protein alignment" is submitted
    Then the server's response status code is 200
    And the user is redirected to the "sequence display" page
    And the current URL is the "sequence display" URL
    And there are protein sequences displayed
    And the sequences are displayed in lines of 80 characters
    And the correct protein sequences are displayed
    And the correct protein sequence metadata are displayed
    And there is a consensus sequence displayed
    And the correct consensus sequence is displayed
    And the correct consensus sequence metadata is displayed
    And there is a "Render" button with "get" method and "align-display" action
    And the action URL of the "Render" button contains a 16 character slug
    And there is a "Formalign.eu" button with "/" href

  Scenario: User submits a custom alignment and renders it
    Given a user visits the URL "/"
    When a protein alignment: "spa protein alignment" is submitted
    And the user clicks the "Render" submit button
    Then the server's response status code is 200
    And the user is redirected to the "alignment display" page
    And there is a "Formalign.eu" button with "/" href

  @pending
  Scenario Outline: User submits a custom alignment with a non default
  consensus selected
    Given a user visits the URL "/"
    When <consensus_type> with <consensus_choice> is selected
    And a custom protein alignment: "spa protein alignment" is submitted
    Then the server's response status code is 200
    And the user is redirected to the "sequence display" page
    And there are sequences displayed
    And there is a consensus sequence named <consensus_choice>
    And the consensus sequence for <consensus_choice> is as expected

  Examples: % identity
    | consensus_type      | consensus_choice |
    | identity            | consensus 70%    |
    | identity            | consensus 100%   |
    | substitution matrix | blosum 62        |
    | substitution matrix | pam 250          |
