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

Feature: Formalign home page
  Checks that the Formalign home page is displayed and contains all expected elements

  Scenario: User checks basic elements of Formalign home page
    Given a user visits the URL "/"
    When the user looks at the page
    Then the user is on the "home" page
    Then the server's response status code is 200
    And the current URL is the "home" URL
    And the brand text says "Formalign.eu"

  Scenario: User checks basic form elements of Formalign home page
    Given a user visits the URL "/"
    When the user looks at the page
    Then the server's response status code is 200
    And there is a form with a label matching "Paste in your alignment:\(FASTA, clustalw, stockholm or phylip\)"
    And there is a text area with a placeholder saying "Alignment (FASTA, clustalw, stockholm or phylip)"
    And there are radio buttons labeled "Input sequence type:"
    And there is a "DNA" input sequence type radio button
    And there is a "Protein" input sequence type radio button
    And the "DNA" button is checkable
    And the "Protein" button is checkable
    And the "DNA" button is checked
    And the "Protein" button is not checked
    And there is a "Demo" button named "custom_data" with the value "demo"
    And there is a "Submit" button named "custom_data" with the value "custom"
    And there is a "Formalign.eu" button with "/" href

  Scenario: User checks advanced form elements of Formalign home page
    Given a user visits the URL "/"
    When the user looks at the page
    Then the server's response status code is 200
    And there is a "Identity" consensus type radio button
    And there is a "Substitution Matrix" consensus type radio button
    And the "Identity" button is checkable
    And the "Substitution Matrix" button is checkable
    And the "Identity" button is checked
    And the "Substitution Matrix" button is not checked

  @pending
  Scenario Outline: User checks consensus choices
    Given a user visits the URL "/"
    When the user looks at the page
    Then the server's response status code is 200
    And there is a "consensus" dropdown menu
    And the "consensus" dropdown menu contains <consensus_choice>

  Examples: % identity
    | consensus_choice |
    | consensus 70%    |
    | consensus 100%   |

  @pending
  Scenario Outline: User does not see hidden consensus choices
    Given a user visits the URL "/"
    When the user looks at the page
    Then the server's response status code is 200
    And the "consensus" dropdown menu contents <hidden_consensus_choice> are hidden

  Examples: substitution matrix
    | hidden_consensus_choice |
    | blosum 62               |
    | pam 250                 |
