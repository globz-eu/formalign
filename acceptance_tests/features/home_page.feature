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

  Scenario Outline: User checks basic elements of Formalign home page
    Given a user visits the URL "/" with "<browser>"
    Then the user is on the "home" page
    And the current URL is the "home" URL
    And the brand text says "Formalign.eu"

  Examples: browsers
     | browser |
     | Chrome  |
     | Firefox |

  Scenario Outline: User checks basic form elements of Formalign home page
    Given a user visits the URL "/" with "<browser>"
    Then there is a form with a label matching "Paste in your alignment:\n\(FASTA, clustalw, stockholm or phylip\)"
    And there is a text area with a placeholder saying "Alignment (FASTA, clustalw, stockholm or phylip)"
    And there are radio buttons labeled "Input sequence type:"
    And there is a "DNA" input sequence type radio button
    And there is a "Protein" input sequence type radio button
    And the "DNA" button is checked
    And the "Protein" button is not checked
    And there is a "Identity" consensus type radio button
    And there is a "Substitution Matrix" consensus type radio button
    And the "Identity" button is checked
    And the "Substitution Matrix" button is not checked
    And there is a "Demo" button named "custom_data" with the value "demo"
    And there is a "Submit" button named "custom_data" with the value "custom"
    And there is a "Formalign.eu" button with "/" href

  Examples: browsers
     | browser |
     | Chrome  |
     | Firefox |

  Scenario Outline: User checks basic functionality of radio buttons
    Given a user visits the URL "/" with "<browser>"
    When the user clicks the "<active>" radio button
    Then the "<active>" button is checked
    And the "<inactive>" button is not checked

  Examples: buttons
    | browser | active              | inactive            |
    | Chrome  | Protein             | DNA                 |
    | Chrome  | DNA                 | Protein             |
    | Chrome  | Identity            | Substitution Matrix |
    | Chrome  | Substitution Matrix | Identity            |
    | Firefox | Protein             | DNA                 |
    | Firefox | DNA                 | Protein             |
    | Firefox | Identity            | Substitution Matrix |
    | Firefox | Substitution Matrix | Identity            |
