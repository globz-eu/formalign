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

@current
Feature: Formalign alignment validation
  Checks that alignments are properly validated on submission

  Scenario Outline: User visits the Formalign home page and submits a series of invalid alignments
    Given a user visits the URL "/" with "<browser>"
    When the user clicks the "<sequence type>" radio button
    And the user pastes a <sequence type> alignment: "<alignment>" in the form text area
    And the user clicks the "Submit" submit button
    Then the user stays on the "home" page
    And the current URL is the "home" URL
    And the user should see the error message: <error message>

  Examples: DNA
    | browser | sequence type | alignment                  | error message                 |
    | Chrome  | DNA           | empty                      | empty error                   |
    | Chrome  | DNA           | invalid characters         | character error               |
    | Chrome  | DNA           | too few sequences          | less than two sequences error |
    | Chrome  | DNA           | different sequence lengths | alignment error               |
    | Chrome  | DNA           | invalid FASTA format       | format error                  |
    | Firefox | DNA           | empty                      | empty error                   |
    | Firefox | DNA           | invalid characters         | character error               |
    | Firefox | DNA           | too few sequences          | less than two sequences error |
    | Firefox | DNA           | different sequence lengths | alignment error               |
    | Firefox | DNA           | invalid FASTA format       | format error                  |

  Examples: protein
    | browser | sequence type | alignment                  | error message                 |
    | Chrome  | protein       | empty                      | empty error                   |
    | Chrome  | protein       | invalid characters         | character error               |
    | Chrome  | protein       | too few sequences          | less than two sequences error |
    | Chrome  | protein       | different sequence lengths | alignment error               |
    | Chrome  | protein       | invalid FASTA format       | format error                  |
    | Firefox | protein       | empty                      | empty error                   |
    | Firefox | protein       | invalid characters         | character error               |
    | Firefox | protein       | too few sequences          | less than two sequences error |
    | Firefox | protein       | different sequence lengths | alignment error               |
    | Firefox | protein       | invalid FASTA format       | format error                  |
