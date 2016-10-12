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

Feature: Formalign alignment validation
  Checks that alignments are properly validated on submission

  Scenario Outline: User visits the Formalign home page and submits a series of invalid alignments
    Given a user visits the URL "/"
    When a custom <sequence type> alignment: "<alignment>" is submitted
    Then the server's response status code is 200
    And the url is the home url
    And the page title is "Formalign.eu Home"
    And the user should see the error message: <error message>

  Examples: DNA
    | sequence type | alignment                  | error message                 |
    | DNA           | empty                      | empty error                   |
    | DNA           | invalid characters         | character error               |
    | DNA           | too few sequences          | less than two sequences error |
    | DNA           | different sequence lengths | alignment error               |
    | DNA           | invalid FASTA format       | format error                  |

  Examples: protein
    | sequence type | alignment                  | error message                 |
    | protein       | empty                      | empty error                   |
    | protein       | invalid characters         | character error               |
    | protein       | too few sequences          | less than two sequences error |
    | protein       | different sequence lengths | alignment error               |
    | protein       | invalid FASTA format       | format error                  |
