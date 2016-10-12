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

@skip
Feature: Supported alignment formats
  Checks functionality of Formalign with different kinds of alignment formats

  Scenario Outline: User visits the Formalign home page and submits a series of alignments in various formats
    Given a user visits the URL "/"
    When a custom <sequence type> alignment: "<alignment>" is submitted as <format>
    Then the server's response status code is 200

  Examples: protein
    | sequence type | alignment | format    |
    | protein       | fasta     | fasta     |
    | protein       | clustal   | clustal   |
    | protein       | ig        | ig        |
    | protein       | nexus     | nexus     |
    | protein       | phylip    | phylip    |
    | protein       | stockholm | stockholm |