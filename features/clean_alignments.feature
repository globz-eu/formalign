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


Feature: Remove old alignments
  Tests that alignments older than the threshold have been removed

  Scenario: User tries to display an alignment older than the threshold
    Given an old alignment with slug "ToOldAlignment01" was present
    And clean_alignments task has been run
    When a user visits the URL "/align-display/ToOldAlignment01/"
    Then the server's response status code is 404

  Scenario: User tries to display sequences of an alignment older than the threshold
    Given an old alignment with slug "ToOldAlignment01" was present
    And clean_alignments task has been run
    When a user visits the URL "/query-sequences/ToOldAlignment01/"
    Then the server's response status code is 404
