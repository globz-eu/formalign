@skip
@pending
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