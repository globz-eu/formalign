Feature: Formalign alignment display page
  Checks that the alignment display page displays the correct sequences
  with the correct alignment formatting

  Scenario: User submits a custom alignment and renders it
    Given a user visits the URL "/"
    When a protein alignment: "spa protein alignment" is submitted
    And the user clicks the "Render" submit button
    Then the server's response status code is 200
    And the user is redirected to the "alignment display" page
    And the alignment is displayed with 80 characters per line in blocks of 10 with sequence IDs
    And the expected alignments are displayed
    And the expected consensus sequence is displayed
    And the sequence elements have the expected color classes
    And there is a "Formalign.eu" button with "/" href
