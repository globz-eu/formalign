Feature: Formalign alignment display page
  Checks that the alignment display page displays the correct sequences
  with the correct alignment formatting

  Scenario Outline: User submits a custom alignment and renders it
    Given a user visits the URL "/" with "<browser>"
    When the user clicks the "Protein" radio button
    And the user pastes a protein alignment: "spa protein alignment" in the form text area
    And the user clicks the "Submit" submit button
    Then the user is redirected to the "sequence display" page
    When the user clicks the "Render" submit button
    Then the user is redirected to the "alignment display" page
    And the alignment is displayed with 80 characters per line in blocks of 10 with sequence IDs
    And the expected alignments are displayed
    And the expected consensus sequence is displayed
    And the sequence elements have the expected color classes
    And there is a "Formalign.eu" button with "/" href

  Examples: browsers
     | browser |
     | Chrome  |
     | Firefox |
