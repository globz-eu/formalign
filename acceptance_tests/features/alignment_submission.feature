Feature: Protein alignment submission
  Checks the submission of a protein alignment

  Scenario Outline: User submits a protein alignment
    Given a user visits the URL "/" with "<browser>"
    When the user clicks the "Protein" radio button
    And the user pastes a protein alignment: "spa protein alignment" in the form text area
    And the user clicks the "Submit" submit button
    Then the user is redirected to the "sequence display" page
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

  Examples: browsers
     | browser |
     | Chrome  |
     | Firefox |

  Scenario Outline: User submits a protein alignment and renders it
    Given a user visits the URL "/" with "<browser>"
    When the user clicks the "Protein" radio button
    And the user pastes a protein alignment: "spa protein alignment" in the form text area
    And the user clicks the "Submit" submit button
    Then the user is redirected to the "sequence display" page
    When the user clicks the "Render" submit button
    Then the user is redirected to the "alignment display" page
    And the current URL is the "alignment display" URL

  Examples: browsers
     | browser |
     | Chrome  |
     | Firefox |
