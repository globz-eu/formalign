Feature: Demo alignment
  Tests submission of the demo alignment

  Scenario: User visits the Formalign home page and submits the demo alignment
    Given a user visits the URL "/"
    When the user clicks the "Demo" submit button
    Then the server's response status code is 200
    And the user is redirected to the "sequence display" page
    And the current URL is the "sequence display" URL
    And there are demo sequences displayed
    And the sequences are displayed in lines of 80 characters
    And the correct demo sequences are displayed
    And the correct demo sequence metadata are displayed
    And there is a "Render" button with "get" method and "align-display" action
    And the action URL of the "Render" button contains a 16 character slug

  Scenario: User submits the demo alignment and renders it
    Given a user visits the URL "/"
    When the user clicks the "Demo" submit button
    And the user clicks the "Render" submit button
    Then the server's response status code is 200
    And the user is redirected to the "alignment display" page
