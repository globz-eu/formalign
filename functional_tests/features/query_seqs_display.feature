Feature: Formalign query sequences display page
  Checks that the query sequences display page displays the sequences with the
  correct formatting

  Scenario: User visits the Formalign home page and submits a protein alignment
    Given a user visits the URL "/"
    When a protein alignment: "spa protein alignment" is submitted
    Then the server's response status code is 200
    And the user is redirected to the "sequence display" page
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
    And there is a "Formalign.eu" button with "/" href

  Scenario: User submits a custom alignment and renders it
    Given a user visits the URL "/"
    When a protein alignment: "spa protein alignment" is submitted
    And the user clicks the "Render" submit button
    Then the server's response status code is 200
    And the user is redirected to the "alignment display" page
    And there is a "Formalign.eu" button with "/" href

  @pending
  Scenario Outline: User submits a custom alignment with a non default
  consensus selected
    Given a user visits the URL "/"
    When <consensus_type> with <consensus_choice> is selected
    And a custom protein alignment: "spa protein alignment" is submitted
    Then the server's response status code is 200
    And the user is redirected to the "sequence display" page
    And there are sequences displayed
    And there is a consensus sequence named <consensus_choice>
    And the consensus sequence for <consensus_choice> is as expected

  Examples: % identity
    | consensus_type      | consensus_choice |
    | identity            | consensus 70%    |
    | identity            | consensus 100%   |
    | substitution matrix | blosum 62        |
    | substitution matrix | pam 250          |
