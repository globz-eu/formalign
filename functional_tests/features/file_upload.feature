@skip
@pending
Feature: Submit alignment by uploading file
  Checks that alignments can be submitted by uploading a file

  Scenario: User submits a custom alignment after uploading it from a file and renders it
    Given a user visits the URL "/"
    When a custom protein alignment: "spa protein alignment" is uploaded
    Then the server's response status code is 200
    And the user is redirected to the "sequence display" page
    And there are sequences displayed
    And the sequences are displayed in lines of 80 characters
    And the correct custom sequences are displayed
    And the correct custom sequence metadata are displayed
