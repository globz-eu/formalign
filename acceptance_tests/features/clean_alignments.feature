Feature: Remove old alignments
  Tests that alignments older than the threshold have been removed

  Scenario: User tries to display an alignment older than the threshold
    Given an old alignment with slug "ToOldAlignment01" was present
    And clean_alignments task has been run
    And a user visits the URL "/align-display/ToOldAlignment01/" with "Chrome"
    Then the user is redirected to the "404" page

  Scenario: User tries to display sequences of an alignment older than the threshold
    Given an old alignment with slug "ToOldAlignment01" was present
    And clean_alignments task has been run
    And a user visits the URL "/query-sequences/ToOldAlignment01/" with "Chrome"
    Then the user is redirected to the "404" page
