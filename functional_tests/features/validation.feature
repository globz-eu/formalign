Feature: Formalign alignment validation
  Checks that alignments are properly validated on submission

  Scenario Outline: User visits the Formalign home page and submits a series of invalid alignments
    Given a user visits the URL "/"
    When a <sequence type> alignment: "<alignment>" is submitted
    Then the user stays on the "home" page
    And the server's response status code is 200
    And the current URL is the "home" URL
    And the user should see the error message: <error message>

  Examples: DNA
    | sequence type | alignment                  | error message                 |
    | DNA           | empty                      | empty error                   |
    | DNA           | invalid characters         | character error               |
    | DNA           | too few sequences          | less than two sequences error |
    | DNA           | different sequence lengths | alignment error               |
    | DNA           | invalid FASTA format       | format error                  |

  Examples: protein
    | sequence type | alignment                  | error message                 |
    | protein       | empty                      | empty error                   |
    | protein       | invalid characters         | character error               |
    | protein       | too few sequences          | less than two sequences error |
    | protein       | different sequence lengths | alignment error               |
    | protein       | invalid FASTA format       | format error                  |
