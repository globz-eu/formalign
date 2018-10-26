Feature: Formalign home page
  Checks that the Formalign home page is displayed and contains all expected elements

  Scenario Outline: User checks basic elements of Formalign home page
    Given a user visits the URL "/" with "<browser>"
    Then the user is on the "home" page
    And the current URL is the "home" URL
    And the brand text says "Formalign.eu"

  Examples: browsers
     | browser |
     | Chrome  |
     | Firefox |

  Scenario Outline: User checks basic form elements of Formalign home page
    Given a user visits the URL "/" with "<browser>"
    And active and inactive buttons
      | active   | inactive            |
      | DNA      | Protein             |
      | Identity | Substitution Matrix |

    Then the user is on the "home" page
    And there is a form with a label matching "Paste in your alignment:\n\(FASTA, clustalw, stockholm or phylip\)"
    And there is a text area with a placeholder saying "Alignment (FASTA, clustalw, stockholm or phylip)"
    When the user looks at the radio buttons
    Then there are radio buttons labeled "Input sequence type:"
    And there is a "DNA" input sequence type radio button
    And there is a "Protein" input sequence type radio button
    And there is a "Identity" consensus type radio button
    And there is a "Substitution Matrix" consensus type radio button
    And the active radio buttons are checked
    And the inactive radio buttons are not checked
    And there is a "Demo" button named "custom_data" with the value "demo"
    And there is a "Submit" button named "custom_data" with the value "custom"
    And there is a "Formalign.eu" button with "/" href

  Examples: browsers
     | browser |
     | Chrome  |
     | Firefox |

  Scenario Outline: User checks basic functionality of radio buttons
    Given a user visits the URL "/" with "<browser>"
    Then the user is on the "home" page
    Given active and inactive buttons
      | active              | inactive            |
      | Protein             | DNA                 |
      | DNA                 | Protein             |
      | Identity            | Substitution Matrix |
      | Substitution Matrix | Identity            |

    When the user clicks the active radio buttons
    Then the active radio buttons are checked
    And the inactive radio buttons are not checked

  Examples: browsers
     | browser |
     | Chrome  |
     | Firefox |

  @pending
  Scenario Outline: User checks consensus choices
    Given a user visits the URL "/" with "<browser>"
    Then the user is on the "home" page
    And there is a "consensus" dropdown menu
    And the "consensus" dropdown menu contains <consensus_choice>

  Examples: % identity
    | browser | consensus_choice |
    | Chrome  | consensus 70%    |
    | Chrome  | consensus 100%   |
    | Firefox | consensus 70%    |
    | Firefox | consensus 100%   |

  @pending
  Scenario Outline: User does not see hidden consensus choices
    Given a user visits the URL "/" with "<browser>"
    Then the user is on the "home" page
    And the "consensus" dropdown menu contents <hidden_consensus_choice> are hidden

  Examples: substitution matrix
    | browser | hidden_consensus_choice |
    | Chrome  | blosum 62               |
    | Chrome  | pam 250                 |
    | Firefox | blosum 62               |
    | Firefox | pam 250                 |

  @pending
  Scenario Outline: User checks consensus choices after choosing consensus type
    Given a user visits the URL "/" with "<browser>"
    When the user clicks the "<active>" radio button
    Then the user is on the "home" page
    And there is a "consensus" dropdown menu
    And the "consensus" dropdown menu contains <consensus_choice>

  Examples: visible
    | browser | active              | consensus_choice |
    | Chrome  | Identity            | consensus 70%    |
    | Chrome  | Identity            | consensus 100%   |
    | Firefox | Identity            | consensus 70%    |
    | Firefox | Identity            | consensus 100%   |
    | Chrome  | Substitution Matrix | blosum 62        |
    | Chrome  | Substitution Matrix | pam 250          |
    | Firefox | Substitution Matrix | blosum 62        |
    | Firefox | Substitution Matrix | pam 250          |

  @pending
  Scenario Outline: User checks consensus choices after choosing consensus type
    Given a user visits the URL "/" with "<browser>"
    When the user clicks the "<active>" radio button
    Then the user is on the "home" page
    And there is a "consensus" dropdown menu
    And the "consensus" dropdown menu contents <hidden_consensus_choice> are hidden

  Examples: hidden
    | browser | active              | hidden_consensus_choice |
    | Chrome  | Identity            | blosum 62               |
    | Chrome  | Identity            | pam 250                 |
    | Firefox | Identity            | blosum 62               |
    | Firefox | Identity            | pam 250                 |
    | Chrome  | Substitution Matrix | consensus 70%           |
    | Chrome  | Substitution Matrix | consensus 100%          |
    | Firefox | Substitution Matrix | consensus 70%           |
    | Firefox | Substitution Matrix | consensus 100%          |
