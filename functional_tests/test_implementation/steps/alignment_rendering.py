from behave import then, use_step_matcher
from formalign.settings import TEST
from helper_funcs.helpers_test import file_to_string
from lxml import html
from io import StringIO
from functional_tests.test_implementation.alignment_rendering import alignment_formatting, get_displayed_seqs


use_step_matcher('re')


@then(r'the alignment is displayed with 80 characters per line in blocks of 10 with sequence IDs')
def check_alignment_formatting(context):
    """
    tests that the alignments are displayed with the correct formatting
    :param context: behave context
    """
    tables = None
    seqs_meta = file_to_string('spa_protein_alignment_meta.txt').splitlines()
    if TEST == 'acceptance':
        align_display_elem = context.browser.find_element_by_css_selector('.align_display')
        align_display = context.browser.execute_script('return arguments[0].innerHTML', align_display_elem)
        align_display_html = html.parse(StringIO(align_display)).getroot()
        tables = align_display_html.find_class('align_table')
    elif TEST == 'functional':
        tables = context.display.find_class('align_table')

    alignment_formatting(seqs_meta, tables)


@then(r'the expected alignments are displayed')
def check_alignment_sequences(context):
    """
    tests that the expected sequences are displayed
    :param context: behave context
    """
    alignment = file_to_string('spa_protein_alignment_seqs.txt')
    alignment_list = [[c for c in a] for a in alignment.split('\n')[:-1]]
    elems = None
    if TEST == 'acceptance':
        align_display_elem = context.browser.find_element_by_css_selector('.align_display')
        align_display = context.browser.execute_script('return arguments[0].innerHTML', align_display_elem)
        align_display_html = html.parse(StringIO(align_display)).getroot()
        elems = align_display_html.cssselect('tr')
    elif TEST == 'functional':
        elems = context.display.cssselect('tr')

    re_seqs = get_displayed_seqs(elems, len(alignment_list))

    for i, al_li in enumerate(alignment_list):
        assert al_li == re_seqs[i], 'expected: %s\n got: %s' % (al_li, re_seqs[i])


@then(r'the expected consensus sequence is displayed')
def check_consensus_sequence(context):
    """
    tests that the expected consensus sequence is displayed
    :param context: behave context
    """
    alignment = file_to_string('spa_protein_alignment_consens.txt')
    alignment_list = [[c for c in a] for a in alignment.split('\n')[:-1]]
    elems = None
    if TEST == 'acceptance':
        align_display_elem = context.browser.find_element_by_css_selector('.align_display')
        align_display = context.browser.execute_script('return arguments[0].innerHTML', align_display_elem)
        align_display_html = html.parse(StringIO(align_display)).getroot()
        elems = align_display_html.cssselect('tr')
    elif TEST == 'functional':
        elems = context.display.cssselect('tr')

    cat_re_seq = get_displayed_seqs(elems, len(alignment_list), cons=True)

    cons_li = alignment_list[-1]
    assert cons_li == cat_re_seq, cat_re_seq


@then(r'the sequence elements have the expected color classes')
def check_alignment_sequences_annotation(context):
    """
    tests that the sequence elements (residues or bases) have the expected
    color classes
    :param context: behave context
    """
    alignment = file_to_string('spa_protein_alignment_seqs_annot.txt')
    alignment_list = [['residue S%s' % a for a in al] for al in alignment.split('\n')[:-1]]
    elems = None
    if TEST == 'acceptance':
        align_display_elem = context.browser.find_element_by_css_selector('.align_display')
        align_display = context.browser.execute_script('return arguments[0].innerHTML', align_display_elem)
        align_display_html = html.parse(StringIO(align_display)).getroot()
        elems = align_display_html.cssselect('tr')
    elif TEST == 'functional':
        elems = context.display.cssselect('tr')

    re_seqs = get_displayed_seqs(elems, len(alignment_list), annot=True)

    for i, al in enumerate(alignment_list):
        assert al == re_seqs[i], 'expected: %s\n got: %s' % (al, re_seqs[i])
