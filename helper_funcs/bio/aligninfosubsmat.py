from Bio.Align import AlignInfo


class AlignInfoSubsMat(AlignInfo.SummaryInfo):
    """
    subclasses SummaryInfo and overrides gap_consensus to use substitution matrices
    """
