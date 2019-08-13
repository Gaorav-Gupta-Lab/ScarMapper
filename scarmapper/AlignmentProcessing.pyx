"""
Cython version of Aligment Processing def.

@author: Dennis A. Simpson
          University of North Carolina
          Lineberger Comprehensive Cancer Center
          450 West Drive
          Chapel Hill, NC  27599-7295
@copyright: 2019
"""

# cython: language_level=3

__author__ = 'Dennis A. Simpson'
__version__ = "0.1.0"
__package__ = 'ScarMapper'

cpdef object alignment_processing(int cutposition, int gap_con_len, int gap_genome_len, object gapped_genomic, object gapped_consensus, object summary_data):
    """

    :param cutposition: 
    :param gap_con_len: 
    :param gap_genome_len: 
    :param gapped_genomic: 
    :param gapped_consensus: 
    :param summary_data: 
    :return: 
    """
    cdef int left_deletion = 0
    cdef int right_deletion = 0
    cdef int left_insertion = 0
    cdef int right_insertion = 0
    cdef int misaligned_insertion = 0

    cdef double lower_limit = 0.10 * gap_genome_len
    cdef double upper_limit = 0.90 * gap_genome_len
    cdef bint indel = False
    cdef bint stop = False
    cdef list microhomology, insertion, results_list, snv
    cdef int position, new_position, left_step, right_step
    cdef str c, g, temp_homology, tmp_snv, tmp_insertion
    microhomology = []
    insertion = []
    snv = []

    for position, (c, g) in enumerate(zip(gapped_consensus, gapped_genomic)):
        temp_homology = ""
        tmp_snv = ''
        tmp_insertion = ""

        if c == g or position < lower_limit or position > upper_limit:
            indel = False

        # Process Deletions and look for microhomology
        elif c == "-":
            if position <= cutposition:
                left_deletion += 1

            elif position > cutposition:
                right_deletion += 1

            if not indel:
                new_position = position

                while new_position < gap_con_len and gapped_consensus[new_position] == "-":
                    new_position += 1

                left_step = position
                right_step = new_position
                stop = False

                # find the microhomology
                while not stop:
                    if right_step >= gap_con_len or left_step >= new_position or gapped_consensus[right_step] == "N":
                        stop = True
                    elif gapped_genomic[left_step] == gapped_consensus[right_step]:
                        temp_homology += gapped_genomic[left_step]
                    else:
                        stop = True

                    right_step += 1
                    left_step += 1

            indel = True
            if temp_homology:
                microhomology.append(temp_homology)

        # Process insertions and capture inserted sequence.
        elif g == "-" and c != "N":

            if position <= cutposition:
                left_insertion += 1

            elif position > cutposition:
                right_insertion += 1

            if not indel:
                new_position = position
                gap_genome_len = len(gapped_genomic)
                while new_position < gap_genome_len and gapped_genomic[new_position] == "-" and \
                        gapped_consensus[new_position] != "N":

                    tmp_insertion += gapped_consensus[new_position]
                    new_position += 1
                if tmp_insertion:
                    insertion.append(tmp_insertion)
            indel = True

        elif c != g and c != "N" and g != "N":
            misaligned_insertion += 1
            if not indel:
                new_position = position

                while new_position < gap_con_len and gapped_genomic[new_position] != gapped_consensus[new_position] \
                        and gapped_genomic[new_position] != "-" and gapped_consensus[new_position] != "N" \
                        and gapped_consensus[new_position] != "-":

                    tmp_snv += gapped_consensus[new_position]
                    new_position += 1
                if tmp_snv:
                    snv.append(tmp_snv)
            indel = True

    micro_homology = ""
    gapped_con = gapped_consensus
    gapped_ref = gapped_genomic
    variants = ""
    insertion_seq = ""
    if microhomology:
        summary_data[8] += 1
        micro_homology = microhomology
    if snv:
        variants = snv
    if insertion:
        insertion_seq = insertion

    results_list = [left_deletion, right_deletion, left_insertion, right_insertion, micro_homology, insertion_seq,
                    variants, gapped_con, gapped_ref]

    if left_deletion > 0:
        summary_data[2] += 1
    if right_deletion > 0:
        summary_data[3] += 1
    if left_insertion > 0:
        summary_data[5] += 1
    if right_insertion > 0:
        summary_data[6] += 1

    return results_list, summary_data