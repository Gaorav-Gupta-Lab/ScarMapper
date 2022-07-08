# cython: language_level=3

"""
Cython version of Sliding Window Module


@author: Dennis A. Simpson
          University of North Carolina
          Lineberger Comprehensive Cancer Center
          450 West Drive
          Chapel Hill, NC  27599-7295
@copyright: 2022

"""

__version__ = "0.5.0"

# from Valkyries import Tool_Box

cpdef sliding_window(str consensus, str target_region, int cutsite, int target_length, int lower_limit, int upper_limit,
                     object summary_data, list left_target_windows, list right_target_windows, str cutwindow,
                     str hr_donor):

    cdef int consensus_length = len(consensus)
    cdef int consensus_lft_junction = 0
    cdef int consensus_rt_junction = 0
    cdef int offset = 15
    cdef double upper_consensus_limit = consensus_length-15
    cdef str ldel, rdel, lft_test, rt_test

    target_lft_junction = cutsite
    target_rt_junction = cutsite

    ldel = ""
    rdel = ""

    cdef bint left_found = False
    cdef bint right_found = False
    cdef bint cut_found = False

    cdef int lft_position, rt_position, consensus_rt_position, consensus_lft_position, i

    '''
    Find the 5' junction.  Start at the cut position, derived from the target region, and move toward the 5'
    end of the read one nucleotide at a time using a 10 nucleotide sliding window.  The 3' position of
    the first perfect match of the window from the query and target defines the 5' junction.
    '''

    consensus_rt_position = consensus_length-10
    consensus_lft_position = consensus_rt_position-10

    while not left_found and consensus_lft_position > lower_limit:
        query_segment = consensus[consensus_lft_position:consensus_rt_position]
        for i, target_segment in enumerate(left_target_windows):
            if query_segment == target_segment:
                query_cutwindow = consensus[consensus_lft_position:consensus_rt_position]

                if query_cutwindow == cutwindow:
                    summary_data[6][1] += 1
                    return [], summary_data

                left_found = True
                target_lft_junction = cutsite-i
                consensus_lft_junction = consensus_rt_position
                ldel = target_region[target_lft_junction:cutsite]
                break

        consensus_lft_position -= 1
        consensus_rt_position -= 1

    '''
    Find the 3' junction.  Start at the cut position, derived from the target region, and move toward the 5'
    end of the read one nucleotide at a time using a 10 nucleotide sliding window.  One plus the 3' position of
    the first perfect match of the window from the query (FASTQ read 2) and target defines the 5' junction.  The
    query and targets have different numbering and the windows move in opposite directions.
    '''

    # Move to the expected cutsite position on the consensus from the 3' end.
    consensus_lft_position = 10
    consensus_rt_position = consensus_lft_position+10
    while not right_found and consensus_rt_position < upper_consensus_limit:
        query_segment = consensus[consensus_lft_position:consensus_rt_position]
        for i, target_segment in enumerate(right_target_windows):
            if query_segment == target_segment:
                right_found = True
                target_rt_junction = cutsite+i
                consensus_rt_junction = consensus_lft_position
                rdel = target_region[cutsite:target_rt_junction]
                break

        # increment consensus window
        consensus_lft_position += 1
        consensus_rt_position += 1

    # No Junction found.
    if consensus_lft_junction < 1 and consensus_rt_junction < 1:
        summary_data[6][0] += 1
        return [], summary_data

    hr_label = ""
    # If requested, do a search for HR Donor
    if hr_donor:
        rt_position = len(hr_donor)+25
        lft_position = 25
        donor_found = False
        while rt_position < len(consensus)-25:
            query_window = consensus[lft_position:rt_position]
            if query_window == hr_donor and not donor_found:
                summary_data[10][0]+=1
                donor_found = True
            elif query_window == hr_donor and donor_found:
                summary_data[10][1]+=1
            rt_position+=1
            lft_position+=1

        if donor_found:
            hr_label = "HR"
        # If HR Donor is found then find but do not score INDELS
        # ToDo: Is there a reason to keep the INDEL search for HR assays?
        bypass = False
        # if donor_found:
        if bypass:
            # extract the insertion
            consensus_insertion = ""
            if 0 < consensus_lft_junction < consensus_rt_junction:
                consensus_insertion = consensus[consensus_lft_junction:consensus_rt_junction]

            consensus_microhomology = ""
            if consensus_lft_junction > consensus_rt_junction > 0:
                consensus_microhomology = consensus[consensus_rt_junction:consensus_lft_junction]

            return [ldel, rdel, consensus_insertion, consensus_microhomology, consensus, consensus_lft_junction, consensus_rt_junction, target_lft_junction, target_rt_junction, "HR"], summary_data

    # extract the insertion
    consensus_insertion = ""
    if 0 < consensus_lft_junction < consensus_rt_junction:
        consensus_insertion = consensus[consensus_lft_junction:consensus_rt_junction]

        # If there is an N in the insertion then don't include read in the analysis.
        if "N" in consensus_insertion:
            return [], summary_data

        cut_found = True
        # Count number of insertions
        summary_data[4] += 1

    # Count left deletions
    if target_lft_junction < cutsite:
        cut_found = True
        summary_data[2] += 1

    # Count right deletions
    if target_rt_junction > cutsite:
        cut_found = True
        summary_data[3] += 1

    # extract the microhomology
    consensus_microhomology = ""
    if consensus_lft_junction > consensus_rt_junction > 0:
        consensus_microhomology = consensus[consensus_rt_junction:consensus_lft_junction]
        if consensus_microhomology:
            cut_found = True
            summary_data[5] += 1

    # No Cut found.
    if not cut_found:
        summary_data[6][1] += 1
        return [], summary_data


    return [ldel, rdel, consensus_insertion, consensus_microhomology, consensus, consensus_lft_junction, consensus_rt_junction, target_lft_junction, target_rt_junction, hr_label], summary_data