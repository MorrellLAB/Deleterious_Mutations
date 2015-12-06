#!/usr/bin/env python
"""Use the output from bedtools genomecov to identify regions that have
sufficient coverage to be considered sequenced regions. This is useful for
identifying regions post-hoc that have been part of a sequence capture
experiment. This script will prodice a BED file, and assumes that all positions
reported are 0-based."""

import sys
#   We set a minimum number of reads across all our samples here to be
#   considered "captured"
coverage_cutoff = 20
#   We also set a tolerance for collapsing features. If two intervals of
#   sufficient coverage are closer than this distance, we will collapse them
#   into one contiguous region.
tolerance = 50
#   This is the minimum length a feature should be
min_dist = 200


def consecutive(x):
    """A function to return start and end as tuples for consecutive integers
    in a list."""
    #   Sometimes x is a single element
    if len(x) == 1:
        return [(x[0], x[0]+1)]
    starts = []
    ends = []
    for i in range(1, len(x)):
        #   Append the first index to starts
        if i == 1:
            starts.append(x[i-1])
        #   If the difference is 1, then they are consecutive numbers and
        #   are part of the same range.
        if x[i] - x[i-1] == 1:
            continue
        else:
            ends.append(x[i-1])
            starts.append(x[i])
    #   Append the final value to the ends
    else:
        ends.append(x[i])
    return zip(starts, ends)


def collapse(intervals, tolerance=tolerance, distance=min_dist):
    """A function to combine intervals that are close together. If they have
    'tolerance' basepairs or less between them, they will be combined."""
    #   If there is only one interval, we return it
    if len(intervals) == 1:
        if intervals[0][1] - intervals[0][0] > distance:
            return intervals
        else:
            return None
    #   We take a similar approach as in consecutive(). Compare the end of one
    #   interval with the start of the next.
    starts = []
    ends = []
    for i in range(1, len(intervals)):
        if i == 1:
            starts.append(intervals[i-1][0])
        #   Check the start of the current interval to see if it is less than
        #   'tolerance' basepairs from the end of the previous
        if intervals[i][0] - intervals[i-1][1] <= tolerance:
            continue
        else:
            ends.append(intervals[i-1][1])
            starts.append(intervals[i][0])
    else:
        #   We have to append the final endpoint, no matter what
        ends.append(intervals[i][1])
    col = []
    for i in zip(starts, ends):
        if i[1] - i[0] > distance:
            col.append(i)
    return col

captured_contigs = {}
with open(sys.argv[1], 'r') as f:
    for line in f:
        tmp = line.strip().split()
        #   Check the number of reads that are mapped to this location. If it
        #   is below the threshold, then we continue
        try:
            cov = int(tmp[-1])
        except ValueError:
            sys.stderr.write('Warning, encountered a non-integer value!\n')
            cov = 0
        if cov < coverage_cutoff:
            continue
        else:
            contig = tmp[0]
            if contig in captured_contigs:
                captured_contigs[contig].append(int(tmp[1]))
            else:
                captured_contigs[contig] = [int(tmp[1])]

#   Now, we print out the BED file
for ctg in captured_contigs.keys():
    consec = consecutive(captured_contigs[ctg])
    collapsed = collapse(consec)
    if not collapsed:
        continue
    for interval in collapsed:
        #   We have to subtract 1, since the coordinates are presented as
        #   1-based positions.
        print '\t'.join([ctg, str(interval[0]-1), str(interval[1]-1)])
