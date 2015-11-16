#!/usr/bin/env python
"""Parses the output from the HyPhy script written by JCF to perform the
likelihood ratio test. Currently, only does a very crude method to pull out
the prediction line."""

import sys
import re
import os

def dir_contents(pred_dir):
    """Returns the contents of the directory, being sure to give each file's
    full path."""
    #   First we expand the tilde, if it was given
    userpath = os.path.expanduser(pred_dir)
    #   Then, remove any ../ or whatever that may be given
    normpath = os.path.normpath(userpath)
    #   And give the full path
    fullpath = os.path.abspath(normpath)
    #   Then get the directory contents
    raw_conts = os.listdir(fullpath)
    #   And convert them to full paths
    conts = [os.path.join(fullpath, x) for x in raw_conts]
    return conts


def get_prediction(fname):
    """Step through the prediction filename and return prediction values. This
    function uses a super simple heuristic to try to find the prediction line.
    It looks for an alignment report line that is different than the surrounding
    lines."""
    with open(fname, 'r') as f:
        for line in f:
            if re.match('^[0-9]', line):
                if len(line.split()) > 4:
                    return line.strip().split()


def get_gene_id(fname):
    """Use the given file name to get the gene ID. This, of course, assumes that
    the gene ID is the filename."""
    #   get the filename without the path
    nodirname = os.path.split(fname)[1]
    #   And then take the first part, without the extension
    gname = nodirname.split('.')[0]
    #   This is for soybena only
    gname = gname.replace('_', ':')
    return gname


def main(pred_dir):
    """Main function"""
    pred_files = dir_contents(pred_dir)
    for x in pred_files:
        #   Check if the file is empty
        if os.stat(x).st_size == 0:
            continue
        gname = get_gene_id(x)
        preds = get_prediction(x)
        print x
        print '\t'.join([gname] + preds)
    return


main(sys.argv[1])

