# ----------------------------------------------------------------------------
#
# TITLE -
# AUTHOR - James Lane
# PROJECT -
# CONTENTS:
#	1.
#
# ----------------------------------------------------------------------------
#
# Docstrings and metadata:
'''Example'''
__author__ = "James Lane"


#Imports
import re

def sort_nicely(l):
    def tryint(s):
            try:
                    return int(s)
            except:
                    return s
    #def
    def alphanum_key(s):
            return [ tryint(c) for c in re.split('([0-9]+)', s) ]
    #def
    l.sort(key=alphanum_key)
#def
