# ENCODING: UTF-8

# This file was created by Kirill Vinnikov on August 10, 2020
# Copyright 2020 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root dircectory of the EXONtools package.

from re import split


def atoi(text):
    return int(text) if text.isdigit() else text


def natural_sort(text):
    return [atoi(c) for c in split('(\d+)', text)]


def mostcommon(lst):
    return max(set(lst), key=lst.count)
