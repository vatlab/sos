#!/usr/bin/env python
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import pickle


def encode_msg(msg):
    return pickle.dumps(msg)


def decode_msg(data):
    return pickle.loads(data)
