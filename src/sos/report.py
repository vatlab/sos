#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import os
from .utils import env, TimeoutInterProcessLock

from contextlib import contextmanager

@contextmanager
def workflow_report(mode='a'):
    if '__workflow_sig__' not in env.sos_dict:
        raise RuntimeError('workflow_sig is not defined.')
    workflow_sig = env.sos_dict['__workflow_sig__']
    with TimeoutInterProcessLock(workflow_sig + '_'):
        with open(workflow_sig, mode) as sig:
            yield sig


def render_report(output_file, workflow_info):
    data = defaultdict(list)
    with open(workflow_info, 'r') as wi:
        for line in wi:
            try:
                entry_type, d = line.split(':', 1)
                entry_data = eval(d)
            except Exception as e:
                env.logger.debug(f'Failed to read report line {line}: {e}')
            data[entry_type].append(entry_data)

    with open(output_file, 'w') as wo:
        wo.write(info)

    env.logger.info(f'Summary of workflow saved to {output_file}')
