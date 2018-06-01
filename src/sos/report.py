#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.
from .utils import env


def render_report(output_file, workflow_info):
    with open(workflow_info, 'r') as wi:
        info = wi.read()

    with open(output_file, 'w') as wo:
        wo.write(info)

    env.logger.info(f'Summary of workflow saved to {output_file}')
