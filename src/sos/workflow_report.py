#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import getpass
import os
import time
import base64
from collections import defaultdict

from .utils import env, format_duration, dot_to_gif
from .controller import request_answer_from_controller
from ._version import __version__


class WorkflowSig(object):

    def __init__(self, workflow_id):
        self.data = defaultdict(lambda: defaultdict(list))

        for entry_type, id, item in request_answer_from_controller(
            ['workflow_sig', 'records', workflow_id]):
            try:
                self.data[entry_type][id].append(item.strip())
            except Exception as e:
                env.logger.debug(f'Failed to read report line: {e}')

    def convert_time(self, info):
        for key in [x for x in info.keys() if x.endswith('_time')]:
            info[key + '_str'] = time.strftime('%Y-%m-%d %H:%M:%S',
                                               time.localtime(info[key]))
        if 'start_time' in info and 'end_time' in info:
            info['duration'] = int(info['end_time'] - info['start_time'])
            info['duration_str'] = format_duration(info['duration'])
        return info

    def workflows(self):
        try:
            # workflows has format
            # workflow  workflow_id dict1
            # workflow  workflow_id dict2
            # workflow  workflow_id2 dict1
            workflows = defaultdict(dict)
            for id, values in self.data['workflow'].items():
                for val in values:
                    workflows[id].update(eval(val))
            for v in workflows.values():
                self.convert_time(v)
                if 'dag' in v:
                    try:
                        v['dag_img'] = dot_to_gif(
                            v['dag'], warn=env.logger.warning)
                    except Exception as e:
                        env.logger.warning(
                            f'Failed to obtain convert dag to image: {e}')
                if 'script' in v:
                    v['script'] = base64.b64decode(v['script']).decode()
            return workflows
        except Exception as e:
            env.logger.warning(f'Failed to obtain workflow information: {e}')
            return {}

    def tasks(self):

        def merge_dict(items):
            ret = {}
            for item in items:
                ret.update(eval(item))
            return ret

        try:
            # there can be multiple task status for each id
            tasks = {
                id: merge_dict(res) for id, res in self.data['task'].items()
            }
            for val in tasks.values():
                self.convert_time(val)
                if 'peak_cpu' in val:
                    val['peak_cpu_str'] = f'{val["peak_cpu"]:.1f}%'
                if 'peak_mem' in val:
                    val['peak_mem_str'] = f'{val["peak_mem"] / 1024 / 1024 :.1f}Mb'
            return tasks
        except:
            return {}

    def steps(self):
        try:
            return {
                wf: [self.convert_time(eval(x)) for x in steps
                    ] for wf, steps in self.data['step'].items()
            }
        except Exception as e:
            env.logger.warning(e)
            return {}

    def transcripts(self):
        try:
            return {
                step: [self.convert_time(eval(x)) for x in items
                      ] for step, items in self.data['transcript'].items()
            }
        except Exception as e:
            env.logger.warning(e)
            return {}

    def placeholders(self):
        try:
            return self.data['placeholder']['file_target']
        except:
            return []


def calc_timeline(info, start_time, total_duration):
    if total_duration == 0 or 'start_time' not in info or 'end_time' not in info:
        info['before_percent'] = 0
        info['during_percent'] = 100
        info['after_precent'] = 0
        return
    info['before_percent'] = int(
        (info['start_time'] - start_time) * 100 / total_duration)
    info['during_percent'] = max(1,
                                 int(info['duration'] * 100 / total_duration))
    info['after_percent'] = 100 - \
        info['before_percent'] - info['during_percent']
    return info


def render_report(output_file, workflow_id):
    data = WorkflowSig(workflow_id)

    from jinja2 import Environment, PackageLoader, select_autoescape
    environment = Environment(
        loader=PackageLoader('sos', 'templates'),
        autoescape=select_autoescape(['html', 'xml']))
    environment.filters['basename'] = os.path.basename
    template = environment.get_template('workflow_report.tpl')

    context = {
        'workflows': data.workflows(),
        'tasks': data.tasks(),
        'steps': data.steps(),
        'transcripts': data.transcripts(),
        'sos_version': __version__,
        'user': getpass.getuser(),
        'time_now_str': time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()),
    }
    # derived context
    context['master_id'] = next(iter(
        context['workflows'].values()))['master_id']
    try:
        # calculate percentage
        start_time = context['workflows'][context['master_id']]['start_time']
        total_duration = context['workflows'][
            context['master_id']]['end_time'] - start_time
        for info in context['workflows'].values():
            calc_timeline(info, start_time, total_duration)
        for steps in context['steps'].values():
            for step in steps:
                calc_timeline(step, start_time, total_duration)
        for info in context['tasks'].values():
            calc_timeline(info, start_time, total_duration)
        with open(output_file, 'w') as wo:
            wo.write(template.render(context))
        env.logger.info(f'Summary of workflow saved to {output_file}')
    except Exception as e:
        env.logger.error(f'Failed to generate report {output_file}: {e}')
