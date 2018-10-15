#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.


import subprocess
import sys
import re

from collections import Iterable, Mapping, Sequence
from typing import Any, Dict, Optional

from .eval import SoS_eval, SoS_exec, accessed_vars
from .parser import SoS_Step
from .targets import (dynamic, remote, sos_targets, sos_step)
from .utils import env, get_traceback, separate_options
from .executor_utils import  __null_func__

def analyze_section(section: SoS_Step, default_input: Optional[sos_targets] = None) -> Dict[str, Any]:
    '''Analyze a section for how it uses input and output, what variables
    it uses, and input, output, etc.'''
    from ._version import __version__

    # these are the information we need to build a DAG, by default
    # input and output and undetermined, and there are no variables.
    #
    # step input and output can be true "Undetermined", namely unspecified,
    # can be dynamic and has to be determined at run time, or undetermined
    # at this stage because something cannot be determined now.
    step_input: sos_targets = sos_targets()
    step_output: sos_targets = sos_targets()
    step_depends: sos_targets = sos_targets([])
    environ_vars = set()
    signature_vars = set()
    changed_vars = set()
    #
    # 1. execute global definition to get a basic environment
    #
    # FIXME: this could be made much more efficient
    if 'provides' in section.options:
        if '__default_output__' in env.sos_dict:
            step_output = env.sos_dict['__default_output__']
    else:
        # env.sos_dict = WorkflowDict()
        env.sos_dict.set('__null_func__', __null_func__)
        # initial values
        env.sos_dict.set('SOS_VERSION', __version__)
        SoS_exec('import os, sys, glob', None)
        SoS_exec('from sos.runtime import *', None)

    env.sos_dict.set('step_name', section.step_name())
    env.logger.trace(
        f'Analyzing {section.step_name()} with step_output {step_output}')

    #
    # Here we need to get "contant" values from the global section
    # Because parameters are considered variable, they has to be
    # removed. We achieve this by removing function sos_handle_parameter_
    # from the SoS_dict namespace
    #
    if section.global_def:
        try:
            SoS_exec('del sos_handle_parameter_\n' + section.global_def)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(e.stderr)
        except RuntimeError as e:
            if env.verbosity > 2:
                sys.stderr.write(get_traceback())
            raise RuntimeError(
                f'Failed to execute statements\n"{section.global_def}"\n{e}')
        finally:
            SoS_exec('from sos.runtime import sos_handle_parameter_', None)

    #
    # 2. look for input statement
    if 'shared' in section.options:
        svars = section.options['shared']
        if isinstance(svars, str):
            changed_vars.add(svars)
            svars = {svars: svars}
        elif isinstance(svars, Sequence):
            for item in svars:
                if isinstance(item, str):
                    changed_vars.add(item)
                elif isinstance(item, Mapping):
                    changed_vars |= set(item.keys())
                else:
                    raise ValueError(
                        f'Option shared should be a string, a mapping of expression, or list of string or mappings. {svars} provided')
        elif isinstance(svars, Mapping):
            changed_vars |= set(svars.keys())
        else:
            raise ValueError(
                f'Option shared should be a string, a mapping of expression, or list of string or mappings. {svars} provided')

    # look for input statement.
    input_statement_idx = [idx for idx, x in enumerate(
        section.statements) if x[0] == ':' and x[1] == 'input']
    if not input_statement_idx:
        input_statement_idx = None
    elif len(input_statement_idx) == 1:
        input_statement_idx = input_statement_idx[0]
    else:
        raise RuntimeError(
            f'More than one step input are specified in step {section.name if section.index is None else f"{section.name}_{section.index}"}')

    # if there is an input statement, analyze the statements before it, and then the input statement
    if input_statement_idx is not None:
        # execute before input stuff
        for statement in section.statements[:input_statement_idx]:
            if statement[0] == ':':
                if statement[1] == 'depends':
                    environ_vars |= accessed_vars(statement[2])
                    key, value = statement[1:3]
                    try:
                        args, kwargs = SoS_eval(f'__null_func__({value})')
                        if any(isinstance(x, (dynamic, remote)) for x in args):
                            step_depends = sos_targets()
                        else:
                            step_depends = sos_targets(*args)
                    except Exception as e:
                        env.logger.debug(
                            f"Args {value} cannot be determined: {e}")
                else:
                    raise RuntimeError(
                        f'Step input should be specified before {statement[1]}')
            else:
                environ_vars |= accessed_vars(statement[1])
        #
        # input statement
        stmt = section.statements[input_statement_idx][2]
        try:
            environ_vars |= accessed_vars(stmt)
            args, kwargs = SoS_eval(f'__null_func__({stmt})')

            if not args:
                if default_input is None:
                    step_input = sos_targets()
                else:
                    step_input = default_input
            elif not any(isinstance(x, (dynamic, remote)) for x in args):
                step_input = sos_targets(*args)
            env.sos_dict.set('input', step_input)

            if 'paired_with' in kwargs:
                pw = kwargs['paired_with']
                if pw is None or not pw:
                    pass
                elif isinstance(pw, str):
                    environ_vars.add(pw)
                elif isinstance(pw, Iterable):
                    environ_vars |= set(pw)
                elif isinstance(pw, Iterable):
                    # value supplied, no environ var
                    environ_vars |= set()
                else:
                    raise ValueError(
                        f'Unacceptable value for parameter paired_with: {pw}')
            if 'group_with' in kwargs:
                pw = kwargs['group_with']
                if pw is None or not pw:
                    pass
                elif isinstance(pw, str):
                    environ_vars.add(pw)
                elif isinstance(pw, Iterable):
                    environ_vars |= set(pw)
                elif isinstance(pw, Iterable):
                    # value supplied, no environ var
                    environ_vars |= set()
                else:
                    raise ValueError(
                        f'Unacceptable value for parameter group_with: {pw}')
            if 'for_each' in kwargs:
                fe = kwargs['for_each']
                if fe is None or not fe:
                    pass
                elif isinstance(fe, str):
                    environ_vars |= set([x.strip() for x in fe.split(',')])
                elif isinstance(fe, Sequence):
                    for fei in fe:
                        environ_vars |= set([x.strip()
                                             for x in fei.split(',')])
                else:
                    raise ValueError(
                        f'Unacceptable value for parameter fe: {fe}')

            if 'from_steps' in kwargs:
                if isinstance(kwargs['from_steps'], str):
                    step_depends.extend(sos_step(kwargs['from_steps']))
                elif isinstance(kwargs['from_steps'], Sequence):
                    step_depends.extend([sos_step(x) for x in kwargs['from_steps']])
                else:
                    raise ValueError(f'Invalid value for input option from {kwargs["from"]}')
        except Exception as e:
            # from parameter is very important so we will need to figure it out even
            # if the statement cannot be parsed.
            if 'from_steps' in stmt:
                for piece in separate_options(stmt):
                    mo = re.match('\s*from_steps\s*=\s*(?P<value>.*)', piece)
                    if mo:
                        try:
                            opt_values = SoS_eval(mo.group('value'))
                        except Exception as e:
                            raise ValueError(f'Invalid value for input option from {piece}: {e}')
                        if isinstance(opt_values, str):
                            step_depends.extend(sos_step(opt_values))
                        elif isinstance(opt_values, Sequence):
                            step_depends.extend([sos_step(x) for x in opt_values])
                        else:
                            raise ValueError(f'Invalid value for input option from {opt_values}')
            # if anything is not evalutable, keep Undetermined
            env.logger.debug(
                f'Input of step {section.name if section.index is None else f"{section.name}_{section.index}"} is set to Undertermined: {e}')
            # expression ...
            step_input = sos_targets(undetermined=stmt)
        input_statement_idx += 1
    else:
        # assuming everything starts from 0 is after input
        input_statement_idx = 0

    # other variables
    for statement in section.statements[input_statement_idx:]:
        # if input is undertermined, we can only process output:
        if statement[0] == '=':
            signature_vars |= accessed_vars('='.join(statement[1:3]))
        elif statement[0] == ':':
            key, value = statement[1:3]
            # if key == 'depends':
            environ_vars |= accessed_vars(value)
            # output, depends, and process can be processed multiple times
            try:
                args, kwargs = SoS_eval(f'__null_func__({value})')
                if not any(isinstance(x, (dynamic, remote)) for x in args):
                    if key == 'output':
                        step_output = sos_targets(*args)
                    elif key == 'depends':
                        step_depends.extend(sos_targets(*args))
            except Exception as e:
                env.logger.debug(f"Args {value} cannot be determined: {e}")
        else:  # statement
            signature_vars |= accessed_vars(statement[1])
    # finally, tasks..
    if section.task:
        signature_vars |= accessed_vars(section.task)
    if 'provides' in section.options and '__default_output__' in env.sos_dict and step_output.valid():
        for out in env.sos_dict['__default_output__']:
            # 981
            if not isinstance(out, sos_step) and out not in step_output:
                raise ValueError(
                    f'Defined output fail to produce expected output: {step_output} generated, {env.sos_dict["__default_output__"]} expected.')

    return {
        'step_name': f'{section.name}_{section.index}' if isinstance(section.index, int) else section.name,
        'step_input': step_input,
        'step_output': step_output,
        'step_depends': step_depends,
        # variables starting with __ are internals...
        'environ_vars': {x for x in environ_vars if not x.startswith('__')},
        'signature_vars': {x for x in signature_vars if not x.startswith('__')},
        'changed_vars': changed_vars
    }
