#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import ast
import subprocess
import sys
import re

from collections import Iterable, Mapping, Sequence
from typing import Any, Dict, Optional

from .eval import SoS_eval, SoS_exec, accessed_vars
from .parser import SoS_Step
from .targets import (dynamic, remote, sos_targets, sos_step)
from .utils import env, get_traceback, separate_options
from .executor_utils import  __null_func__, __sos_groups__, __output_from__


def get_param_of_function(name, stmt, extra_dict={}):
    tree = ast.parse(stmt)
    funcs = [x for x in ast.walk(tree) if x.__class__.__name__ == 'Call' and x.func.id == name]
    params = []
    for func in funcs:
        for arg in func.args:
            try:
                params.append(ast.literal_eval(arg))
            except Exception as e:
                params.append(eval(compile(ast.Expression(body=arg), filename='<string>', mode="eval"), extra_dict))
        for kwarg in func.keywords:
            try:
                params.append(ast.literal_eval(kwarg.value))
            except Exception as e:
                params.append(eval(compile(ast.Expression(body=kwarg.value), filename='<string>', mode="eval"), extra_dict))
    return params

def get_names_of_param(name, param_list, extra_dict={}):
    tree = ast.parse(f'__null_func__({param_list})')
    kwargs = [x for x in ast.walk(tree) if x.__class__.__name__ == 'keyword' and x.arg == name]

    values = []
    for kwarg in kwargs:
        values.extend([x.s for x in ast.walk(kwarg.value) if x.__class__.__name__ == 'Str'])
    return values


def find_statement(section, name):
    stmt_idx = [idx for idx, x in enumerate(section.statements) if x[0] == ':' and x[1] == name]
    if not stmt_idx:
        return None
    elif len(stmt_idx) == 1:
        return stmt_idx[0]
    else:
        raise RuntimeError(
            f'More than one step {name} statement are specified in step {section.step_name()}')

def get_changed_vars(section: SoS_Step):
    '''changed vars are variables that are "shared" and therefore "provides"
    to others '''
    if 'shared' not in section.options:
        return set()

    changed_vars = set()
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
    return changed_vars

def get_environ_vars(section):
    # environ variables are variables that are needed in the step
    # and should be passed to the step, but it should not trigger
    # dependency. The reason why it is needed is because if an environ
    # var is "shared" by another step, then this step should be dependent
    # upon that step. However, this is only implement for "forward-steps"
    environ_vars = set()
    before_input = True
    for statement in section.statements:
        if statement[0] in ('!', '='):
            if before_input:
                environ_vars |= accessed_vars(statement[1])
            continue
        environ_vars |= accessed_vars(statement[2])
        # there is only nasty problem here. With the old paird_with etc
        # they accept parameter name, not parameter, so we will need to
        # get the value of the parameters
        if statement[1] != 'input':
            continue
        else:
            before_input = False
        if 'paired_with' in statement[2]:
            try:
                pws = get_names_of_param('paired_with', statement[2],  extra_dict=env.sos_dict._dict)
                environ_vars |= set(pws)
            except Exception as e:
                raise ValueError(f'Failed to parse parameter paired_with: {e}')
        if 'group_with' in statement[2]:
            try:
                pws = get_names_of_param('group_with', statement[2], extra_dict=env.sos_dict._dict)
                environ_vars |= set(pws)
            except Exception as e:
                raise ValueError(f'Failed to parse parameter group_with: {e}')
        if 'for_each' in statement[2]:
            try:
                pws = get_names_of_param('for_each', statement[2], extra_dict=env.sos_dict._dict)
                for pw in pws:
                    environ_vars |= set(pw.split(','))
            except Exception as e:
                raise ValueError(f'Failed to parse parameter for_each: {e}')
    return {x for x in environ_vars if not x.startswith('__')}

def get_signature_vars(section):
    '''Get signature variables which are variables that will be
    saved with step signatures'''
    signature_vars = set()

    input_idx = find_statement(section, 'input')
    after_input_idx = 0 if input_idx is None else input_idx + 1

    for statement in section.statements[after_input_idx:]:
        if statement[0] == '=':
            signature_vars |= accessed_vars('='.join(statement[1:3]))
        elif statement[0] == '!':
            signature_vars |= accessed_vars(statement[1])
    # finally, tasks..
    if section.task:
        signature_vars |= accessed_vars(section.task)
    return {x for x in signature_vars if not x.startswith('__')}

def get_step_depends(section):
    step_depends: sos_targets = sos_targets([])

    input_idx = find_statement(section, 'input')
    if input_idx is not None:
        # input statement
        stmt = section.statements[input_idx][2]
        if 'output_from' in stmt:
            step_depends.extend([sos_step(x) for x in get_output_from_steps(stmt)])

    depends_idx = find_statement(section, 'depends')
    if depends_idx is not None:
        value = section.statements[depends_idx][2]
        try:
            args, kwargs = SoS_eval(f'__null_func__({value})',
                extra_dict={
                    '__null_func__': __null_func__
                    }
            )
            if any(isinstance(x, (dynamic, remote)) for x in args):
                step_depends = sos_targets()
            else:
                step_depends.extend(sos_targets(*args))
        except Exception as e:
            env.logger.debug(f"Args {value} cannot be determined: {e}")
    return step_depends

def get_step_input(section, default_input):
    '''Find step input
    '''
    step_input: sos_targets = sos_targets()
    # look for input statement.
    input_idx = find_statement(section, 'input')
    if input_idx is None:
        return step_input

    # input statement
    stmt = section.statements[input_idx][2]
    try:
        args, kwargs = SoS_eval(f'__null_func__({stmt})',
            extra_dict={
                '__null_func__': __null_func__,
                'sos_groups': __sos_groups__,
                'output_from': __output_from__
                })
        if not args:
            if default_input is None:
                step_input = sos_targets()
            else:
                step_input = default_input
        elif not any(isinstance(x, (dynamic, remote)) for x in args):
            step_input = sos_targets(*args)
    except Exception as e:
        # if anything is not evalutable, keep Undetermined
        env.logger.debug(
            f'Input of step {section.name if section.index is None else f"{section.name}_{section.index}"} is set to Undertermined: {e}')
        # expression ...
        step_input = sos_targets(_undetermined=stmt)
    return step_input

def get_step_output(section):
    '''determine step output'''
    step_output: sos_targets = sos_targets()
    #
    if 'provides' in section.options:
        if '__default_output__' in env.sos_dict:
            step_output = env.sos_dict['__default_output__']

    # look for input statement.
    output_idx = find_statement(section, 'output')
    if output_idx is None:
        return step_output

    # output statement
    value = section.statements[output_idx][2]
    # output, depends, and process can be processed multiple times
    try:
        args, kwargs = SoS_eval(f'__null_func__({value})',
            extra_dict={
                '__null_func__': __null_func__,
                'sos_groups': __sos_groups__,
                'output_from': __output_from__
                })
        if not any(isinstance(x, (dynamic, remote)) for x in args):
            step_output = sos_targets(*args)
    except Exception as e:
        env.logger.debug(f"Args {value} cannot be determined: {e}")

    if 'provides' in section.options and '__default_output__' in env.sos_dict and step_output.valid():
        for out in env.sos_dict['__default_output__']:
            # 981
            if not isinstance(out, sos_step) and out not in step_output:
                raise ValueError(
                    f'Defined output fail to produce expected output: {step_output} generated, {env.sos_dict["__default_output__"]} expected.')
    return step_output

def get_output_from_steps(stmt):
    '''
    Extract output_from(1), output_from('step_1'), and output_from([1, 2])
    to determine dependent steps
    '''
    opt_values = get_param_of_function('output_from', f'null_func({stmt})', extra_dict=env.sos_dict._dict)

    def step_name(val):
        if isinstance(val, str):
            return val
        elif isinstance(val, int):
            if '_' in env.sos_dict['step_name']:
                return f"{env.sos_dict['step_name'].rsplit('_',1)[0]}_{val}"
            else:
                return str(val)
        else:
            raise ValueError(f'Invalid value {val} for output_from() function')

    res = []
    for value in opt_values:
        if isinstance(value, (int, str)):
            res.append(step_name(value))
        elif isinstance(value, Sequence):
            res.extend([step_name(x) for x in value])
        else:
            raise ValueError(f'Invalid value for input option from {value}')
    return res


def analyze_section(section: SoS_Step, default_input: Optional[sos_targets] = None) -> Dict[str, Any]:
    '''Analyze a section for how it uses input and output, what variables
    it uses, and input, output, etc.'''
    from ._version import __version__

    # initial values
    env.sos_dict.set('SOS_VERSION', __version__)
    SoS_exec('import os, sys, glob', None)
    SoS_exec('from sos.runtime import *', None)

    env.sos_dict.set('step_name', section.step_name())
    env.logger.trace(f'Analyzing {section.step_name()}')

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

    return {
        'step_name': section.step_name(),
        'step_input': get_step_input(section, default_input),
        'step_output': get_step_output(section),
        'step_depends': get_step_depends(section),
        # variables starting with __ are internals...
        'environ_vars': get_environ_vars(section),
        'signature_vars': get_signature_vars(section),
        'changed_vars': get_changed_vars(section)
    }
