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
from .targets import (dynamic, remote, sos_targets, sos_step, named_output)
from .utils import env, get_traceback, separate_options
from .executor_utils import  __null_func__, prepare_env, strip_param_defs
from .syntax import SOS_TARGETS_OPTIONS

def get_param_of_function(name, param_list, extra_dict={}):
    tree = ast.parse(f'__null_func__({param_list})')
    # x.func can be an attribute (e.g. a.b()) and do not have id
    funcs = [x for x in ast.walk(tree) if x.__class__.__name__ == 'Call' and hasattr(x.func, 'id') and x.func.id == name]
    params = []
    for func in funcs:
        for arg in func.args:
            try:
                params.append([ast.literal_eval(arg)])
            except Exception as e:
                params.append([eval(compile(ast.Expression(body=arg), filename='<string>', mode="eval"), extra_dict)])
        for kwarg in func.keywords:
            try:
                params.append([kwarg.arg, ast.literal_eval(kwarg.value)])
            except Exception as e:
                params.append([kwarg.arg, eval(compile(ast.Expression(body=kwarg.value), filename='<string>', mode="eval"), extra_dict)])
    return params

def get_names_of_param(name, param_list, extra_dict={}):
    tree = ast.parse(f'__null_func__({param_list})')
    kwargs = [x for x in ast.walk(tree) if x.__class__.__name__ == 'keyword' and x.arg == name]

    values = []
    for kwarg in kwargs:
        values.extend([x.s for x in ast.walk(kwarg.value) if x.__class__.__name__ == 'Str'])
    return values

def get_num_of_args_and_names_of_kwargs(name, param_list):
    tree = ast.parse(f'__null_func__({param_list})')
    res = []
    for func in [x for x in ast.walk(tree) if x.__class__.__name__ == 'Call' and hasattr(x.func, 'id') and x.func.id == name]:
        res.append([len(func.args), [x.arg for x in func.keywords]])
    return res

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

    # signature vars should contain parameters defined in global section
    # #1155
    signature_vars = set(section.parameters.keys() & accessed_vars(strip_param_defs(section.global_def)))

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
            step_depends.extend([sos_step(x) for x in get_output_from_steps(stmt, section.last_step)])
        if 'named_output' in stmt:
            # there can be multiple named_output calls
            pars = get_param_of_function('named_output', stmt,
                extra_dict=env.sos_dict._dict)
            for par in pars:
                # a single argument
                if len(par) == 1:
                    if not isinstance(par[0], str):
                        raise ValueError(f'Value for named_output can only be a name (str): {par[0]} provided')
                    step_depends.extend(named_output(par[0]))
                else:
                    if par[0] in SOS_TARGETS_OPTIONS:
                        continue
                    elif par[0] == 'name':
                        if not isinstance(par[1], str):
                            raise ValueError(f'Value for named_output can only be a name (str): {par[1]} provided')
                        step_depends.extend(named_output(par[1]))
                    else:
                        raise ValueError(f'Unacceptable keyword argument {par[0]} for named_output()')


    depends_idx = find_statement(section, 'depends')
    if depends_idx is not None:
        value = section.statements[depends_idx][2]
        try:
            args, kwargs = SoS_eval(f'__null_func__({value})',
                extra_dict=env.sos_dict._dict)
            if any(isinstance(x, (dynamic, remote)) for x in args):
                step_depends = sos_targets()
            else:
                step_depends.extend(sos_targets(*args))
        except Exception as e:
            pass
            # env.logger.debug(f"Args {value} cannot be determined: {e}")
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
        env.sos_dict._dict.update({
            'output_from': lambda *args, **kwargs: None,
            'named_output': lambda *args, **kwargs: None
            })
        args, kwargs = SoS_eval(f'__null_func__({stmt})',
            extra_dict=env.sos_dict._dict)
        env.sos_dict._dict.pop('output_from')
        env.sos_dict._dict.pop('named_output')
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

def get_step_output(section, default_output):
    '''determine step output'''
    step_output: sos_targets = sos_targets()
    #
    if 'provides' in section.options and default_output:
        step_output = default_output

    # look for input statement.
    output_idx = find_statement(section, 'output')
    if output_idx is None:
        return step_output

    # output statement
    value = section.statements[output_idx][2]
    # output, depends, and process can be processed multiple times
    try:
        args, kwargs = SoS_eval(f'__null_func__({value})',
            extra_dict=env.sos_dict._dict)
        if not any(isinstance(x, (dynamic, remote)) for x in args):
            step_output = sos_targets(*args)
    except Exception as e:
        pass
        # env.logger.debug(f"Args {value} cannot be determined: {e}")

    if 'provides' in section.options and default_output is not None and step_output.valid():
        for out in default_output:
            # 981
            if not isinstance(out, sos_step) and out not in step_output:
                raise ValueError(
                    f'Defined output fail to produce expected output: {step_output} generated, {default_output} expected.')
    return step_output

def get_output_from_steps(stmt, last_step):
    '''
    Extract output_from(1), output_from('step_1'), and output_from([1, 2])
    to determine dependent steps
    '''
    opt_values = get_param_of_function('output_from', stmt, extra_dict=env.sos_dict._dict)

    def step_name(val):
        if isinstance(val, str):
            return val
        elif isinstance(val, int):
            if val == -1:
                if last_step is None:
                    # there is a case where a regular step is checked as auxiliary step.
                    # we will postpone the decision later because the step might not be
                    # used as such
                    return None
                return last_step
            if '_' in env.sos_dict['step_name']:
                return f"{env.sos_dict['step_name'].rsplit('_',1)[0]}_{val}"
            else:
                return str(val)
        else:
            raise ValueError(f'Invalid value {val} for output_from() function')

    res = []
    for value in opt_values:
        if len(value) == 1:
            # regular argument
            value = value[0]
        elif value[0] == 'steps':
            value = value[1]
        elif value[0] in SOS_TARGETS_OPTIONS:
            continue
        else:
            raise ValueError(f'Unacceptable keyword argument {value[0]} for function output_from')
        if isinstance(value, (int, str)):
            res.append(step_name(value))
        elif isinstance(value, Sequence):
            res.extend([step_name(x) for x in value])
        else:
            raise ValueError(f'Invalid value for input option from {value}')
    return [x for x in res if x is not None]


analysis_cache = {}

def analyze_section(section: SoS_Step, default_input: Optional[sos_targets] = None,
    default_output: Optional[sos_targets] = None,
    vars_and_output_only: bool = False) -> Dict[str, Any]:
    '''Analyze a section for how it uses input and output, what variables
    it uses, and input, output, etc.'''
    from ._version import __version__

    analysis_key = (section.md5, section.step_name(),
        default_input.target_name() if hasattr(default_input, 'target_name') else '',
        default_output.target_name() if hasattr(default_output, 'target_name') else '', vars_and_output_only)
    if analysis_key in analysis_cache:
        return analysis_cache[analysis_key]

    # initialiaze environment, without handling parameter
    prepare_env()

    env.sos_dict.set('step_name', section.step_name())
    env.sos_dict.set('__null_func__', __null_func__)
    env.logger.trace(f'Analyzing {section.step_name()}')

    #
    # Here we need to get "contant" values from the global section
    # Because parameters are considered variable, they has to be
    # removed.
    #
    if section.global_def:
        try:
            SoS_exec(strip_param_defs(section.global_def))
        except subprocess.CalledProcessError as e:
            raise RuntimeError(e.stderr)
        except RuntimeError as e:
            if env.verbosity > 2:
                sys.stderr.write(get_traceback())
            raise RuntimeError(
                f'Failed to execute statements\n"{section.global_def}"\n{e}')

    res = {
        'step_name': section.step_name(),
        'step_output': get_step_output(section, default_output),
        # variables starting with __ are internals...
        'environ_vars': get_environ_vars(section),
        'signature_vars': get_signature_vars(section),
        'changed_vars': get_changed_vars(section)
    }
    if not vars_and_output_only:
        res['step_input'] = get_step_input(section, default_input)
        res['step_depends'] = get_step_depends(section)
    analysis_cache[analysis_key] = res
    return res
