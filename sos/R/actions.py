#!/usr/bin/env python3
#
# This file is part of Script of Scripts (sos), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/bpeng2000/SOS
#
# Copyright (C) 2016 Bo Peng (bpeng@mdanderson.org)
##
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#

import os
import sys
import subprocess
import glob
import shutil
from collections.abc import Sequence

from sos.actions import SoS_Action, SoS_ExecuteScript
from sos.utils import env, natural_keys
from sos.monitor import ProcessMonitor, summarizeExecution


@SoS_Action(run_mode=['prepare', 'run', 'interactive'])
def R(script, args='', **kwargs):
    # > getOption('defaultPackages')
    # [1] "datasets"  "utils"     "grDevices" "graphics"  "stats"     "methods"
    return SoS_ExecuteScript(
        script, 'Rscript --default-packages=datasets,methods,utils,stats,grDevices,graphics ', '.R', args).run(**kwargs)

@SoS_Action(run_mode=['run', 'interactive'])
def Rmarkdown(script=None, output_file=None, **kwargs):
    '''This action can be used in three ways

    Rmarkdown:   output_file='report.html'
      script

    Rmarkdown(filename='report.sos', output_file='report.html')

    Rmarkdown(output_file='report.html')

    '''
    # if in prepare mode, check for Rmarkdown command
    #
    # in run mode, collect report and call Rmarkdown
    sos_script = env.sos_dict['__step_context__'].filename
    # this is the case for stirng input (test only)
    if sos_script is None:
        sos_script = 'string_input'
    if script is not None:
        # get a temporary file with the content
        script_file = '{}.Rmd'.format(os.path.basename(sos_script))
        with open(script_file, 'w') as report:
            report.write(script)
    elif 'filename' in kwargs:
        script_file = kwargs['filename']
    elif env.run_mode == 'interactive' and '__summary_report__' in env.sos_dict:
        script_file = env.sos_dict['__summary_report__']
    else:
        step_reports = glob.glob(os.path.join('.sos', 'report', '*'))
        step_reports.sort(key=natural_keys)
        # merge the files
        script_file = '{}.Rmd'.format(os.path.basename(sos_script))
        env.logger.trace('Gathering reports {} to {}'.format(', '.join(step_reports), script_file))
        with open(script_file, 'w') as combined:
            for step_report in step_reports:
                with open(step_report, 'r') as md:
                    combined.write(md.read())
    #
    arg_output_format = ', output_format={}'.format(kwargs['output_format']) if 'output_format' in kwargs else ''
    if output_file is None:
        raise RuntimeError('Parameter output_file is required for action Rmarkdown')
    elif not isinstance(output_file, str):
        raise RuntimeError('A filename is expected, {} provided'.format(output_file))
    extra_args = ''
    if 'extra_args' in kwargs:
        ea = kwargs['extra_args']
        if isinstance(ea, str):
            extra_args = ea
        elif isinstance(ea, Sequence):
            extra_args = ' '.join(list(ea))
        elif isinstance(ea, dict):
            extra_args = ' '.join('{}={:r}'.format(k,v) for k,v in ea.items())
        extra_args = ', ' + extra_args
    #
    #   render(input, output_format = NULL, output_file = NULL, output_dir = NULL,
    #        output_options = NULL, intermediates_dir = NULL,
    #        runtime = c("auto", "static", "shiny"),
    #        clean = TRUE, params = NULL, knit_meta = NULL, envir = parent.frame(),
    #        run_pandoc = TRUE, quiet = FALSE, encoding = getOption("encoding"))
    command = '''Rscript -e "rmarkdown::render({{}}, output_file={!r} {} {})" '''.format(
        os.path.abspath(output_file), arg_output_format, extra_args)
    try:
        cmd = command.replace('{}', '{!r}'.format(script_file))
        if env.run_mode == 'interactive':
            # need to catch output and send to python output, which will in trun be hijacked by SoS notebook
            p = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            pid = p.pid
            env.register_process(p.pid, 'Runing {}'.format(script_file))
            out, err = p.communicate()
            sys.stdout.write(out.decode())
            sys.stderr.write(err.decode())
            ret = p.returncode
        else:
            p = subprocess.Popen(cmd, shell=True)
            pid = p.pid
            if '__step_sig__' in env.sos_dict and env.sos_dict['__step_sig__'] is not None:
                m = ProcessMonitor(pid, msg=script, sig=env.sos_dict['__step_sig__'])
                m.start()
            env.register_process(pid, 'Runing {}'.format(script_file))
            ret = p.wait()
            if '__step_sig__' in env.sos_dict and env.sos_dict['__step_sig__'] is not None:
                summarizeExecution(pid)
    except Exception as e:
        env.logger.error(e)
    finally:
        env.deregister_process(p.pid)
        # os.remove(script_file)
    if ret != 0:
        temp_file = os.path.join('.sos', 'R_{}.Rmd'.format(os.getpid()))
        shutil.copyfile(script_file, temp_file)
        cmd = command.replace('{}', '{!r}'.format(temp_file))
        raise RuntimeError('Failed to execute script. The script is saved to {}. Please use command "{}" to test it.'
            .format(temp_file, cmd))
    env.logger.info('Report saved to {}'.format(output_file))

