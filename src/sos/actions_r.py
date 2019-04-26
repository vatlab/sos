#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import os
import shutil
import subprocess
import sys
import tempfile

from sos.actions import SoS_Action, SoS_ExecuteScript, collect_input
from sos.eval import interpolate
from sos.targets import sos_targets
from sos.utils import env

from .targets_r import R_library


@SoS_Action(
    acceptable_args=['script', 'args'],
    default_args={
        'default_env': {
            'R_DEFAULT_PACKAGES':
                'datasets,methods,utils,stats,grDevices,graphics'
        }
    })
def R(script, args='', **kwargs):
    '''Execute specified script with command Rscript, with default options
    "--default-packages=datasets,methods,utils,stats,grDevices,graphics". This action accepts
    common action arguments such as input, active, workdir, docker_image and args.
    In particular, content of one or more files  specified by option input would be
    prepended before the specified script.
    '''
    # > getOption('defaultPackages')
    # [1] "datasets"  "utils"     "grDevices" "graphics"  "stats"     "methods"
    return SoS_ExecuteScript(script, 'Rscript', '.R', args).run(**kwargs)


@SoS_Action(acceptable_args=['script', 'args'])
def Rmarkdown(script=None,
              input=None,
              output=None,
              args='{input:r}, output_file={output:ar}',
              **kwargs):
    '''Convert input file to output using Rmarkdown

    The input can be specified in three ways:

    1. instant script, which is assumed to be in md format

    Rmarkdown:   output='report.html'
      script

    2. one or more input files. The format is determined by extension of input file

    Rmarkdown(input, output='report.html')

    3. input file specified by command line option `-r` .
    Rmarkdown(output='report.html')

    If no output is specified, it is assumed to be in html format
    and is written to standard output.

    You can specify more options using the args parameter of the action. The default value
    of args is `${input!r} --output ${output!ar}'
    '''
    if not R_library('rmarkdown').target_exists():
        raise RuntimeError('Library rmarkdown does not exist')

    input = sos_targets(collect_input(script, input))

    output = sos_targets(output)
    if len(output) == 0:
        write_to_stdout = True
        output = sos_targets(
            tempfile.NamedTemporaryFile(
                mode='w+t', suffix='.html', delete=False).name)
    else:
        write_to_stdout = False
    #
    ret = 1
    try:
        #   render(input, output_format = NULL, output_file = NULL, output_dir = NULL,
        #        output_options = NULL, intermediates_dir = NULL,
        #        runtime = c("auto", "static", "shiny"),
        #        clean = TRUE, params = NULL, knit_meta = NULL, envir = parent.frame(),
        #        run_Rmarkdown = TRUE, quiet = FALSE, encoding = getOption("encoding"))
        cmd = interpolate(f'Rscript -e "rmarkdown::render({args})"', {
            'input': input,
            'output': output
        })
        if 'ACTION' in env.config['SOS_DEBUG'] or 'ALL' in env.config['SOS_DEBUG']:
            env.log_to_file('ACTION', f'Running command "{cmd}"')
        if env.config['run_mode'] == 'interactive':
            # need to catch output and send to python output, which will in trun be hijacked by SoS notebook
            p = subprocess.Popen(
                cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            #pid = p.pid
            out, err = p.communicate()
            sys.stdout.write(out.decode())
            sys.stderr.write(err.decode())
            ret = p.returncode
        else:
            p = subprocess.Popen(cmd, shell=True)
            #pid = p.pid
            ret = p.wait()
    except Exception as e:
        env.logger.error(e)
    if ret != 0:
        temp_file = os.path.join('.sos', f'{"Rmarkdown"}_{os.getpid()}.md')
        shutil.copyfile(str(input), temp_file)
        cmd = interpolate(f'Rscript -e "rmarkdown::render({args})"', {
            'input': input,
            'output': sos_targets(temp_file)
        })
        raise RuntimeError(
            f'Failed to execute script. Please use command \n"{cmd}"\nunder {os.getcwd()} to test it.'
        )
    if write_to_stdout:
        with open(str(output[0])) as out:
            sys.stdout.write(out.read())
    else:
        env.logger.info(f'Report saved to {output}')
