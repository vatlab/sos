#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

# A complete rewritten will be needed after spyder officially supports
# third-party kernel.
#
# Note that this kernel is only used by Spyder, not by jupyter notebook
# and qtconsole.
#
import os
import sys
import re
import argparse
from sos.utils import env
from sos.eval import SoS_exec

from .kernel import SoS_Kernel
from spyder.utils.ipython.spyder_kernel import SpyderKernel

class SoS_SpyderKernel(SoS_Kernel, SpyderKernel):
    """Spyder kernel for Jupyter"""

    MAGIC_EDIT = re.compile('^%edit(\s|$)')

    def __init__(self, *args, **kwargs):
        #super(SoS_SpyderKernel, self).__init__(*args, **kwargs)
        SpyderKernel.__init__(self, *args, **kwargs)
        SoS_Kernel.__init__(self, *args, **kwargs)

        # supposedly this should be set by namespacebrowser.py when the browser
        # window starts. no idea why this does not work.
        self.namespace_view_settings = {'check_all': False,
            'exclude_private': True, 'remote_editing': False, 'autorefresh': False,
            'exclude_capitalized': False, 'exclude_uppercase': True,
            'excluded_names': ['nan', 'inf', 'infty', 'little_endian', \
                'colorbar_doc', 'typecodes', '__builtins__', '__main__', '__doc__',\
                'NaN', 'Inf', 'Infinity', 'sctypes', 'rcParams', 'rcParamsDefault', \
                'sctypeNA', 'typeNA', 'False_', 'True_', 'run_mode', 'step_name'] + \
                list(self.original_keys if self.original_keys else []),
            'exclude_unsupported': True, 'minmax': False}
        #
        self.shell.user_ns = env.sos_dict._dict

    def send_sos_msg(self, msg):
        pass

    def get_edit_parser(self):
        parser = argparse.ArgumentParser(prog='%edit',
            description='Edit an existing file in spyder')
        parser.add_argument('filenames', nargs='+')
        parser.add_argument('-c', '--cd', action='store_true', dest='__switch_dir__')
        parser.error = self._parse_error
        return parser

    def handle_magic_edit(self, options):
        import subprocess
        import shlex
        options = self._interpolate_option(options)
        if options is None:
            return
        parser = self.get_edit_parser()
        args = parser.parse_args(shlex.split(options))
        args.filenames = [os.path.expanduser(x) for x in args.filenames]
        for filename in args.filenames:
            if not os.path.isfile(filename):
                self.warn('File does not exist: {}'.format(filename))
                return
        import1 = "import sys"
        import2 = "from spyder.app.start import send_args_to_spyder"
        code = "send_args_to_spyder([{}])".format(','.join('"{}"'.format(x) for x in args.filenames))
        cmd = "{0} -c '{1}; {2}; {3}'".format(sys.executable,
            import1, import2, code)
        subprocess.call(cmd, shell=True)
        if args.__switch_dir__:
            script_dir = os.path.dirname(os.path.abspath(args.filenames[-1]))
            os.chdir(script_dir)
            self.send_response(self.iopub_socket, 'stream',
                  {'name': 'stdout', 'text': 'Current working directory is set to {}\n'.format(script_dir)})

    # add an additional magic that only useful for spyder
    def _do_execute(self, code, silent, store_history=True, user_expressions=None,
                   allow_stdin=False):
        code = self.remove_leading_comments(code)

        if self.MAGIC_EDIT.match(code):
            options, remaining_code = self.get_magic_and_code(code, False)
            self.handle_magic_edit(options)
            # self.options will be set to inflence the execution of remaing_code
            return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
        else:
            return super(SoS_SpyderKernel, self)._do_execute(code, silent, store_history, user_expressions, allow_stdin)

    def _reset_dict(self):
        super(SoS_SpyderKernel, self)._reset_dict()
        # spyder 3 executes commands such as
        #
        # get_ipython().kernel....
        #
        # so we will have to expose the shell to the SoS dictionary
        env.sos_dict.set('__ipython__', self.shell)
        SoS_exec('''
def get_ipython():
    return __ipython__
''', None)
        self.original_keys.add('__ipython__')
        self.original_keys.add('get_ipython')

if __name__ == '__main__':
    from ipykernel.kernelapp import IPKernelApp
    IPKernelApp.launch_instance(kernel_class=SoS_Kernel)
