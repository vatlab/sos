#!/usr/bin/env python3
#
# Copyright (C) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import time

from sos.controller import (request_answer_from_controller, send_message_to_controller)
from sos.eval import interpolate
from sos.targets import path, sos_targets
from sos.utils import env, pexpect_run

 #
 # docker support
 #


class SoS_DockerClient:
     '''A singleton class to ensure there is only one client'''
     _instance = None

     client = shutil.which('docker')
     pulled_images = set()

     def __new__(cls, *args, **kwargs):
         if not cls._instance:
             cls._instance = super(SoS_DockerClient, cls).__new__(cls)
         return cls._instance

     def total_memory(self, image='ubuntu'):
         '''Get the available ram fo the docker machine in Kb'''
         try:
             ret = subprocess.check_output(
                 f'''docker run -t {image} cat /proc/meminfo  | grep MemTotal''',
                 shell=True,
                 stdin=subprocess.DEVNULL)
             # ret: MemTotal:       30208916 kB
             self.tot_mem = int(ret.split()[1])
         except Exception:
             # some system does not have cat or grep
             self.tot_mem = None
         return self.tot_mem

def _is_image_avail(self, image):
         # the command will return ID of the image if it exists
         try:
             return bool(
                 subprocess.check_output(
                     f'''docker images {image} --no-trunc --format "{{{{.ID}}}}"''',
                     shell=True))
         except Exception as e:
             env.logger.warning(f'Failed to check image {image}: {e}')
             return False

 def _run_cmd(self, cmd, **kwargs):
         if env.config['run_mode'] == 'interactive':
             if 'stdout' in kwargs or 'stderr' in kwargs:
                 child = subprocess.Popen(
                     cmd,
                     shell=True,
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE,
                     bufsize=0)
                 out, err = child.communicate()
                 if 'stdout' in kwargs:
                     if kwargs['stdout'] is not False:
                         with open(kwargs['stdout'], 'ab') as so:
                             so.write(out)
                 else:
                     sys.stdout.write(out.decode())

                 if 'stderr' in kwargs:
                     if kwargs['stderr'] is not False:
                         with open(kwargs['stderr'], 'ab') as se:
                             se.write(err)
                 else:
                     sys.stderr.write(err.decode())
                 ret = child.returncode
             else:
                 # need to catch output and send to python output, which will in trun be hijacked by SoS notebook
                 ret = pexpect_run(cmd)
         elif '__std_out__' in env.sos_dict and '__std_err__' in env.sos_dict:
             if 'stdout' in kwargs or 'stderr' in kwargs:
                 if 'stdout' in kwargs:
                     if kwargs['stdout'] is False:
                         so = subprocess.DEVNULL
                     else:
                         so = open(kwargs['stdout'], 'ab')
                 elif env.verbosity > 0:
                     so = open(env.sos_dict['__std_out__'], 'ab')
                 else:
                     so = subprocess.DEVNULL

                 if 'stderr' in kwargs:
                     if kwargs['stderr'] is False:
                         se = subprocess.DEVNULL
                     else:
                         se = open(kwargs['stderr'], 'ab')
                 elif env.verbosity > 1:
                     se = open(env.sos_dict['__std_err__'], 'ab')
                 else:
                     se = subprocess.DEVNULL

                 p = subprocess.Popen(cmd, shell=True, stderr=se, stdout=so)
                 ret = p.wait()

                 if so != subprocess.DEVNULL:
                     so.close()
                 if se != subprocess.DEVNULL:
                     se.close()

             elif env.verbosity >= 1:
                 with open(env.sos_dict['__std_out__'],
                           'ab') as so, open(env.sos_dict['__std_err__'],
                                             'ab') as se:
                     p = subprocess.Popen(cmd, shell=True, stderr=se, stdout=so)
                     ret = p.wait()
             else:
                 p = subprocess.Popen(
                     cmd,
                     shell=True,
                     stderr=subprocess.DEVNULL,
                     stdout=subprocess.DEVNULL)
                 ret = p.wait()
         else:
             if 'stdout' in kwargs:
                 if kwargs['stdout'] is False:
                     so = subprocess.DEVNULL
                 else:
                     so = open(kwargs['stdout'], 'ab')
             elif env.verbosity > 0:
                 so = None
             else:
                 so = subprocess.DEVNULL

             if 'stderr' in kwargs:
                 if kwargs['stderr'] is False:
                     se = subprocess.DEVNULL
                 else:
                     se = open(kwargs['stderr'], 'ab')
             elif env.verbosity > 1:
                 se = None
             else:
                 se = subprocess.DEVNULL

             p = subprocess.Popen(cmd, shell=True, stderr=se, stdout=so)

             ret = p.wait()
             if so is not None and so != subprocess.DEVNULL:
                 so.close()
             if se is not None and se != subprocess.DEVNULL:
                 se.close()
         return ret