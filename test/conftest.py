import pytest
import uuid
import os
import shutil
import tempfile
import textwrap
import pathlib
import yaml


@pytest.fixture
def config_factory():
    filename = tempfile.NamedTemporaryFile(suffix='.yml', delete=False).name

    def get_config(text_or_dict):
        with open(filename, 'w') as conf:
            if isinstance(text_or_dict, str):
                conf.write(textwrap.dedent(text_or_dict))
            elif isinstance(text_or_dict, dict):
                yaml.dump(text_or_dict, conf)
            else:
                raise ValueError(
                    'A text or dictionary is expected for config_factory.')
        return filename

    yield get_config
    os.remove(filename)


@pytest.fixture
def script_factory():
    filename = tempfile.NamedTemporaryFile(suffix='.sos', delete=False).name

    def get_script(text):
        with open(filename, 'w') as conf:
            conf.write(textwrap.dedent(text))
        return filename

    yield get_script
    os.remove(filename)


@pytest.fixture
def temp_factory():

    temp_fds = []
    temp_dirs = []

    def get_tempfiles(*args, **kwargs):
        content = kwargs.get('content', None)
        for names in args:
            if isinstance(names, str):
                names = [names]
            for name in names:
                if content is None:
                    pathlib.Path(name).touch()
                else:
                    with open(name, 'w') as tf:
                        tf.write(content)
                temp_fds.append(name)
        if 'dir' in kwargs:
            if isinstance(kwargs['dir'], str):
                dirs = [kwargs['dir']]
            else:
                dirs = kwargs['dir']
            for dir in dirs:
                if os.path.isdir(dir):
                    shutil.rmtree(dir)
                os.makedirs(dir, exist_ok=True)
                temp_dirs.append(dir)
        return temp_fds, temp_dirs

    yield get_tempfiles

    for temp_file in temp_fds:
        try:
            if os.path.isfile(temp_file):
                os.remove(temp_file)
        except Exception:
            pass

    for temp_dir in temp_dirs:
        try:
            if os.path.isdir(temp_dir):
                shutil.rmtree(temp_dir)
        finally:
            pass


@pytest.fixture
def clear_now_and_after():
    # This fixture takes one or more file names or directory names and will
    # remove them if they exist, and also after the tests are completed.
    temp_fds = []

    def clear_files_and_dirs():
        for temp_fd in temp_fds:
            try:
                if os.path.isfile(temp_fd):
                    os.remove(temp_fd)
                elif os.path.isdir(temp_fd):
                    shutil.rmtree(temp_fd)
            except Exception:
                pass

    def get_names(*args):
        for names in args:
            if isinstance(names, str):
                temp_fds.append(names)
            else:
                temp_fds.extend(names)

        clear_files_and_dirs()
        return temp_fds

    yield get_names

    clear_files_and_dirs()

@pytest.fixture
def sample_workflow():
    with open('sample_workflow.ipynb', 'w') as sn:
        sn.write(r'''{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "metadata": {
    "kernel": "SoS"
   },
   "outputs": [],
   "source": [
    "# this is a test workflow"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {
    "kernel": "SoS"
   },
   "source": [
    "This is a markdown cell"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "metadata": {
    "kernel": "R"
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "This is a cell with another kernel"
     ]
    }
   ],
   "source": [
    "cat('This is a cell with another kernel')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "metadata": {
    "kernel": "SoS"
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "This is a scratch cell\n"
     ]
    }
   ],
   "source": [
    "print('This is a scratch cell')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {
    "kernel": "SoS"
   },
   "outputs": [],
   "source": [
    "# this comment will be included but not shown in help message\n",
    "# because it is for the global\n",
    "[global]\n",
    "a = 1\n",
    "# this comment will become the comment for parameter b\n",
    "parameter: b=2\n",
    "parameter: c=3\n",
    "# this comment will become the comment for parameter d\n",
    "parameter: d='d'"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {
    "kernel": "SoS"
   },
   "outputs": [],
   "source": [
    "# this comment will not be included in exported workflow\n",
    "# because it is not immediately before section\n",
    "\n",
    "# this is a section comment, will be displayed\n",
    "[default]\n",
    "print(f'Hello {a}')"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "SoS",
   "language": "sos",
   "name": "sos"
  },
  "language_info": {
   "codemirror_mode": "sos",
   "file_extension": ".sos",
   "mimetype": "text/x-sos",
   "name": "sos",
   "nbconvert_exporter": "sos_notebook.converter.SoS_Exporter",
   "pygments_lexer": "sos"
  },
  "sos": {
   "default_kernel": "SoS",
   "kernels": [
    [
     "R",
     "ir",
     "R",
     "#DCDCDA"
    ],
    [
     "SoS",
     "sos",
     "",
     ""
    ]
   ],
   "panel": {
    "displayed": true,
    "height": 0,
    "style": "side"
   },
   "version": "0.9.14.11"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}
''')