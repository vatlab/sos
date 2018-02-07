[![PyPI version](https://badge.fury.io/py/sos.svg)](https://badge.fury.io/py/sos)
[![Build Status](https://travis-ci.org/vatlab/SoS.svg?branch=master)](https://travis-ci.org/vatlab/SoS)
[![Build Status](https://ci.appveyor.com/api/projects/status/x092eusa0tta3msw?svg=true
)](https://ci.appveyor.com/project/BoPeng/sos)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/67b766a827fb491fa473032b4f70ebb7)](https://www.codacy.com/app/BoPeng/SoS?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=vatlab/SOS&amp;utm_campaign=Badge_Grade)
[![Coverage Status](https://coveralls.io/repos/github/vatlab/SOS/badge.svg)](https://coveralls.io/github/vatlab/SOS)

**Script of Scripts (SoS)** is a **file format** for the inclusion of
scripts in multiple languages into a single executable script, an
**interactive environment** and **notebook tool** for working with different scripts, and
a **workflow engine** for the execution of workflows consisting with scripts
in different languages. It is designed for data scienticists and bioinformatics who routinely work with scripts in different languages such as R, Python, Perl, and bash.

Please refer to the [SoS homepage](http://vatlab.github.io/SOS) for more information.

### Change Log of SoS and SoS Notebook

SoS 0.9.11.3
* [sos#879](https://github.com/vatlab/SoS/issues/879): Add action options `stdout` and `stderr` to redict output from script-executing actions. 
* [sos-notebook#42](https://github.com/vatlab/sos-notebook/issues/42): Add option `--append` to magic `%capture` .

SoS 0.9.11.2
* [sos-notebook#39](https://github.com/vatlab/sos-notebook/issues/39): Separation installation and deployment and use command `python -m sos_notebook.install` to install `sos` kernel to Jupyter.

SoS 0.9.10.19

* [sos#874](https://github.com/vatlab/SoS/issues/874): Add input option `concurrent=True` to allow parallel execution of input groups.
* [sos#874](https://github.com/vatlab/SoS/issues/874): Optimize task submission of task engines to reduce status checking 

SoS Notebook 0.9.10.17

* [sos-notebook#32](https://github.com/vatlab/sos-notebook/issues/32): Add magic `%capture` to capture output of cells 
