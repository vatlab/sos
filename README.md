**Script of Scripts (SoS)** is a workflow engine that allows you to work
with scripts in different languages, glue them with Python expressions and
variables, and execute them as powerful workflows under different environemtns. It is designed for data scienticists and bioinformatics who routinely work with scripts in different languages such as R, Python, Perl, and bash.

As a **notebook tool** that promotes [literate programming](https://en.wikipedia.org/wiki/Literate_programming), SoS allows you to perform and record your analysis in different languages in a single [Jupyter notebook](http://jupyter.org/), with seamless integration of multiple Jupyter kernels (e.g. python, and [R](https://github.com/IRkernel/IRkernel). Main features of SoS notebook are demonstrated in a [tutorial](https://github.com/BoPeng/SOS/blob/master/examples/NotebookTutorial.ipynb) and an [example](https://github.com/BoPeng/SOS/blob/master/examples/example.ipynb).

As a **workflow engine**, SoS helps you oraganize your commands and scripts in different languages into readable workflows that can be easily understood and modified by others. The workflows can be specified in both forward (step by step), makefile (dependency rules) and even a mixture of both styles. Beacuse of a great emphasis on readability, SoS is an easy-to-use alternative to specialized workflow systems such as [CWL](http://common-workflow-language.github.io/draft-3/) which makes it an ideal tool for the creation and maintainance of workflows that need to be frequently updated and shared with others.

## Installation

SoS is released under [GPL3](http://www.gnu.org/licenses/gpl-3.0.en.html). It supports Linux and Mac OSX systems and requires Python version 3.3 or higher. You can install the latest released version using command

```
% pip3 install sos
```

or compile the latest git version with commands

```
% git clone https://github.com/BoPeng/SOS.git
% cd SOS
% pip3 install sos
% python3 setup.py install
```

Note that

* You might need to use command `pip` and `python` instead of `pip3` and `python3` if you have python 3 as the default python installation.
* If command `sos` is not found after installation, you will need to add paths such as
`/Library/Frameworks/Python.framework/Versions/3.4/bin/` to `$PATH` or
create symbolic links of `sos` and `sos-runner` commands in
`/usr/local/bin`.

If you are using docker, you can run SoS directly using command

```
% docker run -it mdabioinfo/sos:latest /bin/bash
```
or using another container 

```
% docker run -it mdabioinfo/sos-full:latest /bin/bash
```
with [R](https://www.r-project.org/) and [IRkernel](https://github.com/IRkernel/IRkernel), and a running [Jupyter](http://jupyter.org/) server with SoS, python, and R kernels.

## Documentation

Please find more information on **[SoS
wiki](https://github.com/BoPeng/SOS/wiki)**, or use the following links
directly:

* [Quick start guide](https://github.com/BoPeng/SOS/wiki/Quick-Start)
* [Tutorial](https://github.com/BoPeng/SOS/blob/master/examples/NotebookTutorial.ipynb) and [example](https://github.com/BoPeng/SOS/blob/master/examples/example.ipynb) of Jupyter notebooks.
* [A presentation about SoS](https://github.com/BoPeng/SOS/wiki/SoS_March2016.pdf) (updated on Apr. 7th, 2016 for version 0.5.7, already quite outdated given the rapid evolution of SoS)
* [Complete documentation](https://github.com/BoPeng/SOS/wiki/Documentation)
