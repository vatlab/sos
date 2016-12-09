**Script of Scripts (SoS)** is a **file format** for the inclusion of
scripts in multiple languages into a single executable script, an
**interactive environment** and **notebook tool** for working with different scripts, and
a **workflow engine** for the execution of workflows consisting with scripts
in different languages. It is designed for data scienticists and bioinformatics who routinely work with scripts in different languages such as R, Python, Perl, and bash.

As a **file format** that help you record all your analysis into a single
script with minimal effort, SoS allows you to include scripts in different
languages literally into a `.sos` files, and execute it as a single
workflow. The basis of this format is demonstrated in our [Quick start
guide](https://github.com/BoPeng/SOS/wiki/1.-Quick-Start).

As an **interactive environment** and **notebook tool** that promotes [literate
programming](https://en.wikipedia.org/wiki/Literate_programming), SoS
allows you to perform and record your analysis in different languages in a
single [Jupyter notebook](http://jupyter.org/), with seamless integration
of multiple Jupyter kernels (e.g. python, and
[R](https://github.com/IRkernel/IRkernel)). You can check out 
example notebooks and create new ones from [our public Jupyter
server](http://ec2-34-192-184-206.compute-1.amazonaws.com:8000/).

As a **workflow engine**, SoS helps you oraganize your commands and scripts in different languages into readable workflows that can be easily understood and modified by others. The workflows can be specified in both forward (step by step), makefile (dependency rules) and even a mixture of both styles. Beacuse of a great emphasis on readability, SoS is an easy-to-use alternative to specialized workflow systems such as [CWL](http://www.commonwl.org/) which makes it an ideal tool for the creation and maintainance of workflows that need to be frequently updated and shared with others.

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

The following packages are needed if you are using SoS with R

* [R](https://www.r-project.org/) version 3.2 or later (for IRKernel)
* [IRKernel](https://github.com/IRkernel/IRkernel) R kernel for Jupyter
* Python [feather-format](https://github.com/wesm/feather) module for
exchanging data frames between SoS/Python DataFrame and R data.frame. This
module is not installed by default because it need to compile from source
or [install from conda](https://pypi.python.org/pypi/feather-format).

## Running sos in Docker 

If you are using docker, you can run SoS directly using command

```
% docker run -it mdabioinfo/sos:latest /bin/bash
```

to enter a command prompt with sos command. More usefully, you can start a
Jupyter server with [R](https://www.r-project.org/) and [IRkernel](https://github.com/IRkernel/IRkernel),
Julia, Python, and SoS kernels, and many Python and R modules for data sciencists using command


```
% docker run -d -p 8888:8888 mdabioinfo/sos-notebook
```

After the docker is running in the background, you can start a browser and
start working with a complete SoS environment with URL

```
http://localhost:8888
```

You can even use this docker image for your daily data analysis if you make
your local directory available to the Jupyter server using command 

```
% docker run -d -p 8888:8888 -v $HOME:/home/jovyan/work  mdabioinfo/sos-notebook
```

This command mounts your home directory (`$HOME`) to the docker machine but
you can specify any local directory. This container is hosted at [our public Jupyter
server](http://ec2-34-192-184-206.compute-1.amazonaws.com:8000/) from which you can open our sample
notebooks and create your own notebooks without installing anything.

## Documentation

Please find more information on **[SoS
wiki](https://github.com/BoPeng/SOS/wiki)**, or use the following links
directly:

* [Quick start guide](https://github.com/BoPeng/SOS/wiki/1.-Quick-Start)
* [Tutorial](https://github.com/BoPeng/SOS/blob/master/examples/NotebookTutorial.ipynb) and [example](https://github.com/BoPeng/SOS/blob/master/examples/example.ipynb) of Jupyter notebooks.
* [A presentation about SoS](https://github.com/BoPeng/SOS/wiki/SoS_March2016.pdf) (updated on Apr. 7th, 2016 for version 0.5.7, already quite outdated given the rapid evolution of SoS)
* [Complete documentation](https://github.com/BoPeng/SOS/wiki)
* Tutorials on using SOS with
[iPython](https://github.com/BoPeng/SOS/wiki/3.-Using-SoS-with-iPython),
[Jupyter](https://github.com/BoPeng/SOS/wiki/4.-SoS-Notebook-Using-Jupyter),
and [Spyder](https://github.com/BoPeng/SOS/wiki/5.-Using-Spyder-as-SoS-IDE).
* A complex [RNASeq data analysis workflow](https://github.com/BoPeng/SOS/wiki/6.-A-Complete-Example) written in SoS.
