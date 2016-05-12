**Script of Scripts (SoS)** is a workflow engine that allows your to work
with scripts in different languages, glue them with Python expressions and
variables, and execute them as powerful workflows under different
environemtns. It is designed for data scienticists and bioinformatics who
routinely work with scripts in different languages (Python and R as being
most popular).

As an integration tool, SoS allows you to write single scripts with scripts
in different languages and execute them using a single command. It also
allows you to write and debug scripts in different languages in [Jupyter
notebook](http://jupyter.org/), with seamless integration of
multiple-kernels (e.g. python, and
[R](https://github.com/IRkernel/IRkernel) in a single notebook.

As a workflow engine, SoS helps you oraganize your commands and scripts in 
different languages into readable workflows that can be easily understood
and modified by others. It is an easy-to-use alternative to specialized workflow
systems such as [CWL](http://common-workflow-language.github.io/draft-3/) which makes it
an ideal tool for the creation and maintainance of workflows that need to be frequently updated and shared with others.


SoS is released under [GPL3](http://www.gnu.org/licenses/gpl-3.0.en.html). It supports Linux and Mac OSX systems and requires Python version 3.3 or higher. You can install the latest released version using command

```
% pip3 install sos
```

or compile the latest git version with commands

```
% git clone https://github.com/BoPeng/SOS.git
% cd SOS
% python3 setup.py install
```

Note that

* You might need to use command `pip` and `python` instead of `pip3` and `python3` if you have python 3 as the default python installation.
* If command `sos` is not found after installation, you will need to add paths such as
`/Library/Frameworks/Python.framework/Versions/3.4/bin/` to `$PATH` or
create symbolic links of `sos` and `sos-runner` commands in
`/usr/local/bin`.

Please find more information on **[SoS
wiki](https://github.com/BoPeng/SOS/wiki)**, or use the following links
directly:

* [Quick start guide](https://github.com/BoPeng/SOS/wiki/Quick-Start)
* [A presentation about SoS](https://github.com/BoPeng/SOS/wiki/SoS_March2016.pdf) (updated on Apr. 7th for version 0.5.7)
* [Complete documentation](https://github.com/BoPeng/SOS/wiki/Documentation)
