**Script of Scripts (SoS)** is a lightweight workflow system that helps you organize your commands and scripts in different languages into readable workflows that can be easily understood and modified by others. It is an easy-to-use alternative to specialized workflow systems such as [CWL](http://common-workflow-language.github.io/draft-3/) which makes it an ideal tool for the creation and maintenance of workflows that need to be frequently updated and shared with others.

Compared to maintaining multiple scripts, or using more specialized workflow systems such as [YAWL](http://www.yawlfoundation.org/), [CWL](http://common-workflow-language.github.io/), and [Galaxy](https://galaxyproject.org/), or make-file style workflow system such as [snakemake](https://bitbucket.org/johanneskoester/snakemake), 

* **SoS offers a way to organize your scripts in a single file**, which makes it easy to execute and maintain. You can include small and frequently changed commands and scripts in SoS and keep large and stable scripts in separate files.
* **SoS scripts are human readable and writable**. The workflow steps are logically arranged and the scripts and commands are kept mostly in their original form so it is easy for others to read a SoS workflow and modify it if needed. It follows [Python](http://www.python.org) syntax which is easy to understand even for users who do not know Python. In comparison, it would take a lot of time and practice to learn a specialized language to write a [CWL](http://common-workflow-language.github.io/) workflow (see [this CWL tutorial](https://github.com/common-workflow-language/workflows/wiki/Tutorial-DRAFT2) for an example).
* **SoS allows for notebook-, forward- and make-style definitions of workflows**. SoS does not limit you to particular way of defining a workflow. You can use SoS as a notebook to keep track of steps of analysis (notebook-style), transform the notes to a real workflow (forward-style), and add dependency rules if needed (make-style).
* **You can use SoS interactively with [iPython](https://ipython.org/) or [Jupyter](http://jupyter.org/), write and debug SoS scripts in [Spyder](https://pythonhosted.org/spyder/), and execute scripts in batch mode** with advanced workflow features. The workflow features of SoS are easy to use yet very powerful in helping you execute your pipelines efficiently not only locally, but also on cluster and cloud systems.

## SoS as a file format

Please refer to this **[tutorial](https://github.com/BoPeng/SOS/wiki/1.-Quick-Start)** for a quick start on SoS and the rest of the wiki for more details.

## Using SoS interactively

There are a number of ways for you to use SoS interactively. Unless you are familiar with iPython, Jupyter, or Spyder and already have a preferred method of working, the differences between these working environments can be confusing. The following table tries to list the characteristics of these methods to make it easier for you to get started.

Basically, you can use

1. `ipython + sos magic` if you are a diehard iPython user. However, this method is not practical for SoS user because of [a bug in ipython](https://github.com/ipython/ipython/issues/10072).
2. An editor + `qtconsole` if you strongly prefer a certain editor (vim)
3. `qtconsole` if you just want to play with sos for a bit
4. Jupyter if you like the notebook style or if you need to present code and results to others
5. Spyder if you like an IDE with separate editor and console windows. This should be a good choice for most users.

| | ipython kernel + sos magic | qtconsole | editor + qtconsole | jupyter notebook | Spyder|
|:---|:--- |:----- |:----- |:---- |:---- |
|**Description**| Native ipython environments with SoS magic | Console with sos kernel | External editor with help from a console window | Web interface with mixed scripts and results | IDE with separate editor and console | 
|**Good for** | ipython with separate sos environment| Test small pieces of code | Writing of serious SoS workflows with occasional need for testing | Interactive data analysis with mixed scripts and results (notebook style) | Interactive data analysis resulting in a complete script |
|**Command line**| `ipython`, `jupyter qtconsole`, `jupyter notebook` with Python kernel, or `spyder`| `jupyter qtconsole --kernel sos` | `jupyter qtconsole --kernel sos` | `jupyter notebook` with sos kernel| `spyder --kernel sos`| 
|**kernel**| ipython with sos magic | sos | sos | sos | sos |
|**Script editor**| Depends |  None | Editor of your choice (e.g. vim) | jupyter (web) | Spyder (editor) |
|**Enter commands via**| ipython command line | qtconsole (enter command)| qtconsole (copy/paste or `%paste` from editor) | jupyter (web) | Spyder (console)|
|**Subkernel support**| No, but can use other iPython magic such as `%%R`| Yes | Yes | Yes | Yes |
|**Preview support**| No | Yes | Yes | Yes | Yes |
|**Shell command**| Yes (through iPython) | Yes (through sos) | Yes | Yes | Yes |
|**remote access**| If using Jupyter | No | No | Yes through remote Jupyter server | Yes if connects to remote Jupyter server |
|**Magics**| Magics `%sos`, `%sosdict`, `%sospaste`, `%sosget`, `%sosput`, `%sosset` with ipython magics| `%with`, `%use`, `%paste`, `%set`, `%get`, `%preview`, `%run`, `%set`, `%restart`, `%dict`, `%cd`, ...| the same | the same | additional `%edit` magic|

## Writing and executing SoS workflow

