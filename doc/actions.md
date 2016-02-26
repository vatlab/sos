<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [General interface](#general-interface)
- [Execution of scripts](#execution-of-scripts)
- [Utility actions](#utility-actions)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->
### General interface

```python
Action(
	script,
	stdout=None,
	stderr=None,
	interpreter=None,
	workdir=None,
	resource={'ncpu'=1, 'mem'=0},
  )
```

where
* `script`: the script to be executed
* `stdout` and `stderr`: standard output and standard error output of the script. These will be by default
   saved to the runtime directory.
* `interpreter`: command to execute the script. Additional command line options and commands can be used
   to execute the script in background (e.g. `bash {} &`), with additional option (e.g. `python -m module {}`),
   submit jobs (e.g. `cat {} | qsub -N `) or execute the script on another node (e.g. `ssh another bash {}`).
* `workdir`: directory where the command will be executed. Default to current working directory.
* `resource`: a parameter to specify the resource used by the action, which helps the scheduler how to 
   execute the step. Without finer control of the resource, each action is assumed to use 1 cpu. If action sets
   `ncpu=4` it will take 4 cpus from the specified `-j` value. 

Note that some options might be moved to the step level (e.g. resource).

### Execution of scripts

* ``run()``
* ``bash()``
* ``R()``
* ``perl()``
* ``python()``

### Utility actions

* ``fail_if()`` Terminate the execution of the pipeline if some condition is not met. This is SoS's version of  CWL's `requirements` feature.

* ``warn_if()`` Give a warning message if some condition is not met.

* ``show_vars()`` display all or specified pipeline variables

* ``check_cmd()`` check the output of command. Used to check the version and/or existence of command.


