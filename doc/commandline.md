<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Command `sos`](#command-sos)
  - [subcommand `run`](#subcommand-run)
  - [subcommand `show`](#subcommand-show)
  - [subcommand `export`](#subcommand-export)
  - [subcommand `analyze` (or `summarize`)](#subcommand-analyze-or-summarize)
  - [TBD Features](#tbd-features)
- [Command `sos-runner`](#command-sos-runner)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Command `sos`
Command `sos` accepts a number of subcommands (similar to `svn`, `git` etc). Its syntax follows

```bash
sos subcommand [subcommand-options]
```

You can use command

```bash
sos -h
```

to get a list of subcommands with brief descriptions and

```bash
sos subcommand -h
```

to get detailed description of a particular subcommand.


### subcommand `run`

```bash
sos run script [workflow[:steps]] [-f] [-h] [-j NUM_JOBS] [workflow-options]
```

Execute specified `steps` of `workflow` defined in `script` with specified `workflow-options`.

* `script` is a SoS script in [this format](sos_format_v1.md)
* `workflow` is one of the workflows defined in `script`. If left unspecified, the default (unnamed or with name `default`) workflow or the only workflow defined in `script` will be executed.
* `steps` should be prefixed with `:` and can be a comma separated list of `n` (number, specific step), `-n` (up to and include step `n`), `n-` (from step `n`). This option effectively add option `skip` to unspecified steps.
* `workflow-options`: Options defined in the `[default]` section of `script`. variables with a text default value accepts a string input. variables with a list default value (e.g. `names=${[]}`)  accepts a list of strings.
* `-j`: Maximum number of concurrent jobs. A SoS script is by default executed sequentially (`-j 1`) but can have mutliple
  concurrent jobs if a positive number is specified. Please see [work flow control](workflow_control.md) in detail.
* `-f` (fresh or force): ignore runtime signatures and re-execute the workflow even if the workflow has been executed before. This is useful for gathering runtime statistics by executing a workflow repeatedly.
* `-h`: Show usage of command `sos run` or `script` if a script is specified.

The list of acceptable workflow options for `script` can be displayed using command

```bash
sos show script
```
### subcommand `show`

```bash
sos show script [workflow:steps] [-d] [-h]
```
* `-d`: produce processing steps as a directional graph that can be pipelined to [the dot command of graphviz](http://www.graphviz.org/) for plotting.
* `-h`: display usage information

Display details of `workflow` defined in `script`, including command line options defined by the `[parameters]` section of `script`.

### subcommand `export` 

```bash
sos export script [workflow:steps] [workflow-options] [-d OUTPUT_DIR] [-f] [-h]
```

Export `steps` of `workflow` defined in `script` to `script_dir`. This command will write `workflow_step.ext` for each `step` exported with appropriate file extension `ext` (e.g. `.R` for [R](https://www.r-project.org/) script) and `workflow.sh` to execute all the steps. Options of this command include

*  `-d` (output directory): directory to which scripts are written. Default to current directory.
*  `-f` (force): overwrite existing files with different content silently. If unspecified, the command will fail with an error message in such cases. 
* `-h` (help): display help message.

### subcommand `analyze` (or `summarize`)

SoS saves runtime information during the execution of SoS script. Such information includes the start and end time of each step, CPU and memory usage, and the size of the current working directory. Users can learn the details of a previous run by running command

```bash
sos analyze
```

from the directory where the script was excuted. If the script has been executed multiple times, you can use option `-l` (list) to list all available runtime IDs,

```bash
sos analyze -l 
```

and use command

```bash
sos analyze ID
```

to display details of that particular run.

Other options include

* `-s`: summarize runtime information of all previous runs.
* `-o`: write a report to specified file

### TBD Features
* subcommand **`convert`**:
  It miight be useful to convert SoS scripts to other workflow language, at least partially. A complete translation is unlikely to be possible so SoS might simply export the SoS script to shell scripts (`sos export`) and create a pipeline for the master shell script. 

* subcommand **`admin`**:
  Reserved for miscellaneous adminstrative actions such as the setting of user options.
  
* subcommand **`edit`**: potentially a GUI viewer and editor for SoS scripts.
 
These features might be implemented in the future when needs arise.

## Command `sos-runner`

Command `sos-runner` is a short cut for ``sos run`` so

```bash
sos-runner script
```

is equivalent to

```bash
sos run script
```

This allows a SoS script to be executed directly if it is executable with shebang line

```
#!/usr/bin/env sos-runner
```

