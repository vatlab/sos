<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Command `sos`](#command-sos)
  - [subcommand `run` (see also `sos-runner`)](#subcommand-run-see-also-sos-runner)
  - [subcommand `show`](#subcommand-show)
  - [subcommand `export`](#subcommand-export)
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


### subcommand `run` (see also [`sos-runner`](#command-sos-runner))

```bash
sos run script [workflow[:steps]] [--input FILE1 FILE2 ...] [workflow-options] [-j NUM_JOBS]
```

Execute specified `steps` of `workflow` defined in `script` with specified `workflow-options`.

* `script` is a SoS script in [this format](sos_format_v1.md)
* `workflow` is one of the workflows defined in `script`. If left unspecified, the default (unnamed or with name `default`) workflow or the only workflow defined in `script` will be executed.
* `steps` should be prefixed with `:` and can be a comma separated list of `n` (number, specific step), `-n` (up to and include step `n`), `n-` (from step `n`). This option effectively add option `skip` to unspecified steps.
* `--input`: input files. Will be passed as variable `${cmd_input}` to the script.
* `workflow-options`: Options defined in the `[default]` section of `script`. variables with a text default value accepts a string input. variables with a list default value (e.g. `names=${[]}`)  accepts a list of strings.
* `-j`: Maximum number of concurrent jobs. A SoS script is by default executed sequentially (`-j 1`) but can have mutliple
  concurrent jobs if a positive number is specified. Please see [work flow control](workflow_control.md) in detail.

The list of acceptable workflow options for `script` can be displayed using command

```bash
sos show script
```
### subcommand `show`

```bash
sos show script [workflow:steps] 
```

Display details of `workflow` defined in `script`, including command line options defined by the `[default]` section of `script`.

### subcommand `export` 

```bash
sos export script [workflow:steps] [--input FILE1 FILE2 ...] [workflow-options] [-d OUTPUT_DIR] [-f]
```

Export `steps` of `workflow` defined in `script` to `script_dir`. This command will write `workflow_step.ext` for each `step` exported with appropriate file extension `ext` (e.g. `.R` for [R](https://www.r-project.org/) script) and `workflow.sh` to execute all the steps. Options of this command include

* `--input` and `workflow-options`: see command `sos run` for details.
*  `-d` (output directory): directory to which scripts are written. Default to current directory.
*  `-f` (force): overwrite existing files with different content silently. If unspecified, the command will fail with an error message in such cases. 

Please refere to [the export feature](export.md) for detailed examples of this command.

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

This allows a SoS script to be executed directory if it has a shebang line

```
#!/usr/bin/env sos-runner
```

and has appropriate permissions.

