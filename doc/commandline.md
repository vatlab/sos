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


### subcommand `run` (see also `sos-runner`)

```bash
sos run script [workflow[:steps]] [--input FILE1 FILE2 ...] [workflow-options] 
```

Execute specified `steps` of `workflow` defined in `script` with specified `workflow-option`.

* `script` is a SoS script in [this format](doc/sos_format_v1.md)
* `workflow` is one of the workflow defined in `script`. If left unspecified, the default (unnamed or with name `default`) workflow or the only workflow defined in `script` will be executed.
* `steps` should follow be prefixed with `:` and can be a comma separated list of `n` (number, specific step), `-n` (up to and include step `n`), `n-` (from step `n`). This option effectively add option `skip` to unspecified steps.
* `--input`: input files. Will be passed as variable `${cmd_input}` to the script.
* `workflow-options`: Options defined in the `[default]` section of `script`. variables with a text default value accepts a string input. variables with a list default value (e.g. `names=${[]}`)  accepts a list of strings.

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

### Subcommands that might be implemented in the future
* subcommand **`convert`**:
  It miight be useful to convert SoS scripts to other workflow language, at least partially. A complete translation is unlikely to be possible so SoS might simply export the SoS script to shell scripts (`sos export`) and create a pipeline for the master shell script. 

* subcommand **`admin`**:
  Reserved for miscellaneous adminstrative actions such as the setting of user options.
  
* subcommand **`edit`**: potentially a GUI viewer and editor for SoS scripts.
 
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

