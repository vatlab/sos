<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [format of SoS scripts](#format-of-sos-scripts)
  - [shebang and file format lines](#shebang-and-file-format-lines)
  - [Workflow descriptions](#workflow-descriptions)
  - [Global variables](#global-variables)
  - [Command line arguments](#command-line-arguments)
  - [SoS Steps](#sos-steps)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## format of SoS scripts

Although many items can be ignored, a typical SoS script has the following format:

```python
#!/usr/bin/env sos-runner
#fileformat=SOS1.0
#
# Other overal comments such as
#    License
#    change log
#

# Description of the workflows defined in
# script

# workflow1
# Description of workflow1

# workflow2
# Description of workflow2
#

# global variables
var1=value1
var2=value2
...

[default]
par1=default1
   comment1
   
par2=default2
   comment2

[step1]
step1 input, variables, and actions

[step2]
step2 input, variables, and actions
```

### shebang and file format lines

A SoS script usually starts with lines

```python
#!/usr/bin/env sos-runner
#fileformat=SOS1.0
```

The first line allows the script to be executed by command `sos-runner` if it is executed as
an executable script. The second line tells SoS the version of the script. The `#fileformat` 
line does not have to be the first line but should be in the first comment block.

### Workflow descriptions

The following comment blocks are description of workflows defined in the script. The description
is considered to be related to a particular workflow if it starts with only the name of the workflow
in the first line. These descriptions would be displayed in the output of command `sos view script`.

### Global variables

Variables defined before the first section will be evaluated before the execution of the workflow,
but after the determination of variables from command line arguments.

### Command line arguments

Although sections can be arranged in any other, a secion named `[default]` is usually the first
section of a SoS script. This section defined command line arguments of the script, their default
values and descriptions. The variables are defined in the format of


```
[parameters]
# comment 1
par1=default1
   
# comment2
par2=default2
```

where the default values determines the type of variable (string or list of strings). 
Please refer to [Command line options](variables#command-line-options) for more details.

### SoS Steps

A SoS script can define any number of sections for different steps of multiple workflows. 
The name of sections determines the workflow(s) they belong and the order at which 
they are executed. Please refer to [workflow sections](workflow_sections.md) for details.

Although most of the items are options, a complete SoS step follows the following format:


```
[name_step: option1, option2, ...]
#
# description of the step
#
key0=value0

input:
    input files
    : emit options

depends:
    dependent files

key1=value1
key2=value2

action_function

key3=value3
key4=value4

```

where

* **step options** controls how the steps are treated during execution
* **comments** of the step describes what this step does. These comments will
  be displayed in the output of command `sos show script`.
* **input** item specifies the input files of the step
* **input options** specifies how the input files are sent to the step action
* **dependent files** which are not part of the input but are prerequisite of the step. These files will be part of the step signature. The step will not be executed if the files are missing and will be re-executed if these files are changed.
* **pre-action variables** are defined before the action is executed.
* **post-action variables** are defined after the action is executed

Please refer to [step format](step_format.md) for details of these items.
