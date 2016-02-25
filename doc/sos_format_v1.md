<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [shebang and file format lines](#shebang-and-file-format-lines)
- [Workflow descriptions](#workflow-descriptions)
- [Global variables](#global-variables)
- [Command line arguments](#command-line-arguments)
- [Defining workflows](#defining-workflows)
  - [Define single workflow](#define-single-workflow)
  - [Define multiple workflows](#define-multiple-workflows)
  - [Shared steps between workflows](#shared-steps-between-workflows)
  - [Execution of a subset of steps](#execution-of-a-subset-of-steps)
- [Workflow step](#workflow-step)
  - [step options](#step-options)
  - [Input files (`input:`)](#input-files-input)
  - [Input options](#input-options)
    - [Passing input files all at once (default)](#passing-input-files-all-at-once-default)
    - [Passing files of allowed type (option `filetype`)](#passing-files-of-allowed-type-option-filetype)
    - [Passing files in groups (option `group_by`)](#passing-files-in-groups-option-group_by)
    - [Attaching variables to input filenames (option `labels`)](#attaching-variables-to-input-filenames-option-labels)
    - [Looping through values of a SoS variable (Option `for_each`)](#looping-through-values-of-a-sos-variable-option-for_each)
    - [Conditional skip of a step (option `skip`)](#conditional-skip-of-a-step-option-skip)
  - [Dependent files (`depends` (or called `dependent`?))](#dependent-files-depends-or-called-dependent)
  - [Pre-input, pre-action and post-action variables](#pre-input-pre-action-and-post-action-variables)
  - [step actions](#step-actions)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

Although many items can be ignored, a typical SoS script has the following format:

```python
#!/usr/bin/env sos-runner
#fileformat=SOS1.0
#
# Other overall comments such as
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

[parameters]
# comment 1
par1=default1
   
# comment 2
par2=default2

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
but after the determination of variables from command line arguments. These are usually constant
variables such as location of various files.

### Command line arguments

Although sections can be arranged in any other, a secion named `[parameters]` is usually the first
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

### Defining workflows

A SoS script can specify one or more workflows. Each workflow consists of one or more numbered steps. The numbers specify the **logical order** by which the steps are executed, but a later step might be executed before the completion of previous steps if it does not depend on the output of these steps.

#### Define single workflow

A single pipeline can be specified without a name in a SoS script. For example, the following sections specify a pipeline with four steps `0`, `10`, `20`, and `100`. As you can see, the workflow steps can be specified in any order and do not have to be consecutive (which is actually preferred because it allows easy insertion of extra steps).

```
[0]
[20]
[10]
[100]

```

Workflows specified in this way is the `default` step and actually called `default` in SoS output. If you want to give it a meaningful name, you can specify the steps as 

```
[mapping_0]
[mapping_20]
[mapping_10]
[mapping_100]
```

Because the SoS script defines only one workflow, yo do not have to specify the name of workflow from SoS command 

```bash
sos run myscript.sos --input input1.fasta
```

#### Define multiple workflows 

A SoS script can define multiple workflows. For example, the following sections of SoS script defines two workflows named ``mouse`` and ``human``.

```
[mouse_10]
[mouse_20]
[mouse_30]
[human_10]
[human_20]
[human_30]
```

You will have to specify which workflow to execute from the command line, e.g.

```bash
sos run myscript mouse --input input1.fasta
```

If you would like to define a ``default`` and a named workflow, you can define them as

```
[10]
[20]
[30]
[test_10]
[test_20]
[test_30]
```

The `default` workflow will be executed by default using command

```bash
sos run myscript.sos --input input1.fasta
```

The `test` workflow will be executed if its name is specified from the command line

```bash
sos run myscript.sos test --input input1.fasta
```

#### Shared steps between workflows

The most common motivation of defining multiple workflows in a single SoS script is that they share certain processing steps. If this is the case, you can define sections such as

```
[mouse_10,human_10]
[mouse_20]
[human_20]
[mouse_30,human_30]
```

or

```
[*_10]
[mouse_20]
[human_20]
[*_30]
```

or

```
[*_10]
[mouse_20,human_20]
[fly_20]
[*_30,fly_50]
[fly_40]
```

In the last case, step defined by `[*_30,fly_40]` will be expanded to ``mouse_30``, ``human_30``, ``fly_30``, and ``fly_50`` and will be executed twice for the `fly` workflow.

#### Execution of a subset of steps

Although workflows are usually executed in its entirety, there are cases where you would like to execute only a subset of steps. For example, you can execute step 10 of the pipeline mouse using command

```bash
sos run myscript.sos mouse:0
```

Similarly, you can execute step 10 of the default workflow, up to step 20 of workflow `mouse``, steps 10 and 20 of workflow `mouse`, and step 20 and later of workflow `mouse`, respectively, using the following comands

```bash
sos run myscript.sos :10
sos run myscript.sos mouse:-20
sos run myscript.sos mouse:10,20
sos run myscript.sos mouse:20-
```

### Workflow step

Although only a *step action* is required for a SoS step, a complete SoS step can have the following form

```
[name_step: option1, option2, ...]
#
# description of the step
#
key0=value0

input:
    input files, opt1=value1, opt2=value2

depends:
    dependent files

key1=value1
key2=value2

step_action

key3=value3
key4=value4

```

#### step options

**Step options** are specified after step name and controls how the step will be executed. SoS provides the following options

* `skip`: the whole step will be skipped as if it is not defined at all in the script. This option provides a quick method to disable a step.
* `no_input`: this step does not need any input so it is a **root** of the execution tree. This option disconnects the current step with its previous steps so that it can be
    executed before the completion of previous steps. 
* `terminal`: this step is one of the **terminals** or **leafs** of the execution tree. This allows the later steps to be
   executed before the completion of this step. The step can have output but no other step should depend on these output
   files.
* `nonconcurrent`: if the step action will be repeated (using input options `group_by` or `for_each`), the loop actions are assumed to be parallel executable.
  If for some reason this assumption is wrong, you can set option `nonconcurrent` to let the actions execute sequentially.  
* `blocking`: the step can only be executed by one instance of SoS. All other SoS instances will wait until one instance complete this step. This option should be used for actions such as the creation of index and downloading of resources.
* `sigil`: alternative sigil of the step, which should be a string with space. E.g. `sigil='[ ]'` allows the use of expressions such as
  `[input]` in this step.

#### Input files (`input:`)

The input of SoS step follows the following rules:

1. **the input of a SoS step is by default the output of the previous step**.
3. **step option `no_input` specifies that no input file is needed for the current step**
4. **Input of a step can be specified by item `input`**, which should be a list or tupe of either 
    filenames (quoted strings) or SoS variables (use names of variables).  Wildcard characters
	(`*` and `?`) are always expanded. 

Examples of input specification are as follows:

```
input:
	'file1.fasta', 'file2.fasta'

input:
    'file*.fasta', 'filename with space.fasta'
    
input:
    'file*.txt',
    'directory/file2.txt'
   
input:
    aligned_reads

input:
    aligned_reads, reference_genome

input:
    aligned_reads[2:]

```

It does not matter if `aligned_reads` and `reference_genome` are lists of filenames. SoS will expand wildcard characters and
flatten the lists to a single list of filenames.


#### Input options 

The input options of a SoS step control how input files are passed to the step action. To should be keyword arguments appended to list
of input files.

##### Passing input files all at once (default)

This is the default action for a SoS step. That is to say, action in the following step will get all 
`input_files` as its input (variable `input`) to step action.

```
[step]
input:
    input_files

run('echo ${input}')

```

##### Passing files of allowed type (option `filetype`)

Emission options are appended to input file list as comma separated lists. A python parameter syntax is required.

Option `filetype` accepts one or more filetypes (file extension with `.`) or a lambda function. For example,

```
[step]
input:
	input_files, filetype='.fastq'
	
```

passes only files with extension `.fastq`.


```
[step]
input:
	input_files,
	filetype=['.fastq', '.fastq.gz']
	
```

passes only files with extension `.fastq` or `.fastq.gz`.



```
[step]
input:
	input_files,
	filetype=lambda x: open(x).readline().startswith('##fileformat=VCF4.1')
	
```

passes only files with the first line starting with string `##fileformat=VCF4.1``. Here the value of the parameter is a lambda function that will be passed to each input filename.


##### Passing files in groups (option `group_by`)

Option `group_by` pass input files in groups. For example,

```
[step]
input:
	'file1', 'file2', 'file3', 'file4',
	group_by='single'

run('echo ${input}')

```

will execute `run('echo ${input}')` four times, with `input` set to `['file1']`, `['file2']`, `['file3']` and `['file4']` respectively.

Other values of `group_by` includes

* `'pairwise'`: yields two groups `['file1', 'file2']`, `['file2', 'file3']`, and `['file3', 'file4']`.
* `'combinations'`: yields six groups `['file1', 'file2']`, `['file1', 'file3']`,  `['file1', 'file4']`, `['file2', 'file3']`, `['file2', 'file4']` and `['file3', 'file4']`.
* ``pairs``: yields `['file1', 'file3']` and `['file2', 'file4']`. It basically split the input files in half and match files in the first half with files in the second half.

for four input files. Obviously, the output of the `pairs` cases depends on the order of files. If you need to pair files in any particular order, you can control it in input. For example

```
[step]
input:
	sorted([x for x in fastq_files if '_R1_' in x]),
	sorted([x for x in fastq_files if '_R2_' in x]),
	group_by='pairs'

run('echo ${input}')

```

will take all input files and sort them by `_R1_` and `_R2_` and by filename. For example, four files
`FEB_R1_1.txt FEB_R2_2.txt FEB_R1_2.txt FEB_R2_1.txt` will be sorted as
`FEB_R1_1.txt FEB_R1_2.txt FEB_R2_1.txt FEB_R2_2.txt` and be sent to step
action in two groups `['FEB_R1_1.txt', 'FEB_R2_1.txt']` and `['FEB_R1_2.txt', 'FEB_R2_2.txt']`.

##### Attaching variables to input filenames (option `labels`)

There are cases where the command line options or output directories depends on input filename. SoS allows you to label each filename with one or more variables so that they can be used accordingly. For example, if you have input files `bam_files` with values

```python
bam_files = ['case/A1.bam', 'case/A2.bam', 'ctrl/A1.bam', 'ctrl/A2.bam']
```

You might want to define some derived variables such as

```python
mutated = [x.split('/')[0] for x in bam_files]
sample_name = [os.path.basename(x).split('.')[0] for x in bam_files]
```

with values

```python
mutated = ['case', 'case', 'ctrl', 'ctrl']
sample_name = ['A1', 'A2', 'A1', 'A2']
```

Then, if you are processing these files individually, or in pairs, you can attach values of these derived variables as follows:


```
[step]
input:
	bam_files,
	group_by='pairs', labels=['mutated', 'sample_name']

run('process ${input} with variables ${_mutated} and ${_sample_name}')

```

Here `bam_files` will be passed in pairs with file labels set to variables with `_` prefixed to their names. More specifically,

* Group1: `input = ['case/A1.bam', 'ctrl/A1.bam']`, `_mutated=['case', 'ctrl']`, `_sample_name=['A1', 'A1']`
* Group2: `input = ['case/A2.bam', 'ctrl/A2.bam']`, `_mutated=['case', 'ctrl']`, `_sample_name=['A2', 'A2']`

This is equivalent to

```
[step]
input:
	bam_files,
	group_by='pairs'

_mutated = [x.split('/')[0] for x in input]
_sample_name = [os.path.basename(x).split('.')[0] for x in input]

run('process ${input} with variables ${_mutated} and ${_sample_name}')

```

but it is cleaner because you do not have to do this each time when `bam_files` is used.

##### Looping through values of a SoS variable (Option `for_each`)

Option `for_each` allows you to repeat step actions for each value of a variable. For example, if

```python
method=['method1', 'method2']
```

You can repeat the analysis with each method using

```
[step]
input:
	bam_files,
	for_each='method'

run('Analyze ${input} with method ${_method}')

```

The step action will be executed twice with value of parameter `_method` set to `'method1'` and `'method2'` respectively.

Nested loops are also allowed. For example, if

```python
method=['method1', 'method2']
parameters=['-5', '-9']
```

You can execute 

```
[step]
input:
	bam_files,
	for_each=['method', 'parameter']

run('Analyze ${input} with method ${_method} and parameter ${_parameter}')

```

with parameters

* `_methods='method1', _parameter='-5'`
* `_methods='method1', _parameter='-9'`
* `_methods='method2', _parameter='-5'`
* `_methods='method2', _parameter='-9'`

If you would like to loop the action with several parameters, you can put them into a the same level using

```
[step]
input:
	bam_files,
	for_each='method,parameter'

run('Analyze ${input} with method ${_method} and parameter ${_parameter}')

```

The action will then be executed twice with parameters

* `_methods='method1', _parameter='-5'`
* `_methods='method2', _parameter='-9'`

Finally, option `for_each` assumes that the steps can be executed independently and concurrently and will
try to execute the actions in parallel if the script is executed in parallel mode (option `-j`). If for some
reason this assumption is wrong and the step needs to be executed sequentially, you should use the
`nonconcurrent` section option to execute loop actions sequentially.


##### Conditional skip of a step (option `skip`)

Option `skip=True` will make SoS skip the execution of the current step. Using `skip=True` is not very useful so this option is often used with a SoS variable. For example

```python
[10]
input:
	fasta_files,
	skip=len(fasta_failes) == 1
	
run('command to merge multiple fasta files.')
```

Here the `skip` option gets the value of `True` if there is only one input file. The command to merge multiple input files would then be skipped.

One important detail of this option is that the step is actually **executed with a null action** that passes input files to output files
so this step still yields its `step_output`. In comparison, a step is completely ignored if it has step option `skip`. The consequence
of this rule for this particular example is that its next step would get a merged file if there are multiple input files, or
the original file if there is only a single input file.

#### Dependent files (`depends` (or called `dependent`?))

This item specifies files that are required for the step. Although not required, it is a good practice to list resource files and other dependency files for a particular step. For example

```python
[10]
input:
	fasta_files
	
depends:
	reference_seq
```

#### Pre-input, pre-action and post-action variables

* Pre-input variables are defined before `input:` and evaluated before filenames are parsed.

* Pre-action variables are defined between `input:` and step action and evaluated (probably repeatedly) before the executation of 
  step action.

* Post-action variables are defined after step action and are evaluated after the completion of all step actions.

Please refer to [use of SoS variables](variables.md) for details

#### step actions

Please refer to [step actions](actions.md) for details




