<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [step options](#step-options)
- [Input files (`input:`)](#input-files-input)
- [Input options](#input-options)
  - [Passing input files all at once (default)](#passing-input-files-all-at-once-default)
  - [Passing files of allowed type (option `filetype`)](#passing-files-of-allowed-type-option-filetype)
  - [Passing files in groups (option `group_by`)](#passing-files-in-groups-option-group_by)
  - [Attaching variables to input filenames (option `labels`)](#attaching-variables-to-input-filenames-option-labels)
  - [Looping through values of a SoS variable (Option `for_each`)](#looping-through-values-of-a-sos-variable-option-for_each)
  - [Conditional skip of a step (option `skip`)](#conditional-skip-of-a-step-option-skip)
- [Dependent files (`depends:` (`dependent`? `requirements`?)](#dependent-files-depends-dependent-requirements)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

A SoS step generally has the following format

```
[name_step: option1, option2, ...]
#
# description of the step
#
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

### step options

SoS provides the following options

* `skip`: the whole step will be skipped as if it is defined
* `no_input`: specifies that this step does not have any input so the runtime signature does not depend on any input file. 
* `no_output`: specifies that this step does not have any output, so that the later steps can be executed before the end of this step. Because SoS will wait for the completion of a step if it does not have any defined output (through `output` optin of the workflow action), this option allows the following steps to be executed without waiting for the completion of the current step.  
* `input_alias`: a variable will be defined as the input files of the step
* `output_alias`: a variable will be defined as the output files of the step
* `blocking`: the step can only be executed by one instance of SoS. All other SoS instances will wait until one instance complete this step. This option should be used for actions such as the creation of index and downloading of resources.

Question: `no_output` is not very appropriate for its purpose because the action might specify `outpu`. Perhaps `non-blocking` should be used?

### Input files (`input:`)

The input of SoS step follows the following rules:

1. **the input of a SoS step is by default `${cmd_input}` for the first step and the output of the previous step otherwise**.
3. **step option `no_input` specifies that no input file is needed for the current step**
4. **Input of a step can be specified by item `input`**, which would be intepreted as filenames from command line. Wildcard characters (`*` and `?`) are always expanded. White space and other unusual characters in filenames need to be quoted.

Examples of input specification are as follows:

```
input:
	file1.fasta file2.fasta

input:
    file*.fasta filename\ with\ space.fasta
    
input:
    file1.txt
    directory/file2.txt
   
input:
    ${aligned_reads}

input:
    ${aligned_reads}
    ${reference_genome}

input:
    ${aligned_reads[2:]}

```

SoS variables such as `${aligned_reads}` and `${reference_genome}` will be expanded to strings and re-evaluated as input files. For example,

* `'file1.txt file2.txt'` will be translated to `file1.txt file2.txt` and recognized as two files 	`file1.txt` and `file2.txt`.
* `['file1.txt', 'file2 .txt']` will be expanded to `file1.txt file2\ .txt` and recognized as two files `file1.txt` and `file2 .txt`. (Note the space in filename).
* `'file*.txt'` will be transted to `file*.txt` and recognized as all files match pattern `file*.txt`.


### Input options 

The input options of a SoS step control how input files are passed to the step action.

#### Passing input files all at once (default)

This is the default action for a SoS step. That is to say, action in the following step will get all ``input_files`` as its input (variable `${input}`.

```
[step]
input:
    input_files

run('echo ${input}')

```

#### Passing files of allowed type (option `filetype`)

Emission options are appended to input file list as comma separated lists. A python parameter syntax is required.

Option `filetype` accepts one or more filetypes (file extension with `.`) or a lambda function. For example,

```
[step]
input:
	input_files
	: filetype='.fastq'
	
```

passes only files with extension `.fastq`.


```
[step]
input:
	input_files
	: filetype=['.fastq', '.fastq.gz']
	
```

passes only files with extension `.fastq` or `.fastq.gz`.



```
[step]
input:
	input_files
	: filetype=lambda x: open(x).readline().startswith('##fileformat=VCF4.1')
	
```

passes only files with the first line starting with string `##fileformat=VCF4.1``. Here the value of the parameter is a lambda function that will be passed to each input filename.


#### Passing files in groups (option `group_by`)

Option `group_by` pass input files in groups. For example,

```
[step]
input:
	file1 file2 file3 file4
	: group_by='single'

run('echo ${input}')

```

will execute `run('echo ${input}')` four times, with `${input}` set to `['file1']`, `['file2']`, `['file3']` and `['file4']` respectively.

Other values of `group_by` includes

* `'pairwise'`: yields two groups `['file1', 'file2']`, `['file2', 'file3']`, and `['file3', 'file4']`.
* `'combinations'`: yields six groups `['file1', 'file2']`, `['file1', 'file3']`,  `['file1', 'file4']`, `['file2', 'file3']`, `['file2', 'file4']` and `['file3', 'file4']`.
* ``pairs``: yields `['file1', 'file3']` and `['file2', 'file4']`. It basically split the input files in half and match files in the first half with files in the second half.

for four input files. Obviously, the output of the `pairs` cases depends on the order of files. If you need to pair files in any particular order, you can control it in input. For example

```
[step]
input:
	${sorted([x for x in cmd_input if '_R1_' in x])}
	${sorted([x for x in cmd_input if '_R2_' in x])}
	: group_by='pairs'

run('echo ${input}')

```

will take all input files and sort them by `_R1_` and `_R2_` and by filename. For example, four files `FEB_R1_1.txt FEB_R2_2.txt FEB_R1_2.txt FEB_R2_1.txt` will be sorted as `FEB_R1_1.txt FEB_R1_2.txt FEB_R2_1.txt FEB_R2_2.txt` and be sent to step action in two groups `['FEB_R1_1.txt', 'FEB_R2_1.txt']` and `['FEB_R1_2.txt', 'FEB_R2_2.txt']`.

#### Attaching variables to input filenames (option `labels`)

There are cases where the command line options or output directories depends on input filename. SoS allows you to label each filename with one or more variables so that they can be used accordingly. For example, if you have input files `bam_files` with values

```python
bam_files = ['case/A1.bam', 'case/A2.bam', 'ctrl/A1.bam', 'ctrl/A2.bam']
```

You might want to define some derived variables such as

```python
mutated = ${[x.split('/')[0] for x in bam_files]}
sample_name = ${[os.path.basename(x).split('.')[0] for x in bam_files]}
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
	${bam_files}
	: group_by='pairs', labels=['mutated', 'sample_name']

run('process ${input} with variables ${_mutated} and ${_sample_name}')

```

Here `bam_files` will be passed in pairs with file labels set to variables with `_` prefixed to their names. More specifically,

* Group1: `${input} = ['case/A1.bam', 'ctrl/A1.bam']`, `${_mutated}=['case', 'ctrl']`, `${_sample_name}=['A1', 'A1']`
* Group2: `${input} = ['case/A2.bam', 'ctrl/A2.bam']`, `${_mutated}=['case', 'ctrl']`, `${_sample_name}=['A2', 'A2']`

This is equivalent to

```
[step]
input:
	${bam_files}
	: group_by='pairs'

_mutated = ${[x.split('/')[0] for x in input]}
_sample_name = ${[os.path.basename(x).split('.')[0] for x in input]}

run('process ${input} with variables ${_mutated} and ${_sample_name}')

```

but it is cleaner because you do not have to do this each time when `$bam_files}` is used.

#### Looping through values of a SoS variable (Option `for_each`)

Option `for_each` allows you to repeat step actions for each value of a variable. For example, if

```python
method=['method1', 'method2']
```

You can repeat the analysis with each method using

```
[step]
input:
	${bam_files}
	: for_each='method'

run('Analyze ${input} with method ${_method}')

```

The step action will be executed twice with value of parameter `${_method}` set to `'method1'` and `'method2'` respectively.

Nested loops are also allowed. For example, if

```python
method=['method1', 'method2']
parameters=['-5', '-9']
```

You can execute 

```
[step]
input:
	${bam_files}
	: for_each=['method', 'parameter']

run('Analyze ${input} with method ${_method} and parameter ${_parameter}')

```

with parameters

* `_methods='method1', _parameter='-5'`
* `_methods='method1', _parameter='-9'`
* `_methods='method2', _parameter='-5'`
* `_methods='method2', _parameter='-9'`

Finally, if you would like to loop the action with several parameters, you can put them into a the same level using

```
[step]
input:
	${bam_files}
	: for_each='method,parameter'

run('Analyze ${input} with method ${_method} and parameter ${_parameter}')

```

The action will then be executed twice with parameters

* `_methods='method1', _parameter='-5'`
* `_methods='method2', _parameter='-9'`

#### Conditional skip of a step (option `skip`)

Option `skip=True` will make SoS skip the execution of the current step. Using `skip=True` is not very useful so this option is often used with a SoS variable. For example

```python
[10]
input:
	${fasta_files}
	: skip=${len(fasta_failes) == 1}
	
run('command to merge multiple fasta files.')
```

Here the `skip` option gets the value of `True` if there is only one input file. The command to merge multiple input files would then be skipped.

One important detail of this option is that the step is actually **executed with a null action** that passes input files to output files so this step still yields its `${output}`. In comparison, a step is completely ignored if it has step option `skip`. The consequence of this rule for this particular example is that its next step would get a merged file if there are multiple input files, or the original file if there is only a single input file.

### Dependent files (`depends:` (`dependent`? `requirements`?)

This item specifies files that are required for the step. Although not required, it is a good practice to list resource files and other dependency files for a particular step. For example

```
[10]
input:
	${fasta_files}
	
depends:
	${reference_seq}
	
Question: is there a need to test other requirements such as the value of environment variable, existence of other tools, content of file etc?


### Pre-action and post-action variables

Please refer to [use of SoS variables](variables.md) for details

### step actions

Please refer to [step actions](actions.md) for details



