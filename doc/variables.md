<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [SoS Variables](#sos-variables)
  - [Reserved and constant variables](#reserved-and-constant-variables)
  - [Command line options](#command-line-options)
  - [Runtime variables](#runtime-variables)
- [Python Expression with SoS variables](#python-expression-with-sos-variables)
- [Use of SoS variables and expressions (need reconsideration)](#use-of-sos-variables-and-expressions-need-reconsideration)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## SoS Variables

SoS variables are case-sensitive python variables associated with the execution of SoS scripts. For simplicity,
they can only be **string** or **list of strings**. They can be defined almost anywhere in a SoS script.

### Reserved and constant variables

* `home`: (string) home directory
* `workflow_name`: (string) name of workflow being executed


### Command line options

SoS looks for a `[parameters]` section for command line options and their default arguments. The format of each variable is 
```
# comment
var_name=default_value
```

The default value can be string or list of string, which has to be returned as result of an expression. For example

```
# path to tool gatk
gatk_path='GATK'
```

defines a variable `gatk_path` with default value `'GATK'`. 

```
# A list of sample names
sample_names=[]
```

defines a variable `sample_name` with default value `[]`. The comments before variable definition is meaningful
because they will appear in the help message of the script (`sos show script`).

The values of these parameter can be input from command line with `gatk_path` accept a single value, and `sample_name` accept a list of values. E.g.

```bash
sos run myscript.sos --gatk_path /path/to/gatk --sample_names A1 A2 A3
```

The default values will be used if the parameters are not redefined from command line.

### Runtime variables

For a SoS step in the following form

```[workflow_step]
#
# description of step
#
pre-input-variable

input:
	input_files, opt1=value1, opt2=value2

depends:
    dependent_files

pre-action-variable

StepAction()

post-action-variable
```

The following variables will be defined

* Before entering the step:
  * **step_index** (string): step of the current action 
  * **step_input** (list of strings): input from the last completed step, or `[]` if step option `no_input` is set.
* Before the processing of input files
  * **pre-input variables** (string or list of strings): variables defined before `input:`.
* After the processing of input files:
  * **`input`** (list of strings): selected input files. Depending on input options, `StepAction()` might be executed multiple times with different set of `input`. 
  * **file label variables** (list of strings): Labels of files in `input` if `labels` option is defined.
  * **loop variables** (string): value of loop variables if `for_each` or `nc_for_loop` option is defined 
  * **pre-action variables** (string or list of string): variables defined after `input` and before `StepAction`. They will be evalulated with each `input`.
* After the exeuction of `StepAction`
  * **step_output** (list of strings). All output files (if `StepAction` is executed mutliple times) are collected to a variable **step_output**. 
  * **post-action variables** (string or list of strings): variables defined after the execution of `StepAction`. **input** will not be available for use in this step. As a special case, **output** can be redefined to change the output of this step.

`step_output` will become the default input of the next step.


## Python Expression with SoS variables

**Any python expression involving any SoS variables and defined functions can be used to define variables**.
A sequence result will be converted to list of strings and all other result will be converted to string. For example

* `input[0]` gets the first input file
* `os.path.basename(input[0])` get the base name of input file
* `[x for x in input if '_R1_' in x]` returns all files with `_R1_` in its names
* `output_dir + '/' + input[0] + '.bai'` return a file under `output_dir` with `.bai` appended to filename.
* `os.path.getsize(input[0])` return a string (converted from integer) of file size
* `os.path.getsize(input[0]) > 0` return `True` (a string, converted from boolean True) if file is not empty.
* `os.environ['MY_VAR']` return the value of environment `MY_VAR`.
* `glob.glob('*.txt')` return a list of files with extension `.txt` under the current directory.

Note that SoS makes available modules `glob`, `os`, `sys`, and more modules can be imported with workflow action `import_module`.

Because Python expressions are very expressive, the possibility here is almost endless and provide SoS a great source of power
and flexibility.

## Use of SoS variables and expressions (need reconsideration)

A new variable can be created with 

```python
var=value
```

For example
```python
a='blah'
```

assigns value `'blah'` to variable `a`,

```python
a=input[0] + '.bai'
```

assigns value of `input[0]` plus `'.bai'` to vaiable `a`.

SoS variables and their expressions can also be used for text substitution in step actions. To differntiate 
expressions from other parts of the text, the expression should be quoted between `${` and `}`. Also, if the 
result of the expression is a list, it will be converted to a string by joining the strings with a single space.
Therefore,

* `${variable}` will be quivalent to `${variable[0]}` if it is a list of length  1. 
* `run('''process ${files}''')` will be translated to `run('''process file1.txt file2.txt''')` if variable 
  `files` has value `['file1.txt', 'file2.txt']`. This works most of the time but might fail if filename contains
  special characters. This can be avoided by using expressions such as `run('''process '${files[0]}' '${files[1]}'''')`.


If your script is expected to have patterns of `${variable}` that should not be interpreted by SoS, you can
define

```python
_sos_quotation=None
```

to disable variable intepretation or something like

```python
_sos_quotation=['[[', ']]']
```

to use an alternative quoting style ('[[expression]]' in this example). The setting is valid
for the current step only.
