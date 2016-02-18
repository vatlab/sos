## SoS Variables

SoS variables are python variables associated with the execution of SoS scripts. They are **case sensitive** and can only be **string** or **list of strings**. 

### Reserved variables

* `home`: (string) home directory
* `workflow_name`: (string) name of workflow being executed


### Default command line option

* `cmd_input`: (list of strings) option to command line argument `--input`

### Command line options

SoS looks for a `[default]` section for command line options and their default arguments. The format of each variable is 
```
var_name=default_value
    description
```

The default value can be string or list of string, which has to be returned as result of an expression. For example

```
gatk_path=GATK
    path to tool gatk
```

defines a variable `gatk_path` with default value `'GATK'`. 

```
sample_names=${[]}
    A list of sample names
```

defines a variable `sample_name` with default value `[]`. 

The values of these parameter can be input from command line with `gatk_path` accept a single value, and `sample_name` accept a list of values. E.g.

```bash
sos run myscript.sos --gatk_path /path/to/gatk --sample_names A1 A2 A3
```

The default values will be used if the parameters are not redefined from command line.

### Runtime variables

For a SoS step in the following form

```[workflow_step: input_alias=ia, output_alias=oa]
input:
	input_files
	: emit options

depends:
    dependent_files

pre-action-variable

StepAction()

post-action-variable
```

The following variables will be defined


* Before entering the step:
  * **workflow_step** (string): step of the current action 
  * **input alias** (list of strings): if step option `input_alias` is set, tthis variable will include all input files of the step.
* After the processing of input files:
  * **`${input}`** (list of strings): selected input files. Depending on emit options, `StepAction()` might be executed multiple times with different set of `input`. 
  * **file labels** (list of strings): Labels of files in `input` if `labels` option is defined.
  * **loop variables** (string): value of loop variables if `for_each` option is defined 
  * **pre-action variables** (string or list of string): variables defined after `input` and before `StepAction`. They will be evalulated with each `input`.
* After the exeuction of `StepAction`
  * **${output}** (list of strings). All output files (if `StepAction` is executed mutliple times) are collected to a variable **${output}**. 
  * **post-action variables** (string or list of strings): variables defined after the execution of `StepAction`. **${input}** will not be available for use in this step. As a special case, **${output}** can be redefined to change the output of this step.
	* **output alias** (list of strings): if step option `output_alias` is set, this variable will be set to `${output}`.

`${output}` will become the default input of the next step.


## Python Expression with SoS variables

**Any python expression involving any SoS variables and defined functions can be used within `${}`**. A sequence result will be converted to list of strings and all other result will be converted to string. For example

* `${input[0]}` gets the first input file
* `${os.path.basename(input[0])}` get the base name of input file
* `${[x for x in input if '_R1_' in x]}` returns all files with `_R1_` in its names
* `${output_dir + '/' + input[0] + '.bai'}` return a file under `output_dir` with `.bai` appended to filename.
* `${os.path.getsize(input[0])}` return a string (converted from integer) of file size
* `${os.path.getsize(input[0]) > 0}` return `True` (a string, converted from boolean True) if file is not empty.
* `${os.environ['MY_VAR']}` return the value of environment `MY_VAR`.
* `${glob.glob('*.txt')}` return a list of files with extension `.txt` under the current directory.

Note that SoS makes available modules `glob`, `os`, `sys`, and more modules can be imported with workflow action `import_module`.

Because Python expressions are very expressive, the possibility here is almost endless and provide SoS a great source of power and flexibility.

## Creation of new variables

A new variable can be created with 

```python
var=value
```

For example
```python
a=blah
```

assigns value `'blah'` to variable `a`,

```python
a=${input[0] + '.bai'}
```

assigns value of `input[0]` plus `'.bai'` to vaiable `a`.

The only exception here is that `var` will be a list of string if `value` is an expression with list type. For example,

```python
a = ${input}
b = ${[os.path.basename(x) for x in input]}
c = ${[]}
```
will produce variables of list types.

## Use of SoS variables and expressions

SoS variables (and their expressions) are used slightly differently in different places, but the rules are largely intuitive.

1. In `input:` specification of SoS steps, the variables are treated as filenames or filenames with wildcard characters. An error will be generated if any of the resulting files do not exist.
2. As mentioned above, direct variable assignment (`var_name=${expression}`) will keep the type of expression.
3. In all other cases, `${expression}` will be converted to a string. List of strings will be joint by a space. 

Because of the third rule

* `${variable}` will be quivalent to `${variable[0]}` if it is a list of length  1. 
* `run('''process ${files}''')` will be translated to `run('''process file1.txt file2.txt''')` if variable `files` has value `['file1.txt', 'file2.txt']`. This works most of the time but might fail if filename contains special characters. This can be avoided by using expressions such as `run('''process '${files[0]}' '${files[1]}'''')`.


