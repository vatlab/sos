<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Data model](#data-model)
  - [Unifying SoS dictionary and list?](#unifying-sos-dictionary-and-list)
  - [Read only SoS variable?](#read-only-sos-variable)
- [File format](#file-format)
  - [Default input of step?](#default-input-of-step)
  - [design of `for_each`](#design-of-for_each)
  - [Enforce naming convention?](#enforce-naming-convention)
- [Workflow features](#workflow-features)
  - [Runtime control](#runtime-control)
  - [Resource control](#resource-control)
  - [Nested workflow](#nested-workflow)
  - [Portability of runtime signatures](#portability-of-runtime-signatures)
  - [Libraries](#libraries)
- [Actions](#actions)
  - [Session info?](#session-info)
- [External support](#external-support)
  - [Support for docker](#support-for-docker)
- [Resolved](#resolved)
  - [Other requirements of steps](#other-requirements-of-steps)
  - [Handling of filenames with spaces and other special characters](#handling-of-filenames-with-spaces-and-other-special-characters)
  - [Section option `concurrent`](#section-option-concurrent)
  - [Default parameter `--input`](#default-parameter---input)
  - [variable definition](#variable-definition)
  - [A more pythonic approach?](#a-more-pythonic-approach)
  - [Potential conflict of `${}` with other languages](#potential-conflict-of--with-other-languages)
  - [Format of input options](#format-of-input-options)
  - [Use of '${}' in places other than script string substitution](#use-of--in-places-other-than-script-string-substitution)
  - [trouble with sigil usage](#trouble-with-sigil-usage)
  - [Moving parameters of action to section level](#moving-parameters-of-action-to-section-level)
  - [Alternative or configurable global sigil?](#alternative-or-configurable-global-sigil)
  - [Section option `no_input`](#section-option-no_input)
  - [Section option `terminal` and `starting`](#section-option-terminal-and-starting)
  - [Use of dictionary variable](#use-of-dictionary-variable)
  - [Use `None` instead of `[]` for no input or output](#use-none-instead-of--for-no-input-or-output)
  - [Dynamic output](#dynamic-output)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Data model

### Unifying SoS dictionary and list?

Snakemake has some magic that allows `input` to be accessed as both list and object with attributes. For example, 
it allows

```
input:  reference='seq.fasta'
        bamfiles=['a.bam', 'b.bam']
```

and expressions such as `input.reference`, `input.bamfile`, and `input`. This can be easily achieved with some Python magic
but I am not sure if I should follow it and make SoS more difficult to understand.

```Python
class SoS_List(Collections.OrderedDict):
    def __getitem__(self, key):
	    # access as dictionary if key is a name
		#
		# input.rawfile
		#
		# access as a list if key is a number
		#
		# input[key]
		#
	
	def __repr__(self):
	    # connect all values() as a list ...
		#
		# ${input}
		#
	
```

In the end, we can claim that SoS list is 

* A list of strings
* Single or slice of the list can have a name
* Items can be access by integer index, attribute, or dictionary key

A SoS list can be initialized with

* A list: `['seq.fasta', 'a.bam', 'b.bam']`
* A dictionary: `{'reference': 'seq.fasta`, 'bamfiles': ['a.bam', 'b.bam']}`
* Parameters: `reference='seq.fasta', bamfiles=['a.bam', 'b.bam']`

Elements could be accessed as

* `a[1]`, `a[2:]`
* `a.reference`, `a.bamfiles` (if names are given)
* `a['reference']`, `a['bamfiles']`

Cons:
* `a.reference` should be a string or list of strings (for consistency) 
* It is difficult to change this datatype (we can make such variables readonly)
* It does not work well with input options such as `group_by`.
* Unlike snakemake, SoS already allows definition of variables, so users can do

```python
reference='seq.fasta'
bamfiles=['a.bam', 'b.bam']

input: reference, bamfiles

```
and use `input`, `reference`, and `bamfiles` separately. Here `reference` more likely 
belongs to `depends`. 

For simplicity and mostly the last reason, I incline not to use complicated data structure.


### Read only SoS variable?

A SoS can be easier to understand if we make most SoS variables readonly. That is to say, a SoS variable can not be changed 
after it has been initialized. Exceptions to this rule can be system variables and temporary looping variables. 


A less stringent version of this rule is that variables can only be replaced in their entirety, which disallows
changing part of the file list.


## File format
### Default input of step?

Would it clearer to require explicit input files? Right now the a step's `step_input` is 
the output of previous step. The problem with no default input is that output alias are 
almost always required and make the script a bit cumbersome.

Pros and cons of using default input:

pro:

* Avoid specifying input repeated
* Avoid using output alias. This could be avoided by adding runtime variables such as `input100`, but `input: output200` is not
  a good idea because the current step will be broken if step indexes are changed. A properly named variable is much better 
  for this purpose.

cons:

* Difficult to differentiate between 'No input' and 'None input' because SoS will have to look at its previous step
  to determine the current status.



### design of `for_each`

I am not sure the current design is intuitive. It requires the variables to be defined before and uses a derived variable (`_name` for variable `name`), but the usual `for_each` loop does not fit into the existing parameter structure.

Another concern is that we are using `input` for both regular and looped cases. It might
make sense to use `_input` for the latter.


### Enforce naming convention?

Does it make sense to force some name conversion so that users immediately know the type and nature of variables? This might
make the script a bit more readable. For example

1. Constant being string and all captical letters. (e.g. `RESOURCE_DIR='/path/to/resource'`) 
2. Derived or temporary variables having leading underscore. (e.g. `_loop_variable`, `_label`, `_input`, right now `input` is used)
3. String and list of string having different name conversion? I hate perl's `$` and `@` symbols though.

By 'enforce', I mean SoS can give warning even error if a variable's usage does not match its name convention.

## Workflow features

### Runtime control

Limiting the files that the directory the step action can write to? 

### Resource control

Limiting or monitoring the RAM and CPU (cores) the step uses?


### Nested workflow

Not sure how this would work, but in SoS one can certainly define a step with

```python
run('sos-runner anotherworkflow ${input}'
```

to execute another workflow.


### Portability of runtime signatures

Runtime signatures are usually not portable because filenames and pathes would change when the project
is moved to another location. However, because it is possible for a project to be archived, and restored
for some further analysis, saving the runtime signatures with the project might be usable. In addition,
if runtime signatures include standard and error output of the commands, it can be useful to keep them
for later reference. The interface can be

```bash
sos admin --pack_runtime runtime.db
```

for packing runtime information to a file,


```bash
sos admin --unpack_runtime runtime.db
```

for unpacking runtime information of the current project, and

```bash
sos admin --unpack_runtime runtime.db
```

The `sos view` or `sos edit` commands might be used to show standard and error output of commands.

### Libraries

Libraries would be python modules with defined SoS actions, but how to maintain and import these modules require further investigation. Furthermore, extensive use of libraries somehow beats the purpose of SoS (readability) because libraries hide the details of actions.


## Actions

### Session info?

There is a need to save the version of program or packages of R ...

## External support

### Support for docker 

What sure what is involved/required to support docker.

## Resolved

### Other requirements of steps

CWL has the requirement keyword which handles all sorts of requirement. SoS solves this by

* `depends:` item in pipeline step, which lists all dependent files that will go into step signature.
* `fail_if()` action, which stops the workflow if the requirements are not satisfied.

### Handling of filenames with spaces and other special characters

Because there is no way SoS knows how individual scripts and commands handle filenames with spaces, SoS would need to

1. Read filenames correctly from command line. Users might use `'file name with space.txt'` (with single quote) or `file\ name\ with\ space.txt` (with backslash, without quote) but SoS will get the correct name as `'file name with space.txt'`.
2. Filenames will be used internally so that it can be treated correct by SoS variables.
3. It is users to responsibility to use these variables correctly in the script. For example

```python
run('cat ${input}')
```

will fail with `input=['file name with space.txt']` because the command will be translated to 

```
cat file name with space.txt
```

Users should use 

```python
run(''' cat '${input}' ''')
```

or

```python
run(''' cat ${shlex.quote(input[0])} ''')
```

if they expect filenames with special characters.

On the other hand,

```python
R('''
open("${input}")

''')
```

should work correctly and it would be wrong if SoS mangles `${input}` during variable substitution.


### Section option `concurrent`


There are some options to allow concurrent execution of step actions.

1. Default to concurrent but allow section option `nonconcurrent`, because the actions should be safe 
  to execute in parallel most of the time (processing input files one by one or in pair, or
  with different options).

2. Do not use section option, but specify this in input parameters. E.g. `for_each` as
  concurrent for each, and `nc_for_each` for nonconcurrent for each. Perhaps 
  `for_all` as nonconcurrent? The problem is that `group_by` would also generate
  multiple actions so we might also need `nc_group_by`.

Decided to use section option `nonconcurrent`)

### Default parameter `--input`

We do not have to allow a default parameter `--input`. It is easier to use 
but the parameter itself is not documented like other command line parameters. It might  
be better to force the definition of all parameters in the `[default]` section.

Decided to remove `input`.

### variable definition

The original design uses

```
path=/default/path
sample_names=${['a', 'b']}
```

without quotation marks. Which can be confusing because we sometimes require qutation marks
(e.g. `group_by='single' in input options)` and sometimes do not. 

The current design uses all python syntax

```
path='/default/path'
sample_names=['a', 'b']
```

which is more consistent because the right hand side are always valid python expressions.

### A more pythonic approach?

Currently there are pices of the script that is not python, most notably the input
specification. It might make sense to turn them all to python syntax. Even wrap everything
to huge function calls, such as

```
step(index=10,
	input=...,
	requires=...,
	action=...,
)
```

This too ugly to adopt.

### Potential conflict of `${}` with other languages 

`${}` is used in shell, `` ` ` `` is used in python (Python 2 only) and perl, `{}` is commonly used in many languages, and the list goes on. 

SoS uses regular expression to identify expressions and text sustitution to compose script so it should be flexible
in which quotation style to use, we can use a section option `sigil` to define qutations for a particular section.

e.g.

```
[10: sigil='[[ ]]' ]

run('command with [[input]]')

```

Here the value to `quotation` should be a string with a single separting two pieces of quote.



### Format of input options

Currently input options are appended to input specification. Perhaps we can ignore the second `:` and use `,` all over. I mean, we can change

```python
input:  fasta_files: group_by='single'
```

to 

```python
input:  fasta_files, group_by='single'
```

because filenames usually do not have this form. This also make the parsing of `input:` easier because we input options can be parsed with input. In addition

```python
input: group_by='single'
```

is easier to understand than

```python
input: : group_by='single'
```

Decided to use `,` instead of `:`

### Use of '${}' in places other than script string substitution 

For example, replace `${resource_path}` with the string in `input` specification:

```
input:
	'${resource_path}/hg19/wholeGenome.fasta'
```

Using proper python expression would require

```
input:
	resource_path + '/hg19/wholeGenome.fasta'
```

or

```
input:
	os.path.join(resource_path, '/hg19/wholeGenome.fasta')
```

which can be less readable.

Decision: it makes sense to allow this for consistency purposes.


### trouble with sigil usage

Python 3.6 introduced formated string, which uses `{}`. It also accepts expressions (with addition of conversion `!` and format specification `:`)
but so I am very attempted to declear that all SoS string literals are python format string. However, the use of `{}` forces the doubling of regular
`{` and `}`, which is disastrous for scripts such as shell and R.

I guess I will have to specify that SoS string is different and has to be configurable sigil for different languages.

Rules:

1. Python string literals are acceptable (e.g. `'string'`, `"string"`, `'''string'''`, `"""string"""`.
2. Single quote versions are automatically treated as `r` string literals. E.g. `'\n'` is `r'\n'` in Python form.
  This format is generally recommended.
3. String interpolation with configurable SoS sigil is handled before the underlying function. Further processing
  by underlying Python expression is allowed.

  
### backward dependency rules?

There can be a need for dependency rules. For example, if a bam index file is needed (dependent upon) and
the bam file exists, then `samtools index` would be automatically called. This does not sound like a good
idea because `samtools index` can always be put as a regular part of the workflow. On the other hand, adding
such rules can help if these common steps are not always needed, or needed for multiple steps of the pipeline.
I am not sure how useful such magic is. Implementation wise it might not be too difficult, some syntax
like the following could be used

```
# section name does not have index so it would be called if and only if the 
# step is needed. pattern of output is defined as section option.
[ index_bam : output=*.bam.bai]

# input is defined from output
input = output[0][:-4]

# the action is defined as usual, but output does not need to be defined
# again.
run('samtools index ${input}')

```

could potentially work. This allows gnumake style definition of pipelines that can be used to construct complete workflows. It seems to contradict 
the design of SoS though.

Decision: auxiliary steps.


### Moving parameters of action to section level

Currently we have

```python
[step]

action( ... 
	wordir= ..., 
	resource= ...,
)
```

It might be easier to move them out side of the function

```python
[step]

workdir:
resource:

action( ... 
)
```

The major problem here is that a step can execute multiple actions, loop or sequential action
but some options might be action specific (e.g. resource)

```python
run('command1'), run('command2')
```

Decided to move output to the step level. Other options will be added one by one to either action
or step level.


### Alternative or configurable global sigil?

It is very tempting to use `{ }` as default sigil because this is what snakemake uses and what
python format string uses. However, because SoS tends to include scripts in R etc, `{}` would
be a bad choice.

It is easy enough to allow a global variable `default_sigil` to change default sigil to a user
specified one, but it seems to be an overkill.

Decision: not worry about it now because default sigil should be good enough for most cases.



### Section option `no_input`

This can be replaced by something like `input: None` so no separate option is needed.

Decided to use `input: None`


### Section option `terminal` and `starting`

These can be determined automatically if input and output files are specified, right? 
These options are therefore not needed. In addition, adding the option `terminal` does not prevent
other steps to depend on the ouput files of a step, leading to potential errors.

Decided to use no special option and rely on input and output specification of steps.


### Use of dictionary variable

It is easy to implement and can be very useful, but I am not sure if we should complicate SoS with another dictionary type. A potentially useful user case is to replace `labels`, e.g.

```python
sample_name = {x:os.path.basename(y) for x in x in sam_files}
```

and use the variable as

```python
${sample_name[input[0]]}
```

Decision: temporarily decided to add it because dictionary types can be really useful.



### Use `None` instead of `[]` for no input or output

Right now `input` and `output` directives can only accept a list of filenames so `[]` is used
for no input. Using `input: None` can potentially be more readable.

Decide: None is reserved for unknown input so `[]` should be used for no input.


### Dynamic output

If the output of step is dynamically determined, for example, by running `glob.glob` from output directory, the output
might be empty or wrong at the planning stage. For example

```python
[100]
output: glob.glob('*.bam')

[200]

run('samtools sort')
```

Step 200 might get an empty list at the planning stage and will be treated as a leaf/starting step with no input, and 
causes an error.

Solution:

* A `dynamic` section option to tell SoS that the output of step 100 is dynamic?
* Holding off executing step 200 if it does not have explicit input instead of relying
  on output of previous step?
* A `dynamic` property to related variables.

Decision: a dynamic function is added to SoS expressions. This essentially converts

```python
dynamic(glob.glob('*.bam'))
```

to a DynamicExpression object with `expr="glob.glob('*.bam')"` with the following definition,

```
class DynamicExpression:
	def __init__(self, expr):
		self.expr = expr
	
	def __repr__(self):
		return __repr__(eval(self.expr))

```

and will evaluate the result each time when it is used.

