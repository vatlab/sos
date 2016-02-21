<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Problems and ideas](#problems-and-ideas)
  - [Support for docker](#support-for-docker)
  - [Runtime control](#runtime-control)
  - [Resource control](#resource-control)
  - [Libraries](#libraries)
  - [Nested workflow](#nested-workflow)
  - [Requirement of steps (solved)](#requirement-of-steps-solved)
  - [Handling of filenames with spaces and other special characters (solved)](#handling-of-filenames-with-spaces-and-other-special-characters-solved)
  - [Section option `concurrent`](#section-option-concurrent)
  - [Default parameter `--input`](#default-parameter---input)
  - [variable definition (decided to use Python syntax)](#variable-definition-decided-to-use-python-syntax)
  - [A more pythonic approach?](#a-more-pythonic-approach)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Problems and ideas

### Support for docker 

### Runtime control

Limiting the files that the directory can the step can write to??? 


### Resource control

Limiting or monitoring the RAM and CPU (cores) the step uses???


### Libraries

Libraries would be python modules with defined SoS actions, but how to maintain and import these modules require further investigation. Furthermore, extensive use of libraries somehow beats the purpose of SoS (readability) because libraries hide the details of actions.

### Nested workflow

Not sure how this would work, but in SoS one can certainly define a step with

```python
run('sos-runner anotherworkflow ${input}'
```

to execute another workflow.

### Requirement of steps (solved)

CWL has the requirement keyword which handles all sorts of requirement. SoS solves this by

* `depends:` item in pipeline step, which lists all dependent files that will go into step signature.
* `fail_if()` action, which stops the workflow if the requirements are not satisfied.

### Handling of filenames with spaces and other special characters (solved)

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

if they expect filenames with special characters.

On the other hand,

```python
R('''
open("${input}")

''')
```

shoudl work correctly and it would be wrong if SoS mangles `${input}` during variable substitution.


### Section option `concurrent`

There are some other options to allow concurrent execution of step actions.

1. Default to concurrent but allow option `nonconcurrent`, because the actions should be safe 
  to execute in parallel most of the time (processing input files one by one or in pair, or
  with different options).

2. Do not use section option, but specify this in input parameters. E.g. `for_each` as
  concurrent for each, and `nc_for_each` for nonconcurrent for each. Perhaps 
  `for_all` as nonconcurrent?

### Default parameter `--input`

We do not have to allow a default parameter `--input`. It is easier to use 
but the parameter itself is not documented like other command line parameters. It might  
be better to force the definition of all parameters in the `[default]` section.

### variable definition (decided to use Python syntax)

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
specification. It might make sense to turn them all to python syntax.

