<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Problems and ideas](#problems-and-ideas)
  - [Support for docker (still investigating)](#support-for-docker-still-investigating)
  - [Runtime control (still investigating)](#runtime-control-still-investigating)
  - [Resource control (still investigating)](#resource-control-still-investigating)
  - [Libraries (still investigating)](#libraries-still-investigating)
  - [Nested workflow (still investigating)](#nested-workflow-still-investigating)
  - [Requirement of steps (solved)](#requirement-of-steps-solved)
  - [Handling of filenames with spaces and other special characters (solved)](#handling-of-filenames-with-spaces-and-other-special-characters-solved)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Problems and ideas

### Support for docker (still investigating)

### Runtime control (still investigating)

Limiting the files that the directory can the step can write to??? 


### Resource control (still investigating)

Limiting or monitoring the RAM and CPU (cores) the step uses???


### Libraries (still investigating)

Libraries would be python modules with defined SoS actions, but how to maintain and import these modules require further investigation. Furthermore, extensive use of libraries somehow beats the purpose of SoS (readability) because libraries hide the details of actions.

### Nested workflow (still investigating)

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

