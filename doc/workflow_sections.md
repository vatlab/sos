<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Specification and execution of workflow](#specification-and-execution-of-workflow)
  - [Define a single workflow](#define-a-single-workflow)
  - [Define multiple workflows](#define-multiple-workflows)
  - [Shared steps between workflows](#shared-steps-between-workflows)
  - [Execution of a subset of steps](#execution-of-a-subset-of-steps)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

# Specification and execution of workflow

A SoS script can specify one or more workflows. Each workflow consists of one or more numbered steps. The numbers specify the **logical order** by which the steps are executed, but a later step might be executed before the completion of previous steps if it does not depend on the output of these steps.

## Define a single workflow

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

## Define multiple workflows 

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

## Shared steps between workflows

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

## Execution of a subset of steps

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
