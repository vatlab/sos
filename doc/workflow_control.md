<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Sequential execution mode](#sequential-execution-mode)
- [Parallel execution in sequential mode](#parallel-execution-in-sequential-mode)
- [Parallel execution of SoS scripts](#parallel-execution-of-sos-scripts)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Sequential execution mode

By default, a SoS script is executed sequentially with steps executed on numeric order. This is the case when
you use SoS as a container of scripts and do not inform SoS input or output of steps.

## Parallel execution in sequential mode

Limited parallel execution is supported if a step explicitly runs its action in background or in another machine
(e.g. through the submission of a job). The workflow still executes sequentially and proceed until the output of
background jobs is needed, in this case it will wait indefinitely as long as the background job is still running,
and resume the execution as soon as the files become available.

To execute SoS script in this mode, you would need to define the output of the step being executed in background,
and the input of steps that rely on the output. You will also need to specify the `submitter` parameter of the background
step (e.g. `sh {} &`) to execute the actions in background.

This technique is suited for time-consuming quality control steps that are not needed for the data analysis steps.


## Parallel execution of SoS scripts

If parameter `-j` (jobs) is specified with a number more than 1, the SoS script will be executed in parallel
up to the specified number of concurrent processes. There are three sources of concurrent execution:

1. **multiple branches of the execution tree**: If the SoS script has multiple starting points (defined by step options
	`no_input` and `starting`), these branches will be executed in parallel.

2. **concurrent actions of single steps**: If a single SoS step executes actions repeatedly (with input option 
	`for_each`) and can be executed concurrently (with step option `concurrent`), these actions will be 
	executed in parallel.

3. **no wait** execution of sequential steps: If a step specifies output files, SoS will not wait for the
   completion of the step, execute steps with explicit input (or `no_input`) and proceed until the output file
   is needed, in which case SoS will have to wait till the output file is available.

SoS takes an conservative approach and execute actions in parallel only when it is safe to do so. That is 
to say, if no step option, input, or output of a step is specified (e.g. using SoS as a container of scripts),
the step will only be started with the completion of the previous step, and will block the execution of the
workflow until it is completed. For a SoS script to be executely safely in parallel mode, you should

1. Use option `no_input` if the step does not rely on any input. (e.g. action `download`)
2. Use option `concurrent` if the step contains loop (input parameter `for_each`) and the looped
  actions can be executed in parallel (e.g. processing each input files).
3. Use option `blocking` if the step should not be run by multiple processes. (e.g. action `download`)
4. Specify `input` and `depends` files for each step so that they would not be executed
  without needed input or dependent files.
5. Specify `output` files so that SoS knows what files to expect from a step and wait for the
  completion of the step if necessary.

Failure to specify these options correctly will make your script executed in sequential mode, or
fail to execute in parallel mode.


