## Format of SoS step


```
[name_step: option1, option2, ...]
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

### step input

See [step input](input_spec.md)

### step output and variables

See [use of SoS variables](variables.md)

### step actions

See [step actions](actions.md)



