<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [format of SoS scripts](#format-of-sos-scripts)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## format of SoS scripts

```python
#!/usr/bin/env sos-runner
##fileformat=SOS1.0
#
# Other overal comments such as
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

[default]
par1=default1
   comment1
   
par2=default2
   comment2

[step1]
step1 input, variables, and actions

[step2]
step2 input, variables, and actions
```
