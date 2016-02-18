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
