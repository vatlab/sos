## Export a SoS script to several scripts.

```
sos export myscript.sos --input filenames 
```

Limitations:

* The scripts might have a number of SoS variables that depends on your input. For example, an parameter might be `${len(input)}$` which is difficult to export without knowning the input. SoS can ensure that the exported scripts work for the provided parameters, but not more.


