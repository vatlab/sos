## Design decision 1:

Which is better?

```bash
sos run myscript.sos 
```

or 

```bash
sos myscript.sos --run 
```

The latter has the advantage of executing as

```bash
myscript --run
```

if the script is named ``myscript`` with appropriate shebang line

```
#!/usr/bin/env sos
```

## Design decision 2:

```bash
sos run myscript.sos mouse:2
```

or 

```bash
sos run myscript.sos --workflow mouse:2 
```

The latter is more verbose (or a rarely used case) and has the advantage (?) that allowing positional argument of ``--input``, namely

```bash
sos run myscript.sos file1 file2
```

instead of 

```bash
sos run myscript.sos --input file1 file2
```


