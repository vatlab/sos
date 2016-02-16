## Command `sos-runner`

Command `sos-runner` is a short cut for ``sos run`` so

```bash
sos-runner myscript.sos
```

is equivalent to

```bash
sos run myscript.sos
```

The value of this alias is to allow the execution of SoS scripts directly if it has shebang line

```
#!/usr/bin/env sos-runner

```

and appropriate permissions set.


## Command `sos`

Command `sos` accept a number of subcommands (similar to `svn`, `git` etc). Its syntax follows

```bash
sos subcommand [subcommand-options]
```

### subcommand `run` (`sos-runner`)

```bash
sos run specfile [workflow:steps] [workflow-options] [--dryrun]
```

### subcommand `export`

```bash
sos export specfile [workflow:steps] [workflow-options] [--portable]
```

### subcommand `show`

```bash
sos show specfile [workflow:steps] 
```

Display details of all or specified workflow defined in `specfile`

### subcommand `admin` (set option??)


### subcommand `edit` (GUI editor??)

