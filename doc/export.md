## Export a SoS script to a series of scripts

```
sos export script.sos [workflow:steps] [workflow options] [--portable] [--script_dir DIR]
```

This command exports steps of a SoS script into a series of scripts, and a master shell script to call them. 

## An example 
For example, for the SoS script defined in the [tutorial](../README.md),

```python
#!/usr/bin/env sos-runner
##fileformat=SOS1.0

# This script aligns raw reads of a control and a mutated sample 
# to the reference genome and compare the expression values
# of the samples at genes A, B and C.

[1: no_input]
# create a index for reference genome
run('''
STAR --runMode genomeGenerate --genomeFastaFile human38.fasta --genomeDir STAR_index
''', output='STAR_index/chrName.txt')
    
[2]
# align the reads to the reference genome
input:
	${cmd_input}

run('''
STAR --genomeDir STAR_index --outSAMtype BAM SortedByCoordinate  --readFilesIn ${input[0]}  \
    --quantMode GeneCounts --outFileNamePrefix aligned/control
STAR --genomeDir STAR_index --outSAMtype BAM SortedByCoordinate --readFilesIn ${input[1]}  \
    --quantMode GeneCounts --outFileNamePrefix aligned/mutated
''', output=['aligned/control.out.tab', 'aligned/mutated.out.tab'])

[3]
# compare expression values
R('''
control.count = read.table('${input[0]}')
mutated.count = read.table('${input[1]}')
# normalize, compare, output etc, ignored.
pdf('myfigure.pdf')
# plot results
dev.off()
''', output='myfigure.pdf')

```

Running command

```bash
sos export myscript.sos --input case.fasta ctrl.fasta --script_dir export
```

will generate four files

```
export/default_1.sh
export/default_2.sh
export/default_3.R
export/myscript.sh
```
with contents

```
!/usr/bin/env bash
# index reference genome
STAR --runMode genomeGenerate --genomeFastaFile human38.fasta --genomeDir STAR_index
```

```
!/usr/bin/env bash
STAR --genomeDir STAR_index --outSAMtype BAM SortedByCoordinate  --readFilesIn case.fasta \
    --quantMode GeneCounts --outFileNamePrefix aligned/control
STAR --genomeDir STAR_index --outSAMtype BAM SortedByCoordinate --readFilesIn ctrl.fasta  \
    --quantMode GeneCounts --outFileNamePrefix aligned/mutated
```

```
control.count = read.table('aligned/control.out.tab')
mutated.count = read.table('aligned/mutated.out.tab')
# normalize, compare, output etc, ignored.
pdf('myfigure.pdf')
# plot results
dev.off()
```

and

```
#!/usr/bin/env bash

bash default_1.sh
bash default_2.sh
Rscript default_3.R
```

repectively.

You can repeat the analysis by running

```bash
bash export/myscript.sh
```

although you will lose all the workflow control features of Sos.

## Export selected workflow and steps

Similar to the execution of SoS workflows, you can specify the workflow and steps you would like to export, for example you can export up to step 50 of the mouse workflow defined in ``myscript1.sos`` using command

```bash
sos export myscript1.sos mouse:-50
```

Please refer to [execution of SoS workflows](execution.md) for more details.

## Export portable workflows (needed feature?)

The workflow exported by default expands all SoS variables with specified input files and parameters. The scripts are guaranteed to be executed correctly but are not flexible enough to work with, for example, other input files. The `--portable` option tries not to expand all parameters and add themm as language-specific variables. For example, 

```bash
sos export myscript.sos --input case.fasta ctrl.fasta --portable --script_dir export
```

will generate four files

```
export/default_1.sh
export/default_2.sh
export/default_3.R
export/myscript.sh
```
with contents

```
!/usr/bin/env bash
# index reference genome
STAR --runMode genomeGenerate --genomeFastaFile human38.fasta --genomeDir STAR_index
```

```
!/usr/bin/env bash

STAR --genomeDir STAR_index --outSAMtype BAM SortedByCoordinate  --readFilesIn $1 \
    --quantMode GeneCounts --outFileNamePrefix aligned/control
STAR --genomeDir STAR_index --outSAMtype BAM SortedByCoordinate --readFilesIn $2 \
    --quantMode GeneCounts --outFileNamePrefix aligned/mutated
```

```
args <- commandArgs()

control.count = read.table(args[1])
mutated.count = read.table(args[2])
# normalize, compare, output etc, ignored.
pdf('myfigure.pdf')
# plot results
dev.off()
```

and

```
#!/usr/bin/env bash

bash default_1.sh 
bash default_2.sh case.fasta ctrl.fasta
Rscript default_3.R aligned/control.out.tab aligned/mutated.out.tab
```

The generated scripts usually works but variables such as `${len(input)}` (number of input files) will be evaluated and might not work for other input files. Because there is no easy way for SoS to translate complex Python expressions to native languages, scripts exported with option `--export` might need some additional modifications to execute sucessfully. 
