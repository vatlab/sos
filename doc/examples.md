<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Input specification](#input-specification)
  - [Process paired files one by one](#process-paired-files-one-by-one)
  - [Arrange input files by filenames](#arrange-input-files-by-filenames)
- [Command `sos export`](#command-sos-export)
  - [Export a SoS script to a series of scripts](#export-a-sos-script-to-a-series-of-scripts)
  - [Export selected workflow and steps](#export-selected-workflow-and-steps)
  - [Export portable workflows (TBD feature)](#export-portable-workflows-tbd-feature)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->


## Input specification
### Process paired files one by one

```
input:
    [x for x in bam_files if '_R1_' in x],
    [x for x in bam_files if '_R2_' in x],
    group_by='pairs'

run('')
```

### Arrange input files by filenames

```
run('''STAR

    INPUT ${','.join([x for x in bam_files if '_R1_' in x]}
    ${','.join([x for x in bam_files if '_R2_' in x]}
''')
```

## Command `sos export`

### Export a SoS script to a series of scripts

```
sos export script.sos [workflow:steps] [workflow options] [--portable] [--script_dir DIR]
```

This command exports steps of a SoS script into a series of scripts, and a master shell script to call them. 

For example, for the SoS script defined in the [tutorial](../README.md),

```python
#!/usr/bin/env sos-runner
##fileformat=SOS1.0

# This script aligns raw reads of a control and a mutated sample 
# to the reference genome and compare the expression values
# of the samples at genes A, B and C.


[parameters]
# Two input files in .fasta formats. The first one for control sample
# and the second one for mutated sample.
fasta_files=['control.fasta', 'mutated.fasta']

[1]
# create a index for reference genome
run('''
STAR --runMode genomeGenerate --genomeFastaFile human38.fasta --genomeDir STAR_index
''', output='STAR_index/chrName.txt')

ref_index=step_output
    
[2]
sample_type=['control', 'mutated']

# align the reads to the reference genome
input:
	fasta_files,
	group_by='single', for_each='sample_type'

depends:
   ref_index

output:
	['aligned/control.out.tab', 'aligned/mutated.out.tab']

run('''
STAR --genomeDir STAR_index --outSAMtype BAM SortedByCoordinate  --readFilesIn ${input[0]}  \
    --quantMode GeneCounts --outFileNamePrefix aligned/${_sample_type}
''')

[3]
# compare expression values
output:	'myfigure.pdf'

R('''
control.count <- read.table('${input[0]}')
mutated.count <- read.table('${input[1]}')
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

export/default_1.sh:

```
!/usr/bin/env bash
# index reference genome
STAR --runMode genomeGenerate --genomeFastaFile human38.fasta --genomeDir STAR_index
```

export/default_2.sh:

```
!/usr/bin/env bash
STAR --genomeDir STAR_index --outSAMtype BAM SortedByCoordinate  --readFilesIn case.fasta \
    --quantMode GeneCounts --outFileNamePrefix aligned/control
STAR --genomeDir STAR_index --outSAMtype BAM SortedByCoordinate --readFilesIn ctrl.fasta  \
    --quantMode GeneCounts --outFileNamePrefix aligned/mutated
```

export/default_3.R:

```
control.count <- read.table('aligned/control.out.tab')
mutated.count <- read.table('aligned/mutated.out.tab')
# normalize, compare, output etc, ignored.
pdf('myfigure.pdf')
# plot results
dev.off()
```

and

export/myscript.sh:

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

### Export selected workflow and steps

Similar to the execution of SoS workflows, you can specify the workflow and steps you would like to export, for example you can export up to step 50 of the mouse workflow defined in ``myscript1.sos`` using command

```bash
sos export myscript1.sos mouse:-50
```

Please refer to [execution of SoS workflows](execution.md) for more details.

### Export portable workflows (TBD feature)

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

export/default_1.sh:

```
!/usr/bin/env bash
# index reference genome
STAR --runMode genomeGenerate --genomeFastaFile human38.fasta --genomeDir STAR_index
```

export/default_2.sh:

```
!/usr/bin/env bash

STAR --genomeDir STAR_index --outSAMtype BAM SortedByCoordinate  --readFilesIn $1 \
    --quantMode GeneCounts --outFileNamePrefix aligned/control
STAR --genomeDir STAR_index --outSAMtype BAM SortedByCoordinate --readFilesIn $2 \
    --quantMode GeneCounts --outFileNamePrefix aligned/mutated
```

export/default_3.R:

```
args <- commandArgs()

control.count <- read.table(args[1])
mutated.count <- read.table(args[2])
# normalize, compare, output etc, ignored.
pdf('myfigure.pdf')
# plot results
dev.off()
```

and

export/myscript.sh:

```
#!/usr/bin/env bash

bash default_1.sh 
bash default_2.sh case.fasta ctrl.fasta
Rscript default_3.R aligned/control.out.tab aligned/mutated.out.tab
```

The generated scripts usually works but variables such as `${len(input)}` (number of input files) will be evaluated and might not work for other input files. Because there is no easy way for SoS to translate complex Python expressions to native languages, scripts exported with option `--export` might need some additional modifications to execute sucessfully. 
