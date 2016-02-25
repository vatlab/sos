<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Design of SoS (in progress, comments welcome)](#design-of-sos-in-progress-comments-welcome)
- [Tutorial](#tutorial)
  - [Organize your scripts as a SoS script](#organize-your-scripts-as-a-sos-script)
  - [Make the script work for other input files](#make-the-script-work-for-other-input-files)
  - [Ignore steps that do not need to be rerun](#ignore-steps-that-do-not-need-to-be-rerun)
  - [Execute steps in parallel](#execute-steps-in-parallel)
- [Unique features and limitations](#unique-features-and-limitations)
- [Summary](#summary)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

**Script of Scripts (SoS)** is a lightweight workflow system that helps you turn your scripts in shell, R, Python, Perl, and other languages into readable pipelines that can be easily understood and modified by others. It is also an easy-to-use alternative to workflow systems such as [CWL](http://common-workflow-language.github.io/draft-3/) with an emphasis on readability. 

## Design of SoS (in progress, comments welcome)
The core of SoS has mostly been implemented in another project but we are re-designing and re-implementing it to make SoS more user-friendly and powerful. Your involvement and suggestions are very welcome.

* User interface: [commandline arguments](doc/commandline.md)
* Format of SoS script: [format specification](doc/sos_format_v1.md)
* Features:
	* [proposed actions](doc/actions.md)
	* [export SoS scripts](doc/export.md)
	* [define and use of variables](doc/variables.md)
	* [workflow control](doc/workflow_control.md)
* Implementation details:
	*  [runtime signature](doc/runtime_signature.md)
	*  [Pending problems](doc/problems.md)
* Cookbook: [recipes for common and uncommon scenarios](doc/cookbook.md)
 

## Tutorial

A SoS script consists of one or more scripts, comments and optional SoS directives. In its simplest form, a sos script is simply a series of scripts that can be executed sequentially by different intepreters.

Let us assume that you are a bioinformaticist needed to compare the expression levels between two samples. After reading some online tutorials, you ended up with some working commands

```bash
# index reference genome
STAR --runMode genomeGenerate --genomeFastaFile human38.fasta --genomeDir STAR_index
# align reads to the reference genome
STAR --genomeDir STAR_index --outSAMtype BAM SortedByCoordinate  --readFilesIn control.fasta \
    --quantMode GeneCounts --outFileNamePrefix aligned/control
STAR --genomeDir STAR_index --outSAMtype BAM SortedByCoordinate --readFilesIn mutated.fasta \
    --quantMode GeneCounts --outFileNamePrefix aligned/mutated
```

The first command builds an index of the reference genome, the second command aligns reads from the first sample to the reference genome, and the third command aligns reads from the second sample to the reference genome. Do not panic if you do not know what these commands are doing, this is just an example.

These commands generate, among other files, two files named ``aligned/control.out.tab`` and ``aligned/mutated.out.tab`` with expression counts of all genes. You then wrote a ``R`` script to analyze the results, something like

```R
control.count = read.table('aligned/control.out.tab')
mutated.count = read.table('aligned/mutated.out.tab')
# normalize, compare, output etc, ignored.
pdf('myfigure.pdf')
# plot results
dev.off()
```

### Organize your scripts as a SoS script
The project completed successfully and you needed to archive the scripts for later reference. Instead of having two files lying around with perhaps another ``README`` file to describe what you have done, you can write a single SoS script named ``myanalysis.sos``

```python
#!/usr/bin/env sos-runner
##fileformat=SOS1.0

# This script aligns raw reads of a control and a mutated sample 
# to the reference genome and compare the expression values
# of the samples at genes A, B and C.

[1]
# index reference genome
run('''
STAR --runMode genomeGenerate --genomeFastaFile human38.fasta --genomeDir STAR_index
''')
    
[2]
# align reads to the reference genome
run('''
STAR --genomeDir STAR_index --outSAMtype BAM SortedByCoordinate  --readFilesIn control.fasta \
    --quantMode GeneCounts --outFileNamePrefix aligned/control
STAR --genomeDir STAR_index --outSAMtype BAM SortedByCoordinate --readFilesIn mutated.fasta \
    --quantMode GeneCounts --outFileNamePrefix aligned/mutated
''')

[3]
# compare expression values
R('''
control.count = read.table('aligned/control.out.tab')
mutated.count = read.table('aligned/mutated.out.tab')
# normalize, compare, output etc, ignored.
pdf('myfigure.pdf')
# plot results
dev.off()
''')

```

You can execute two shell scripts and a R script defined in this SoS script sequentially by running command

```bash
sos run myanalysis.sos
```

or simply

```bash
myanalysis.sos
```

if you give `myanalyis.sos` executable permission (`chmod +x myanalysis.sos`). 

### Make the script work for other input files
After a while, before you almost forgot about this analysis, you needed to analyze another pair of samples.
You could copy ``myanalysis.sos`` to ``myanalysis2.sos``, change filenames and run it, but an easier way is to
change your SoS file to accommodate other input files. This can be done by
defining a command line argument and passing files name to a **SoS variable**:

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
# index reference genome
run('''
STAR --runMode genomeGenerate --genomeFastaFile human38.fasta --genomeDir STAR_index
''')
    
[2]
# align reads to the reference genome
run('''
STAR --genomeDir STAR_index --outSAMtype BAM SortedByCoordinate --readFilesIn ${fasta_files[0]}  \
    --quantMode GeneCounts --outFileNamePrefix aligned/control
STAR --genomeDir STAR_index --outSAMtype BAM SortedByCoordinate --readFilesIn ${fasta_files[1]}  \
    --quantMode GeneCounts --outFileNamePrefix aligned/mutated
''')

[3]
# compare expression values
R('''
control.count = read.table('aligned/control.out.tab')
mutated.count = read.table('aligned/mutated.out.tab')
# normalize, compare, output etc, ignored.
pdf('myfigure.pdf')
# plot results
dev.off()
''')

```

A command line argument `fasta_files` is defined in the `[parameters]`
section of the script. With this definition, you can pass two filenames to
variable `fasta_files` from command line

```bash
sos run myanalysis.sos --fasta_files control1.fasta control2.fasta
```

`${fasta_files[0]}` and `${fasta_files[1]}` in command `STAR --genomeDir
...`  are replaced with their values before the commands are executed.

### Ignore steps that do not need to be rerun

Although the SoS script now accepts command line arguments, it is still no more than a compilation of scripts and you immediately realized that it is a waste of time to execute the first command each time. To solve this problem, you can convert the SoS script to a real workflow by telling SoS some more details of the commands:

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


[1: no_input]
# create a index for reference genome
run('''
STAR --runMode genomeGenerate --genomeFastaFile human38.fasta --genomeDir STAR_index
''', output='STAR_index/chrName.txt')
    
[2]
# align the reads to the reference genome
input:
	fasta_files

depends:
	'STAR_index/chrName.txt'

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

Here we
 
- Use **step option** ``no_input`` to tell SoS the first step does not need any input.
- Use `output='STAR_index/chrName.txt'` to specify the expected output of step 1.
- Use **input directive** to specify the input of step 2.
- Use **depends directive** to let step 2 depend on the output of step 1
- Use `${input[0]}` and `${input[1]` to use whatever files specified by the input directive. This is not required for this particular example but it makes the script a bit more general.
- Use `output=['aligned/control.out.tab', 'aligned/mutated.out.tab']` to indicate the expected output of step 2.
- Use ``${input[0]}`` and ``${input[1]}`` to present the input of step 3, which is the output of step 2. This effectively *connects* step 2 and step 3.

With such information, when you run the same command

```bash
sos run myanalysis.sos --input control1.fasta control2.fasta
```

SoS will ignore step 1 if this step has been run with output `STAR_index/chrName.txt`. The same happens to step 2 and 3 so all steps will be ignored if you run the script repeatedly with the same input and processing scripts. SoS uses **runtime signature** for each step and will re-run the step if and only if the content or filename of input files or the processing scripts are changed.

### Execute steps in parallel

SoS allows execution of workflows in parallel. Although the three steps of this example has to be
executed sequentially, the second step runs the `STAR` command twice on two input files, and can
be executed in parallel. You can tell SoS by modifying the script as follows

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

[1: no_input]
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
	group_by='single', labels='sample_type'

depends:
   ref_index

run('''
STAR --genomeDir STAR_index --outSAMtype BAM SortedByCoordinate  --readFilesIn ${input}  \
    --quantMode GeneCounts --outFileNamePrefix aligned/${_sample_type}
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

Here we 

1. Although not needed, save the output of the first step as `ref_index`
and use the variable instead of filename as a dependency of step 2. 
2. Use option `group_by='single'` to pass input one by one to action. The
action will be executed twice with `input` set to the first and second
input file respectively.
3. Define a variable `sample_type` and mark it as the labels for input
files (option `labels`). This will generate a variable `_sample_type` for
each input file so `_sample_type` will be `control` for the first input
file, and `mutated` for the second.

Now, if you execute the script with option `-j 2` (2 concurrent processes),

```bash
sos run myanalysis.sos --input control1.fasta control2.fasta -j 2
```

the second step would be run in parallel.

We have showed you four versions of the same SoS script, each using more features of SoS. This actually demonstrates one of the advantages of the SoS system, namely you can start using SoS in minutes without knowing any of its advanced features, and can gradually improve your script when needs arise.


## Unique features and limitations

The biggest difference between SoS and GNU make style workflow systems such as
[snakemake](https://bitbucket.org/johanneskoester/snakemake) is that SoS is
`input driven` instead of `output driven`. SoS will execute all steps in a
workflow even if some steps do not produce any output (but can for example
check input files and stop the execution if an error occurs). It can also
execute steps conditionally, repeatedly, or generate different outputs with
different input files. In contrast, make-style systems focus on output files
and ignore steps that are irrelevant to the final outcome.

Compared to XML or system-specific workflow lanaguges such as
[CWL](http://common-workflow-language.github.io/), SoS emphasizes greatly
on readability. The workflow steps are logically arranged (although actual
execution order might differ) and contain explict commands and scripts.
Most comments are significant and provide usable information in output such
as `sos show script`. This makes it easy for readers to read the workflow
and modify it if needed. This is very useful for fields such as
bioinformatics where workflows need to be frequently updated to accomodate
new tools, reference genomes, and annotations.

SoS is essentially a tool to compose and execute commands and scripts. It
uses workflow variables (and their derived forms) and string substitution
to create scripts (from user-provided templates) and execute them with
their own interpreters. Therefore, 

* SoS only supports command line tools. There is no plan to extend it for interactive or GUI usages.
* SoS does not understand the underlying tools and scripts and have no
control over exactly what they do. Incorrect use of workflow variables
might result in non-executable scripts.
* SoS workflow system is largely file based (although you can check environments such as environmental variables through SoS variables). It does not support features such as piping between steps (neither do systems such as CWL support it).

## Summary

The above example only shows a small fraction of what SoS can offer, but should be enough to demonstrate the unique features of SoS. Compared to maintaining multiple scripts or using more specifilized workflow systems such as [YAWL](http://www.yawlfoundation.org/), [CWL](http://common-workflow-language.github.io/), and [Galaxy](https://galaxyproject.org/),

* **SoS offers a way to organize your scripts in a single file**, which makes it easy to execute and maintain. You can include small and freqently changed commands and scripts in SoS and keep large and stable scripts as separate scripts.
* **SoS scripts are human readable and writable**. A SoS script consists of only commands, scripts, descriptions, and a small amount of SoS directives. It is highly readble and can be easily modified by users with little or no knowledge of SoS. This is especially important for fields such as bioinformatics where workflows need to be constantly changed to reflect new reference genomes, annotation sources, and new or newer versions of tools. In comparison, it would take a lot of time and practice to learn its syntax write a [CWL](http://common-workflow-language.github.io/) workflow (see [this CWL tutorial](https://github.com/common-workflow-language/workflows/wiki/Tutorial-DRAFT2) for an example). 
* **SoS help you execute your scripts with advanced workflow features**. The workflow features of SoS is easy to use yet very powerful in helping you execute your pipelines efficiently not only locally, but also on cluster and cloud systems. For example, using appropriate parameters, step 2 in the above example can be executed in parallel or be submitted to different computing nodes of a cluster system. step 3 will automatically start once step 2 is completed.

If you are afraid of being tied to a new workflow tool, rest assured that SoS allows you to **[export SoS scripts](doc/export.md)** to a series of scripts called by a master bash (or windows .bat) script. This would allow you to execute your workflow in an environement without SoS installed.

Please refer to the SoS documentation for more details and feel free to [contact me](mailto:ben.bob@gmail.com) if you have any comment on this project.

