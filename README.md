# Script of Scripts (SoS)
Script of Scripts (SoS) is a lightweight workflow system that helps you turn your scripts in shell, R, Python, Perl, and other languages into readable pipelines that can be easily understood and modified by others. 


## Status
The core of SoS has mostly been implemented in another project but we are re-designing and re-implementing it to make SoS more user-friendly and powerful. Your involvement and suggestions are very welcome.

* **in progress**: Command line interface: [commandline arguments](doc/commandline.md)
* **in progress**: Format of SoS script: [file format/workflow sections](doc/workflow_sections.md)

## A trivial example
In its simplest form, a sos script consits of a series of scripts that can be executed sequentially by different intepreters.

Let us assume that you are a bioinformaticist needed to compare the expression levels between two samples. After reading some online tutorials, you ended up with some working commands such as

```bash
# index reference genome
STAR --runMode genomeGenerate --genomeFastaFile human38.fasta --genomeDir STAR_index
# align reads to the reference genome
STAR --genomeDir STAR_index --outSAMtype BAM SortedByCoordinate  --readFilesIn control.fasta \
    --quantMode GeneCounts --outFileNamePrefix aligned/control
STAR --genomeDir STAR_index --outSAMtype BAM SortedByCoordinate --readFilesIn mutated.fasta \
    --quantMode GeneCounts --outFileNamePrefix aligned/mutated
```

The first command builds an index, the second command aligns reads from the first sample, and the third command aligns reads from the second sample. Do not panic if you do not know what these commands are doing, this is just an example. 

These commands generate, among other files, two files named ``aligned/control.out.tab`` and ``aligned/mutated.out.tab`` with expression counts of all genes. You then wrote a ``R`` script to analyze the results, something like

```R
control.count = read.table('aligned/control.out.tab')
mutated.count = read.table('aligned/mutated.out.tab')
# normalize, compare, output etc, ignored.
pdf('myfigure.pdf')
# plot results
dev.off()
```

The project completed successfully and you needed to archive the scripts for later reference. Instead of having two files lying around with perhaps another ``README`` file to describe what you have done, you can write a single SoS script named ``myanalysis.sos``

```python
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

If you run command

```bash
sos run myanalysis.sos
```

It will execute two shell scripts and a R script in three steps sequentially. 

## Making the script work for other input files
After a while, before you almost forgot about this analysis, you needed to analyze another pair of samples. You could copy ``myanalysis.sos`` to ``myanalysis2.sos``, change filenames and run it, but an easier way is to change your SoS file to accommodate other input files. This can be done by replacing input filenames in ``analysis.sos`` with **SoS variables** `${cmd_input}` (command line input):

```python
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
STAR --genomeDir STAR_index --outSAMtype BAM SortedByCoordinate  --readFilesIn ${cmd_input[0]}  \
    --quantMode GeneCounts --outFileNamePrefix aligned/control
STAR --genomeDir STAR_index --outSAMtype BAM SortedByCoordinate --readFilesIn ${cmd_input[1]}  \
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

and execute the script as

```bash
sos run myanalysis.sos --input control1.fasta control2.fasta
```

Basically, command line parameters are passed to SoS as variable `cmd_input`, which is a Python list with value `['control.fasta', 'control2.fasta']` in this example. Because these two files are processed separately, you use `${cmd_input[0])` and `${cmd_input[1]}` to return two filenames.

## Convert the SoS script to a real pipeline

Although the SoS script now accepts command line arguments, it is still no more than a compilation of scripts and you immediately realized that it is a waste of time to execute the first command each time. To solve this problem, you can convert the SoS script to a real pipeline by telling SoS some more details of the commands:

```python
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

Here we
 
- Use **step option** ``no_input`` to tell SoS the first step does not need any input.
- Use `output='STAR_index/chrName.txt'` to specify the expected output of step 1.
- Use **input directive** to specify the input of step 2.
- Use `${input[0]}` and `${input[1]` to use whatever files specified by the input directive. This is not required for this particular example but it makes the script a bit more general.
- Use `output=['aligned/control.out.tab', 'aligned/mutated.out.tab']` to indicate the expected output of step 2.
- Use ``${input[0]}`` and ``${input[1]}`` to present the input of step 3, which is the output of step 2. This effectively *connects* step 2 and step 3.

With such information, when you run the same command

```bash
sos run myanalysis.sos --input control1.fasta control2.fasta
```

SoS will ignore step 1 if this step has been run with output `STAR_index/chrName.txt`. The same happens to step 2 and 3 so all steps will be ignored if you run the script repeatedly with the same input and processing scripts. SoS uses **runtime signature** for each step and will re-run the step if and only if the content or filename of input files or the processing scripts are changed.

## Summary

The above example only shows a small fraction of what SoS can offer, but shouls be enough to demonstrate the unique features of SoS. Compared to maintaining multiple scripts or using more specifilized workflow systems such as [YAWL](http://www.yawlfoundation.org/), [CWL](http://common-workflow-language.github.io/), and [Galaxy](https://galaxyproject.org/),

* SoS organizes all scripts in a single file, which makes it easy to execute and maintain. In the meantime, SoS allows the creation of shared actions to organize common actions into modules.
* SoS script consists of only scripts, descriptions, and a small amount of SoS directives. It is highly readble and can be easily modified by you and others. This is especially important for fields such as bioinformatics where workflows need to be constantly changed to reflect new reference genomes, annotation sources, and new or newer versions of tools.
* The pipeline features of SoS is easy to use yet very powerful in helping you execute your pipelines efficiently not only locally, but on cluster and cloud systems. For example, using appropriate parameters, step 2 in the above example can be executed in parallel or be submitted to different computing nodes of a cluster system. step 3 will automatically start once step 2 is completed.

Please refer to the SoS documentation for more details and feel free to [contact me](mailto:ben.bob@gmail.com) if you have any comment on this project.


