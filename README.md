# Script of Scripts (SoS)
Script of Scripts (SoS) is a lightweight workflow system that helps you turn your scripts in shell, R, Python, Perl, and other languages into readable pipelines that can be easily understood and modified by others. 

## A trial example
In its simplest form, a sos script consits of a series of scripts that can be executed sequentially by different intepreters.

Let us assume that you are an bioinformaticist that need to compare the expression levels of two samples for some genes of interest. After reading some online tutorial, you ends up with some working commands such as

```shell
STAR --runMode genomeGenerate --genomeFastaFile human38.fasta \
    --genomeDir STAR_index
STAR --genomeDir STAR_index --outSAMtype BAM SortedByCoordinate \
    --readFilesIn control.fasta  --quantMode GeneCounts \
    --outFileNamePrefix aligned/control
STAR --genomeDir STAR_index --outSAMtype BAM SortedByCoordinate \
    --readFilesIn mutated.fasta  --quantMode GeneCounts \
    --outFileNamePrefix aligned/mutated
```

The first command build an index, the second command aligns reads from the first sample, and the third command aligns reads from the second sample. Do not panic if you do not know what these commands are doing, this is just an example. 

These commands generates, among other files, two files named ``aligned/control.out.tab`` and ``aligned/mutated.out.tab`` with expression counts of all genes.  You then write a simple ``R`` script to analyze the results, something like

```R
control.count = read.table('aligned/control.out.tab')
mutated.count = read.table('aligned/mutated.out.tab')
# normalize, compare, output etc, ignored.
```

The real analysis is much more complicated than this but you get the idea. Instead of having two files lying around with perhaps another ``README`` file to describe what they are and how they are used, you can write a single SoS script named ``myanalysis.sos``

```python
##fileformat=SOS1.0

# This script aligns raw reads of a control and a mutated sample 
# to the reference genome and compare the expression values
# of the samples at genes A, B and C.

[1]
# create a index for reference genome
run('''
STAR --runMode genomeGenerate --genomeFastaFile human38.fasta \
    --genomeDir STAR_index
    ''')
    
[2]
# align the reads to the reference genome
run('''
STAR --genomeDir STAR_index --outSAMtype BAM SortedByCoordinate \
    --readFilesIn control.fasta  --quantMode GeneCounts \
    --outFileNamePrefix aligned/control
STAR --genomeDir STAR_index --outSAMtype BAM SortedByCoordinate \
    --readFilesIn mutated.fasta  --quantMode GeneCounts \
    --outFileNamePrefix aligned/mutated
''')

[3]
@ compare expression values
R('''
control.count = read.table('aligned/control.out.tab')
mutated.count = read.table('aligned/mutated.out.tab')
# normalize, compare, output etc, ignored.
''')

```

If you run command

```shell
sos run myanalysis.sos
```

It will execute a shell command, a shell script, and a R script in three steps. 

## Making the script work for other input files
After a while, before you almost forgot about this analysis, you need to analyze another pair of samples. You can of course save ``myanalysis.sos`` to ``myanalysis2.sos``, change filenames and run it. However, an easier way is to change your sos file to accommodate other files. This can be done by 

```python
##fileformat=SOS1.0

# This script aligns raw reads of a control and a mutated sample 
# to the reference genome and compare the expression values
# of the samples at genes A, B and C.

[1]
# create a index for reference genome
run('''
STAR --runMode genomeGenerate --genomeFastaFile human38.fasta \
    --genomeDir STAR_index
    ''')
    
[2]
# align the reads to the reference genome
run('''
STAR --genomeDir STAR_index --outSAMtype BAM SortedByCoordinate \
    --readFilesIn ${cmd_input[0]}  --quantMode GeneCounts \
    --outFileNamePrefix aligned/control
STAR --genomeDir STAR_index --outSAMtype BAM SortedByCoordinate \
    --readFilesIn ${cmd_input[1]}  --quantMode GeneCounts \
    --outFileNamePrefix aligned/mutated
''')

[3]
@ compare expression values
R('''
control.count = read.table('aligned/control.out.tab')
mutated.count = read.table('aligned/mutated.out.tab')
# normalize, compare, output etc, ignored.
''')

```

and execute the script as

```shell
sos run myanalysis.sos --input control1.fasta control2.fasta
```

Here the actual filenames are replaced by *SoS variables*. Basically, command line parameters are passed to SoS as variable ``cmd_input`` (a Python list with value ``['control.fasta', 'control2.fasta']`` in this example), and you can use the variable in ``${ }`` format. SoS variables are very powerful in that you can define and use many variables, and you can use any Python expression in ``${}``. This might sounds scary if you do not know Python well, but trust me, it is easy.

## Control the execution of SoS script

The script currently executes sequentially and runs all commands as if you enter them from command line (or within R). However, using some simple rules, you can convert your script to a real pipeline. Here is how the magic happens:

```python
##fileformat=SOS1.0

# This script aligns raw reads of a control and a mutated sample 
# to the reference genome and compare the expression values
# of the samples at genes A, B and C.

[1: no_input]
# create a index for reference genome
run('''
STAR --runMode genomeGenerate --genomeFastaFile human38.fasta \
    --genomeDir STAR_index
    ''', output='STAR_index/chrName.txt')
    
[2]
# align the reads to the reference genome
run('''
STAR --genomeDir STAR_index --outSAMtype BAM SortedByCoordinate \
    --readFilesIn control.fasta  --quantMode GeneCounts \
    --outFileNamePrefix aligned/control
STAR --genomeDir STAR_index --outSAMtype BAM SortedByCoordinate \
    --readFilesIn mutated.fasta  --quantMode GeneCounts \
    --outFileNamePrefix aligned/mutated
''', output=['aligned/control.out.tab', 'aligned/mutated.out.tab'])

[3]
@ compare expression values
R('''
control.count = read.table('${input[0]}')
mutated.count = read.table('${input[1]}')
# normalize, compare, output etc, ignored.
''')

```

Here we
 
- Use *step option* ``no_input`` to tell SoS the first step does not need any input.
- Use ``output='STAR_index/chrName.txt'`` to tell SoS the expected output of step 1.
- Use '' output=['aligned/control.out.tab', 'aligned/mutated.out.tab']`` to indicate the expected output of step 2.
- Use ``${input[0]}`` and ``${input[1]}`` to present the input of step 3, which is actually the output of step 2.

Now, when you run the same command
```shell
sos run myanalysis.sos --input control1.fasta control2.fasta
```

SoS will

- Ignore step 1 because this step has been run before with expected output.
- Ignore step 2 if you run the script multiple times because the output file already exist. However, this step will be re-run if the content of input files change, or if you change the command options in the script.
- Run step 2 in parallel if you specify a single command and tell SoS to apply the step to two input files. More on this later. 

These examples show roughly how SoS works. In the later sections, I will explain features of SoS in detail and demonstrates how to write powerful, yet readable scripts in SoS.