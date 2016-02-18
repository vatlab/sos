## Cookbook

### Process paired files one by one

```
input:
    ${[x for x in bam_files if '_R1_' in x]}
    ${[x for x in bam_files if '_R2_' in x]}
    : group_by='pairs'

run('')

### Arrange input files by filenames

run('''STAR

    INPUT ${','.join([x for x in bam_files if '_R1_' in x]}
    ${','.join([x for x in bam_files if '_R2_' in x]}
''')