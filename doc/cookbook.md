<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Cookbook](#cookbook)
  - [Process paired files one by one](#process-paired-files-one-by-one)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

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