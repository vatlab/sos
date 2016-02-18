## Runtime signatures of SoS

### Content of runtime signature

Runtime signature of a SoS step consists of

1. Exact command or script
2. MD5 signature of all input files
3. MD5 signature of all dependent files
4. MD5 signature of all output files
5. (For information only): standard output from command execution
6. (For information only): standard error output from command execution

The completion of a SoS step will create a step signature. The step will be ignored if the signature is not changed, or re-run if the signature is changed when it is executed the second time.

### Storage of runtime signature

The name and location of a runtime signature is determined by the first output file of the step. 

There are two types of output files

1. If the output file is outside of the current working directory, it is a global output. Such output files are typically downloaded resource files and files (e.g. indexes) derived from them. These files are typically stored in a global resource directory.

2. If the ouput file is under the current working directory, it a local project-specific output.

The global signatures are stored under `$HOME/.sos/.runtime/absolute_path_to_file`. For example, if the output file is `/reference/hg18/wholegenome/STAR_index/chrname.txt`, the signature would be `$HOME/.sos/.runtime/reference/hg18/wholegenome/STAR_index/chrname.txt.exe_info`.

The local signatures are stored under the current working directory with name `.sos/.runtime/relative_path_to_file`. For example, if the output file is `aligned/accepted_hits.bam`, the signature file would be `.sos/.runtime/aligned/accepted_hits.bam.exe_info`.

### Portability of runtime signatures: discouraged

Runtime signatures are usually not portable because filenames and pathes would change when the project is moved to another location. It is possible to make the signature more portable, e.g., try to use relative path but there is no guarantee of portability of runtimes signatures. I would suggest that not making the portability of runtime signatures one of the design goals.

