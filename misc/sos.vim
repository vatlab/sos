" Vim syntax file
" Language:	SoS (extended from python.vim)
" Maintainer:	Bo Peng (bpeng@mdanderson.org)
"
" Usage
"
" copy this file to $HOME/.vim/syntax directory and add
"
"     let sos_fold=1
"     autocmd BufNewFile,BufRead *.sos set syntax=sos foldcolumn=2
"
" to your $HOME/.vimrc (or .gvimrc) file.
"
" Remove sos_fold=1 line and foldcolumn=2 if you do not want to fold scripts.

" disable python builtin keyword
let python_no_builtin_highlight = 1

" load settings from system python
source $VIMRUNTIME/syntax/python.vim

" then add back EXCEPT for input, which is a SOS directive

syn keyword pythonBuiltin	False True None
syn keyword pythonBuiltin	NotImplemented Ellipsis __debug__
" built-in functions
syn keyword pythonBuiltin	abs all any bin bool chr classmethod include
syn keyword pythonBuiltin	compile complex delattr dict dir divmod
syn keyword pythonBuiltin	enumerate eval filter float format
syn keyword pythonBuiltin	frozenset getattr globals hasattr hash
syn keyword pythonBuiltin	help hex id int isinstance
syn keyword pythonBuiltin	issubclass iter len list locals map max
syn keyword pythonBuiltin	min next object oct open ord pow print
syn keyword pythonBuiltin	property range repr reversed round set
syn keyword pythonBuiltin	setattr slice sorted staticmethod str
syn keyword pythonBuiltin	sum super tuple type vars zip __import__
" Python 2.6 only
syn keyword pythonBuiltin	basestring callable cmp execfile file
syn keyword pythonBuiltin	long raw_input reduce reload unichr
syn keyword pythonBuiltin	unicode xrange
" Python 3.0 only
syn keyword pythonBuiltin	ascii bytearray bytes exec memoryview
" non-essential built-in functions; Python 2.6 only
syn keyword pythonBuiltin	apply buffer coerce intern
highlight link pythonBuiltin	Function

" SoS rule:
" directive starts with  "name:"
" section has the format of  "[section: options]"
" script has "action: ...."

" parenthetical part of def and class
syn match sos_section_head "^\[\s*\w\+.*\]\s*$"
syn match sos_directive "^\(input\|output\|depends\|parameter\|task\)\s*:"
syn match sos_preprocessor "^%\(if\|elif\|else\|endif\|cell\|include\|from\|set_options\).*$"
syn match sos_magic "^%\(with\|use\|set\|get\|put\|sandbox\|paste\|\|dict\|restart\|dict\|preview\|run\|rerun\).*$"

" match a line with non input/output/depends/task, and ends before a section
" or another directive
" we cannot put ms=e+1 to start because it will mess up folding...
" see http://stackoverflow.com/questions/8693721/vim-how-do-i-start-a-syntax-fold-on-the-line-after-a-regexp-match-python-func
" for details.
if exists("sos_fold")
  syn region script
     \ start="^\(input\|output\|depends\|parameter\|task\)\@!\w\+\s*:"hs=e+1
     \ end="^\S"me=s-1,he=s-1,re=s-1
     \ fold contains=TOP containedin=ALL
else
  syn region script
     \ start="^\(input\|output\|depends\|parameter\|task\)\@!\w\+\s*:"ms=e+1,hs=e+1
     \ end="^\S"me=s-1,he=s-1,re=s-1
     \ contains=TOP containedin=ALL
endif

highlight sos_section_head guibg='Purple' gui=none
highlight sos_directive guifg='LightBlue' gui=bold
highlight sos_preprocessor guifg='LightRed' gui=bold
highlight sos_magic guifg='Orange' gui=bold
highlight script guifg='Gray' gui=none

syn sync fromstart
setlocal foldmethod=syntax
setlocal foldminlines=4
" vim:set tabstop=8 softtabstop=0 expandtab shiftwidth=4 smarttab:

