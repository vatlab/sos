" Vim syntax file
" Language:	SoS (extended from python.vim)
" Maintainer:	Bo Peng (bpeng@mdanderson.org)
"
" Usage
"
" copy to $HOME/.vim/syntax directory and add:
"
" au BufNewFile,BufRead SoS set syntax=sos
" au BufNewFile,BufRead *.sos set syntax=sos
"
" to your $HOME/.vimrc file
"
" force coloring in a vim session with:
"
" :set syntax=sos
"

" disable python builtin keyword
let python_no_builtin_highlight = 1

" load settings from system python
source $VIMRUNTIME/syntax/python.vim

" then add back EXCEPT for input, which is a SOS directive

syn keyword pythonBuiltin	False True None
syn keyword pythonBuiltin	NotImplemented Ellipsis __debug__
" built-in functions
syn keyword pythonBuiltin	abs all any bin bool chr classmethod
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
syn match sos_directive "^\(input\|output\|depends\|task|\w\+\)\s*:.*"
" match a line with non input/output/depends/task, and ends before a section
" or another directive
syn region script 
   \ start="^\(input\|output\|depends\|task\)\@!\w\+\s*:"ms=e+1,hs=e+1
   \ end="\(^\[\s*\w\+.*\]\s*$\|^\w\+\s*:\)"me=s-1,he=s-1,re=s-1

highlight sos_section_head guibg='Purple' gui=none
highlight sos_directive guifg='LightBlue' gui=none
highlight script guifg='Gray' gui=none


let b:current_syntax = "sos"

" vim:set tabstop=8 softtabstop=0 expandtab shiftwidth=4 smarttab:

