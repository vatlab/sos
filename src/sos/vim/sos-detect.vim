function! s:isSoS()
	let shebang = getline(1)
	if shebang =~# '^#!.*/bin/env\s\+sos-runner\>' | return 1 | en
	if shebang =~# '^#!.*/bin/sos-runner\>' | return 1 | en
	return 0
endfunction

au BufRead,BufNewFile * if !did_filetype() && s:isSoS() | setf sos | en

