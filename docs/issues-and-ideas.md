
Issues and ideas:


-- Documenter:
using FinEtoolsAcoustics
using DocumenterTools
Travis.genkeys(user="PetrKryslUCSD", repo="https://github.com/PetrKryslUCSD/FinEtoolsAcoustics.jl")

using Pkg; Pkg.add("DocumenterTools");                                 
using DocumenterTools                                                  
DocumenterTools.genkeys(user="PetrKryslUCSD", repo="git@github.com:PetrKryslUCSD/FinEtoolsAcoustics.jl.git")                                                  
using Pkg; Pkg.rm("DocumenterTools");  

- How to go back to the code?
- 
{ // This build system simply opens a new interactive Julia REPL 
    // 2023 Petr Krysl
    // Depends on AbbreviatedStackTraces
    "target": "terminus_open",
    "auto_close": false,
    "title": "Julia REPL",
    "shell_cmd": "julia", // Assume that PATH is set before invoking the editor
    "cwd": "${file_path:${folder}}",
    "selector": "source.julia",
    "focus": true,
    // "file_regex": "^in expression starting at (\\S.+.jl):([0-9]+)\\s*"
    "file_regex": "^\\s*[@]\\s\\S.*\\s(\\S.+.jl):([0-9]+)\\s*"
    // "file_regex": "^\\s*[@]\\s\\S.*\\s(\\S.+.jl):([0-9]+)\\s*"
}
// @ C:\Users\pkonl\Documents\00WIP\GAexplorations.jl\example1.jl:15 


{ // This build system simply opens a new interactive Julia REPL 
    // 2023 Petr Krysl
    // Depends on AbbreviatedStackTraces
    "target": "terminus_open",
    "auto_close": false,
    "title": "Julia REPL",
    "shell_cmd": "julia", // Assume that PATH is set before invoking the editor
    "cwd": "${file_path:${folder}}",
    "selector": "source.julia",
    "focus": true,
    "file_regex": "^\\s*[@]\\s\\S.*\\s(\\S.+.jl):([0-9]+)\\s*"
}


def displaymatch(match):
    if match is None:
        return None
    return '<Match: %r, groups=%r>' % (match.group(), match.groups())

import re
valid = re.compile(r"^\s*[@]\s\S.*\s(\S.+.jl):([0-9]+)\s*")
displaymatch(valid.match(r" @ something C:\Users\pkonl\Documents\00WIP\GAexplorations.jl\example1.jl:15 "))
