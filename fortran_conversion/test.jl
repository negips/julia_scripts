println("Testing FortranTranspiler.jl")

#include("./FortranTranspiler.jl")

#using .FortranTranspiler

nekpath = "/home/prabal/workstation/git/Nek5000_playing_around/core/"

fname = "reader_re2"

rd = readdir(nekpath)

f  = open("$nekpath$fname.f",read=true)

code = read(f,String)
jcode = convert_fortran(code)
close(f)

f = open("$fname.jl",write=true)
write(f,jcode)
close(f)

