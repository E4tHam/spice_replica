1. To use the parser in MATLAB, keep all exsisting files and benchmarks in the same directory. The command is:
2. [LINELEM, NLNELEM, INFO, NODES, LINNAME, NLNNAME, PRINTNV, PRINTBV, PRINTBI, PLOTNV, PLOTBV, PLOTBI] = parser(file);
   file is the file name of benchmark, for example: [LINELEM, NLNELEM, INFO, NODES, LINNAME, NLNNAME, PRINTNV, PRINTBV, PRINTBI, PLOTNV, PLOTBV, PLOTBI] = parser('./Benchmarks/rc_line.ckt');
3. The parser can work on a 64-bit windows system with MATLAB R2022b
