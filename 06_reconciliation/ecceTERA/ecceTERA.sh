mkdir standard_soxBYZ_cir
./ecceTERA_linux64 species.file=Cons_cir.chronogram gene.file=04_soxBYZ_standard.ufboot.rename amalgamate=1  compute.TD=false print.reconciliations=1 print.newick=1 output.dir=standard_soxBYZ_cir output.prefix=standard_soxBYZ_cir &>standard_soxBYZ.log

mkdir standard_soxAX_cir
./ecceTERA_linux64 species.file=Cons_cir.chronogram gene.file=04_soxAX_standard.ufboot.rename amalgamate=1  compute.TD=false print.reconciliations=1 print.newick=1 output.dir=standard_soxAX_cir output.prefix=standard_soxAX_cir &>standard_soxAX.log

mkdir standard_soxCD_cir
./ecceTERA_linux64 species.file=Cons_cir.chronogram gene.file=04_soxCD_standard.ufboot.rename amalgamate=1  compute.TD=false print.reconciliations=1 print.newick=1 output.dir=standard_soxCD_cir output.prefix=standard_soxCD_cir &>standard_soxCD.log

mkdir complex_soxBYZ_cir
./ecceTERA_linux64 species.file=Cons_cir.chronogram gene.file=04_soxBYZ_complex.ufboot.rename amalgamate=1  compute.TD=false print.reconciliations=1 print.newick=1 output.dir=complex_soxBYZ_cir output.prefix=complex_soxBYZ_cir &>complex_soxBYZ.log

mkdir complex_soxAX_cir
./ecceTERA_linux64 species.file=Cons_cir.chronogram gene.file=04_soxAX_complex.ufboot.rename amalgamate=1  compute.TD=false print.reconciliations=1 print.newick=1 output.dir=complex_soxAX_cir output.prefix=complex_soxAX_cir &>complex_soxAX.log

mkdir complex_soxCD_cir
./ecceTERA_linux64 species.file=Cons_cir.chronogram gene.file=04_soxCD_complex.ufboot.rename amalgamate=1  compute.TD=false print.reconciliations=1 print.newick=1 output.dir=complex_soxCD_cir output.prefix=complex_soxCD_cir &>complex_soxCD.log

mkdir standard_dsrAB_cir
./ecceTERA_linux64 species.file=Cons_cir.chronogram gene.file=04_dsrAB_standard.ufboot.rename amalgamate=1  compute.TD=false print.reconciliations=1 print.newick=1 output.dir=standard_dsrAB_cir output.prefix=standard_dsrAB_cir &>standard_dsrAB.log

mkdir complex_dsrAB_cir
./ecceTERA_linux64 species.file=Cons_cir.chronogram gene.file=04_dsrAB_complex.ufboot.rename amalgamate=1  compute.TD=false print.reconciliations=1 print.newick=1 output.dir=complex_dsrAB_cir output.prefix=complex_dsrAB_cir &>complex_dsrAB.log
