libs = ../matlib

CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:, $(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(CODE)
TIME = time --output=$@ --append --format "%e %U"

default : Out.txt Out.Higgs.svg

Out.txt Out.data.txt : main.exe higgs.data.txt
	mono $< < higgs.data.txt 1> Out.txt 2> Out.data.txt


Out.Higgs.svg : higgs.data.txt Out.data.txt
	echo '\
	set terminal svg size 1200,700 background "white" ;\
	set out "$@" ;\
	unset key ;\
	set key right top ;\
	set xlabel "Energy [GeV]" ;\
	set ylabel "Signal [Certain units]" ;\
	set tics out ;\
	set grid ;\
	set title "Higgs boson data fitted to the Breit-Wigner function" ;\
	plot [90:170] \
	"$<" with errorbars pointtype 6 lc rgb "red" title "Data" \
	,"Out.data.txt" index 0 using 1:2 with lines lt 8 lw 2 dashtype 2 title "Breit-Wigner fit" ;\
	' | tee log.gpi | gnuplot	

matrix.dll : ../linear_equations/main.cs $(libs)/matrix/matrix.cs $(libs)/matrix/vector.cs $(libs)/genlist.cs ; $(MKLIB)

main.exe : main.cs matrix.dll ; $(MKEXE)

clean:
	$(RM) *.exe *.dll *.log [Oo]ut* [Ll]og* [Ee]rr
