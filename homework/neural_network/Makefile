libs = ../matlib

CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:, $(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(CODE)
TIME = time --output=$@ --append --format "%e %U"

default : Out.txt Out.data.svg 

Out.txt Out.data.txt : main.exe
	mono $< 1> Out.txt 2> Out.data.txt

Out.data.svg : Out.data.txt
	echo '\
	set terminal svg size 800,500 background "white" ;\
	set out "$@" ;\
	set key left ;\
	set xlabel "x" ;\
	set ylabel "y" ;\
	set tics out ;\
	set grid ;\
	set title "Network vs training data" ;\
	plot [-1:1]\
	"$<" index 0 using 1:2 with lines title "Linear interpolation of network values" \
	,"$<" index 1 using 1:2 with points pointtype 6 lc rgb "red" title "Training points" ;\
' | tee log.gpi | gnuplot	

extra.dll : $(libs)/matrix/matrix.cs $(libs)/matrix/vector.cs minimization.cs QRGS.cs; $(MKLIB)

main.exe : main.cs extra.dll; $(MKEXE)

clean:
	$(RM) *.exe *.dll *.log [Oo]ut* [Ll]og* [Ee]rr
