libs = ../matlib/matrix

CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:, $(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(CODE)
TIME = time --output=$@ --append --format "%e %U"

default : Out.txt Out.erf.svg

Out.txt Out.erf.txt : main.exe
	mono $< 1>Out.txt 2>Out.erf.txt

Out.erf.svg : Out.erf.txt erf.data.txt
	echo '\
	set terminal svg background "white" ;\
	set out "$@" ;\
	set key left ;\
	set xlabel "x" ;\
	set ylabel "erf(x)" ;\
	set tics out ;\
	set grid ;\
	set title "Error function" ;\
	plot [-3:3]\
	"$<" with lines title "Erf integrator values" \
	,"erf.data.txt" using 1:2 with points pointtype 6 title "Erf exact values" ;\
' | tee log.gpi | gnuplot	


#vector.dll : $(libs)/vector.cs ; $(MKLIB)

main.exe : main.cs ; $(MKEXE)

clean:
	$(RM) *.exe *.dll *.log [Oo]ut* [Ll]og* [Ee]rr
