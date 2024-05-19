libs = ../matlib

CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:, $(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(CODE)
TIME = time --output=$@ --append --format "%e %U"

default : Out.txt Out.PartA.svg Out.PartB.svg Out.PartC.svg

Out.txt Out.data.txt: main.exe
	mono $< 1>Out.txt 2>Out.data.txt


Out.PartA.svg : Out.data.txt
	echo '\
	set terminal svg background "white" ;\
	set out "$@" ;\
	set key left ;\
	set xlabel "Time" ;\
	set tics out ;\
	set grid ;\
	set title "ODE IVP" ;\
	plot [0:10]\
	"$<" using 1:2 index 0 with points pointtype 6 title "u´´ = -u" \
	,"$<" using 1:2 index 1 with points pointtype 6 title "LV prey" \
	,"$<" using 1:3 index 1 with points pointtype 6 title "LV predators" ;\
' | tee log.gpi | gnuplot

Out.PartB.svg : Out.data.txt
	echo '\
	set terminal svg background "white" ;\
	set out "$@" ;\
	set key left ;\
	set xlabel "x" ;\
	set ylabel "y" ;\
	set tics out ;\
	set grid ;\
	set title "Equatorial motion" ;\
	plot "$<" using ((1/$$2)*cos($$1)):((1/$$2)*sin($$1)) index 2 with lines title "{/Symbol e} = 0, Circular" \
	, "$<" using ((1/$$2)*cos($$1)):((1/$$2)*sin($$1)) index 3 with lines title "{/Symbol e} = 0, Elliptical"  \
	,  "$<" using ((1/$$2)*cos($$1)):((1/$$2)*sin($$1)) index 4 with lines title "{/Symbol e} = 0.01, Elliptical" ;\
' | tee log.gpi | gnuplot	

Out.PartC.svg : Out.data.txt
	echo '\
	set terminal svg background "white" ;\
	set out "$@" ;\
	set key left ;\
	set xlabel "x" ;\
	set ylabel "y" ;\
	set tics out ;\
	set grid ;\
	set title "Three-body problem - Movement over time" ;\
	plot "$<" using 2:3 index 5 with points pointtype 6 title "Body 1" \
	, "$<" using 4:5 index 5 with points pointtype 6 title "Body 2" \
	, "$<" using 6:7 index 5 with points pointtype 6 title "Body 3" ;\
' | tee log.gpi | gnuplot	


vector.dll : $(libs)/matrix/vector.cs ; $(MKLIB)
genlist.dll : $(libs)/genlist.cs ; $(MKLIB)

main.exe : main.cs vector.dll genlist.dll; $(MKEXE)

clean:
	$(RM) *.exe *.dll *.log [Oo]ut* [Ll]og* [Ee]rr
