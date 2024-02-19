CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:, $(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(CODE)

default: Out.lngamma.gnuplot.svg Out.gamma.gnuplot.svg Out.erf.gnuplot.svg

Out.lngamma.gnuplot.svg : Out.lngamma.txt Makefile
	echo '\
	set terminal svg background "white" ;\
	set out "$@" ;\
	set key left ;\
	set xlabel "x" ;\
	set ylabel "lngamma(x)" ;\
	set tics out ;\
	set grid ;\
	set title "LnGamma function" ;\
	plot [0:5]\
	"$<" with lines notitle \
	,"gamma.data.txt" using 1:(log($$2)) with points pointtype 6 notitle ;\
' | tee log.gpi | gnuplot


Out.gamma.gnuplot.svg : Out.gamma.txt Makefile
	echo '\
	set terminal svg background "white" ;\
	set out "$@" ;\
	set key left ;\
	set xlabel "x" ;\
	set ylabel "gamma(x)" ;\
	set tics out ;\
	set grid ;\
	set title "Gamma function" ;\
	plot [0:5]\
	"$<" with lines notitle \
	,"gamma.data.txt" using 1:2 with points pointtype 6 notitle ;\
' | tee log.gpi | gnuplot

Out.erf.gnuplot.svg : Out.erf.txt Makefile
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
	"$<" with lines notitle \
	,"erf.data.txt" using 1:2 with points pointtype 6 notitle \
	,"erf.data.txt" using (-($$1)):(-($$2)) with points pointtype 6 notitle ;\
' | tee log.gpi | gnuplot

Out.erf.txt: mainerf.exe
	mono $< > $@

Out.gamma.txt: maingamma.exe
	mono $< > $@

Out.lngamma.txt: mainlngamma.exe
	mono $< > $@

sfuns.dll : sfuns.cs ; $(MKLIB)

mainerf.exe : mainerf.cs sfuns.dll ; $(MKEXE)
	
maingamma.exe : maingamma.cs sfuns.dll ; $(MKEXE)

mainlngamma.exe : mainlngamma.cs sfuns.dll ; $(MKEXE)

clean:
	$(RM) *.exe *.dll [Oo]ut* [Ll]og*
