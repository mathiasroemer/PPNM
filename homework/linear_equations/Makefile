libs = ../matlib/matrix

CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:, $(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(CODE)
TIME = time --output=$@ --append --format "%e %U"

default : Out.txt Out.n3fit.svg


Out.txt : main.exe
	mono $< > $@

Out.n3fit.svg : out.times.data
	echo '\
		set terminal svg background "white" ;\
		set out "$@" ;\
		set xlabel "Size N of matrix";\
		set ylabel "Elapsed time (s)";\
		set grid ;\
		set title "Elapsed time as a function of N";\
		f(x)=A*x+B;\
		fit f(x) "$<" via A,B;\
		plot "$<" with points pointtype 6 title "Data" \
		,f(x) title "A*N^3 + C" \
	'|tee log.gpi |gnuplot

out.times.data : time.exe
	>$@
	for N in $$(seq 100 20 200); do \
		time --format "$$N %e" --output $@ --append \
		mono $< -size:$$N 1>out 2>err ;\
	done 

matrix.dll : $(libs)/matrix.cs $(libs)/vector.cs ; $(MKLIB)

main.exe : main.cs matrix.dll ; $(MKEXE)
time.exe : time.cs matrix.dll ; $(MKEXE)

clean:
	$(RM) *.exe *.dll *.log [Oo]ut* [Ll]og* [Ee]rr
