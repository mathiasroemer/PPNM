CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:, $(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(CODE)
TIME = time --output=$@ --append --format "$$nthreads %e %U"

N=2e8

default: out.times.svg out.for.txt out.for-threadlocal.txt

out.times.svg : out.times.txt
	echo '\
		set terminal svg background "white";\
		set out "$@" ;\
		set xlabel "number of threads" ;\
		set ylabel "real time" ;\
		set title "Running time as function of number of threads" ;\
		plot "$<" with lp ;\
		'| tee log.gpi | gnuplot
#lp is linepoints

out.times.txt: main.exe
	for nthreads in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do\
		$(TIME) mono main.exe -nthreads:$$nthreads -nterms:$(N);\
	done

out.for-threadlocal.txt : main-for-threadlocal.exe
	$(TIME) mono $< -nterms:$(N) >>$@

out.for.txt : main-for.exe
	$(TIME) mono $< -nterms:$(N) >>$@

main.exe : main.cs ; $(MKEXE)
main-for.exe : main-for.cs ; $(MKEXE)
main-for-threadlocal.exe : main-for-threadlocal.cs ; $(MKEXE)

test: out1 out2 out3 out4

out1: main.exe ; mono $< -nthreads:1 -nterms:2e8 > $@
out2: main.exe ; mono $< -nthreads:1 -nterms:2e8 > $@
out3: main.exe ; mono $< -nthreads:1 -nterms:2e8 > $@
out4: main.exe ; mono $< -nthreads:1 -nterms:2e8 > $@

test2: main.exe
	mono $< -nthreads:1 -nterms:2e8 > out1 & 
	mono $< -nthreads:1 -nterms:2e8 > out2 & 
	mono $< -nthreads:1 -nterms:2e8 > out3 & 
	mono $< -nthreads:1 -nterms:2e8 > out4 & 

clean:
	$(RM) *.exe *.dll [Oo]ut* [Ll]og*
