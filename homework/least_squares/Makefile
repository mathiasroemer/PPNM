libs = ../matlib/matrix

CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:, $(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(CODE)
TIME = time --output=$@ --append --format "%e %U"

default : Out.txt Out.fit.svg

Out.fit.svg : data.data
	echo '\
		set terminal svg background "white" ;\
		set out "$@" ;\
		set xlabel "Time t (days)";\
		set ylabel "Activity y of ThX (relative units)";\
		set grid ;\
		set title "Least squares fit of radioactive decay";\
		f(x) = a*exp(-l*x);\
		fp(x) = (a+da)*exp(-(l+dl)*x);\
		fm(x) = (a-da)*exp(-(l-dl)*x);\
		a=142.010145492731; l=0.171046463563344; da=1.21398361; dl=0.0237065392;\
		plot "$<" with errorbars pointtype 6 title "Data" \
		,f(x) title "a*exp({/Symbol l}*x)" \
		,fp(x) title "(a+da)*exp(-({/Symbol l}+d{/Symbol l})*x)" \
		,fm(x) title "(a-da)*exp(-({/Symbol l}-d{/Symbol l})*x)" \
	'|tee log.gpi |gnuplot


data.data :
	>$@
	echo "1 117 6">>$@
	echo "2 100 5">>$@
	echo "3 88 4">>$@
	echo "4 72 4">>$@
	echo "6 53 4">>$@
	echo "9 29.5 3">>$@
	echo "10 25.2 3">>$@
	echo "13 15.2 2">>$@
	echo "15 11.1 2">>$@

Out.txt: main.exe
	mono $< > $@

QRGS.dll : QRGS.cs $(libs)/matrix.cs $(libs)/vector.cs ; $(MKLIB)

main.exe : main.cs QRGS.dll ; $(MKEXE)

clean:
	$(RM) *.exe *.dll *.log [Oo]ut* [Ll]og* [Ee]rr
