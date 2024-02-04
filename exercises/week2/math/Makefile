Out.txt : math.exe
	mono math.exe > Out.txt

sfuns.dll : sfuns.cs
	mcs -target:library -out:sfuns.dll sfuns.cs

math.exe : math.cs sfuns.dll
	mcs -out:math.exe -target:exe -reference:sfuns.dll math.cs
	#mcs -out:math.exe -target:exe math.cs

.PHONEY : clean
clean:
	rm -f math.exe Out.txt sfuns.dll
