Out1.txt Out2.txt Out3.txt: main.exe main2.exe main3.exe inputfile input.txt input3.txt
	mono main.exe $$(cat inputfile) > Out1.txt
	mono main2.exe < input.txt 2> Out2.txt
	mono main3.exe -input:input3.txt -output:Out3.txt

main.exe: main.cs
	mcs -target:exe -out:$@ $<

main2.exe: main2.cs
	mcs -target:exe -out:$@ $<

main3.exe: main3.cs
	mcs -target:exe -out:$@ $<

clean:
	rm -f main.exe Out1.txt main2.exe Out2.txt main3.exe Out3.txt
