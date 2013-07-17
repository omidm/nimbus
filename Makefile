all: main
main: main.o app.o
	g++ main.o app.o
main.o: worker/main.cc
	g++ -c worker/main.cc -I .
app.o: application/fluid-simulation/app.cc
	g++ -c application/fluid-simulation/app.cc -I .

