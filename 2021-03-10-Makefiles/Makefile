all: datos.pdf

# implicit rules
# variables
# $^ : todas las depencias
# $< : Solamente la primera dependencia
# $@ : target

CXX_FLAGS = -O2 -g

%.x : %.cpp
	g++ $(CXX_FLAGS) $< -o $@

datos.txt: updown.x
	./$< 1000 > $@

datos.pdf: updown.gp datos.txt
	gnuplot $<
	open $@

clean:
	rm -f datos.txt updown.x

clean-all:
	rm -f datos.txt updown.x datos.pdf
