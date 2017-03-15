CC := clang

fast:
	$(CC) -O3 -ffast-math  mainTracing.c 000pixel.h 000pixel.o -lglfw -framework OpenGL
accurate:
	$(CC) -O3 mainTracing.c 000pixel.h 000pixel.o -lglfw -framework OpenGL
debug:
	$(CC) -g mainTracing.c 000pixel.h 000pixel.o -lglfw -framework OpenGL