CC     = gcc
CFLAGS = -I./ -O3 -Wall `pkg-config --cflags gtk+-3.0 epoxy`


all:
	$(CC) $(CFLAGS) src/pd_io.c src/main.c -o pd `pkg-config --libs gtk+-3.0 epoxy`
