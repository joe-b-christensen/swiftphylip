CC ?= gcc
CFLAGS ?= -O2 -g -std=gnu11 -Wall -Wextra -Wno-unused-parameter
LDLIBS = -lz -ldl -lpthread

build/swiftphylip: swiftphylip.c
	mkdir -p build
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -rf build
