CC ?= cc
CFLAGS ?= -std=c99 -Wall -Wextra -Werror -pedantic

TEST_BINARY := tests/test_ephemeris

.PHONY: all test clean

all: $(TEST_BINARY)

$(TEST_BINARY): src/ephemeris.c include/ephemeris.h tests/test_ephemeris.c
	$(CC) $(CFLAGS) -Iinclude src/ephemeris.c tests/test_ephemeris.c -o $@

test: $(TEST_BINARY)
	./$(TEST_BINARY)

clean:
	rm -f $(TEST_BINARY)
