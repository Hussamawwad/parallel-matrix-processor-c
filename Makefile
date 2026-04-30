CC      = gcc
CFLAGS  = -std=c11 -Wall -Wextra -O2 -Iinclude -fopenmp
LDFLAGS = -lm -fopenmp

SRC = \
  src/main.c \
  src/matrix.c \
  src/util.c \
  src/ops_addsub.c \
  src/ops_mul.c \
  src/ops_det.c \
  src/ops_eigen.c

OBJ = $(SRC:.c=.o)
BIN = matrix_menu

.PHONY: all run clean debug


all: $(BIN)

$(BIN): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $@ $(LDFLAGS)


src/%.o: src/%.c
	$(CC) $(CFLAGS) -c -o $@ $<

run: $(BIN)
	./$(BIN) config.txt


debug: clean
	$(CC) -std=c11 -Wall -Wextra -g -O0 -Iinclude -fopenmp $(SRC) -o $(BIN) $(LDFLAGS)
	@echo "✅ Build finished with debugging symbols (OpenMP enabled)."
	@echo "💡 Run: gdb ./$(BIN)"
	@echo "💡 Then use: run | bt | print var"


clean:
	rm -f $(BIN) $(OBJ)
