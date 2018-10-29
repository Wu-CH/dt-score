CPLUS = g++
CFLAGS = -O3

all: DT-Score coverage

DT-Score: DT-Score/dt-score.o
	${CPLUS} $(CFLAGS) DT-Score/dt-score.o -o DT-Score/dt-score
dt-score.o: DT-Score/dt-score.cpp
	${CPLUS} $(CFLAGS) -c DT-Score/dt-score.cpp

coverage: coverage/coverage.o
	${CPLUS} $(CFLAGS) coverage/coverage.o -o coverage/coverage
coverage.o: coverage/coverage.cpp
	${CPLUS} $(CFLAGS) -c coverage/coverage.cpp

clean:
	-rm -f DT-Score/dt-score.o DT-Score/dt-score
	-rm -f coverage/coverage.o coverage/coverage
