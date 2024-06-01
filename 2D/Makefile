main: *.cpp *.h
	g++ -fopenmp -foffload=nvptx-none -fno-stack-protector -fcf-protection=none -foffload=-lm -fno-fast-math -fno-associative-math -lm -ldl -O3 -g *.cpp -o main.exe