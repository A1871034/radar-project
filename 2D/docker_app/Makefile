main: *.cpp Dockerfile
	g++ -fopenmp -foffload=nvptx-none -fno-stack-protector -fcf-protection=none -g gpu_test.cpp -o gpu_test.exe
