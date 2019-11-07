all: tools.so

tools.so: toolsf.f Makefile
	f2py3  -c -m tools only: committor evecs npcommittor : toolsf.f -llapack


