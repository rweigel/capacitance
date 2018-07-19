all: pushcode.c
	gcc -O3 pushcode.c -lm -o pushcode

test:
	- mkdir tmp
	printf '%s\n%s\n%s' "-1" "0" "1" > tmp/pushcode_xo.txt
	./pushcode 3 10

