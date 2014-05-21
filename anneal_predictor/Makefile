CC=gcc

FRAMEWORKS:= -framework Foundation
LIBRARIES:= -lobjc -lgnustep-base

SOURCE=main.m

CFLAGS=-Wall -Werror -g -v $(SOURCE)
LDFLAGS=$(LIBRARIES) $(FRAMEWORKS)
OUT=-o main

all:
	gcc -lobjc -o main main.m -std=c11 -fconstant-string-class=NSConstantString -lgnustep-base -g 
