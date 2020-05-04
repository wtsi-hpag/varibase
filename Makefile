CC = cc
CFLAGS = -O4

LIB = -lm

OBJS=	\
		varibase.o \
                fast.o \

SRCS= $(OBJS:.o=.c)

INCS=   \
                fasta.h \

PROGRAM = varibase

$(PROGRAM): $(OBJS) $(INCS)
	$(CC) $(CFLAGS) -o $(PROGRAM) $(OBJS) $(LIB)

clean:
	rm -f $(PROGRAM) *.o
