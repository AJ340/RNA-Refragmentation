##
# @file      Makefile
# @author    Andres Quinones
# @Copyright None
# @brief     Builds RNA-Refragmentation Solution


GCC      = g++ -std=c++11

SRCS       = main.cpp
OBJS       = $(patsubst %.c,%.o,$(SRCS))
EXECS      = RNA-Refragmentation-AQ 

CFLAGS     = 


%.o : %.c
	$(GCC) -c $< -o $@ $(CFLAGS) 

all:
	make $(EXECS)

RNA-Refragmentation-AQ: $(OBJS)
	$(GCC) -o $@ $(OBJS) $(CFLAGS) 

clean:
	-rm -f $(EXECS) 

cleanall:
	-rm -f $(EXECS) 

