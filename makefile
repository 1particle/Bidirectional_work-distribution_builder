 CC = g++
 CFLAGS  = -g -Wall
 main = src/Bidirectional_work-distribution_builder
 bin = bin_Linux_x86_64/Bidirectional_work-distribution_builder__Linux_x86_64
 util = src/distributionsCalculation



$(bin) :	Bidirectional_work-distribution_builder.o distributionsCalculation.o
	$(CC) $(CFLAGS)	-o	$(bin)	Bidirectional_work-distribution_builder.o	distributionsCalculation.o

Bidirectional_work-distribution_builder.o :  $(main).cpp
	$(CC) $(CFLAGS) -c $(main).cpp

distributionsCalculation.o : $(util).cpp  $(util).h 
	$(CC) $(CFLAGS) -c $(util).cpp

clean: 
	$(RM) count *.o *~
