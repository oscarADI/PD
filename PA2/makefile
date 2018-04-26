CC = g++
CFLAGS = -c -std=c++11
OPTFLAGS = -O3

##############
# Executable #
##############

Floorplan: main.o floorplan.o      
	@$(CC) $(OPTFLAGS) main.o floorplan.o -o Floorplan
################
# Object files #
################
main.o: src/main.cpp
	@$(CC) $(CFLAGS) $(OPTFLAGS) $< -o $@ 
floorplan.o: src/floorplan.cpp src/floorplan.h
	@$(CC) $(CFLAGS) $(OPTFLAGS) $< -o $@

#########
# Clean #
#########
clean:
	@rm -rf *.o Floorplan
