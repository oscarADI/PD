CC = g++
CFLAGS = -c -std=c++11
OPTFLAGS = -O3

##############
# Executable #
##############

Floorplan: main.o       
	@$(CC) $(OPTFLAGS) main.o -o Floorplan
################
# Object files #
################
main.o: src/main.cpp
	@$(CC) $(CFLAGS) $(OPTFLAGS) $< -o $@ 

#########
# Clean #
#########
clean:
	@rm -rf *.o
