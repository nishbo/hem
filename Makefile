
CC            = gcc
CXX           = g++ #-m64
CFLAGS        = -g -Wall

DEL_FILE      = rm
CREATE_DIR    = mkdir -p

SOURCE_DIR    = src
OBJECT_DIR    = obj
SOURCES       = $(SOURCE_DIR)/main.cpp \
    $(SOURCE_DIR)/neuron.cpp \
    $(SOURCE_DIR)/synapse.cpp \
    $(SOURCE_DIR)/topology.cpp \
    $(SOURCE_DIR)/simulation.cpp \
    $(SOURCE_DIR)/vfdistributions.cpp \
    $(SOURCE_DIR)/vfdiscrete.cpp \
    $(SOURCE_DIR)/vfrandom.cpp \
    $(SOURCE_DIR)/inout.cpp
OBJECTS       = $(OBJECT_DIR)/main.o \
    $(OBJECT_DIR)/neuron.o \
    $(OBJECT_DIR)/synapse.o \
    $(OBJECT_DIR)/topology.o \
    $(OBJECT_DIR)/simulation.o \
    $(OBJECT_DIR)/vfdistributions.o \
    $(OBJECT_DIR)/vfdiscrete.o \
    $(OBJECT_DIR)/vfrandom.o \
    $(OBJECT_DIR)/inout.o
DESTDIR_TARGET = hem.exe

first: all

all: Makefile $(DESTDIR_TARGET)
	@echo ' '
	@echo 'Program succesfully built!'

$(DESTDIR_TARGET): dir $(OBJECTS)
	$(CXX) $(CFLAGS) -o $(DESTDIR_TARGET) $(OBJECTS)

dir: 
	$(CREATE_DIR) obj
	$(CREATE_DIR) data
	$(CREATE_DIR) export
	$(CREATE_DIR) import

clean:
	-$(DEL_FILE) $(OBJECTS)

$(OBJECT_DIR)/main.o: $(SOURCE_DIR)/main.cpp $(SOURCE_DIR)/cad.h \
    $(SOURCE_DIR)/neuron.h \
    $(SOURCE_DIR)/vfdistributions.h \
    $(SOURCE_DIR)/synapse.h \
    $(SOURCE_DIR)/vfdiscrete.h \
    $(SOURCE_DIR)/topology.h \
    $(SOURCE_DIR)/simulation.h \
    $(SOURCE_DIR)/vfrandom.h
	$(CXX) -c $(CFLAGS) -o $(OBJECT_DIR)/main.o $(SOURCE_DIR)/main.cpp

$(OBJECT_DIR)/neuron.o: $(SOURCE_DIR)/neuron.cpp $(SOURCE_DIR)/neuron.h \
    $(SOURCE_DIR)/cad.h \
    $(SOURCE_DIR)/vfdistributions.h
	$(CXX) -c $(CFLAGS) -o $(OBJECT_DIR)/neuron.o $(SOURCE_DIR)/neuron.cpp

$(OBJECT_DIR)/synapse.o: $(SOURCE_DIR)/synapse.cpp $(SOURCE_DIR)/synapse.h \
    $(SOURCE_DIR)/cad.h \
    $(SOURCE_DIR)/vfdistributions.h \
    $(SOURCE_DIR)/vfdiscrete.h
	$(CXX) -c $(CFLAGS) -o $(OBJECT_DIR)/synapse.o $(SOURCE_DIR)/synapse.cpp

$(OBJECT_DIR)/topology.o: $(SOURCE_DIR)/topology.cpp $(SOURCE_DIR)/topology.h \
    $(SOURCE_DIR)/cad.h \
    $(SOURCE_DIR)/vfdistributions.h \
    $(SOURCE_DIR)/vfdiscrete.h
	$(CXX) -c $(CFLAGS) -o $(OBJECT_DIR)/topology.o $(SOURCE_DIR)/topology.cpp

$(OBJECT_DIR)/simulation.o: $(SOURCE_DIR)/simulation.cpp $(SOURCE_DIR)/simulation.h \
    $(SOURCE_DIR)/cad.h \
    $(SOURCE_DIR)/neuron.h \
    $(SOURCE_DIR)/vfdistributions.h \
    $(SOURCE_DIR)/synapse.h \
    $(SOURCE_DIR)/vfdiscrete.h \
    $(SOURCE_DIR)/topology.h \
    $(SOURCE_DIR)/vfrandom.h
	$(CXX) -c $(CFLAGS) -o $(OBJECT_DIR)/simulation.o $(SOURCE_DIR)/simulation.cpp

$(OBJECT_DIR)/vfdistributions.o: $(SOURCE_DIR)/vfdistributions.cpp $(SOURCE_DIR)/vfdistributions.h
	$(CXX) -c $(CFLAGS) -o $(OBJECT_DIR)/vfdistributions.o $(SOURCE_DIR)/vfdistributions.cpp

$(OBJECT_DIR)/vfdiscrete.o: $(SOURCE_DIR)/vfdiscrete.cpp $(SOURCE_DIR)/vfdiscrete.h \
    $(SOURCE_DIR)/cad.h
	$(CXX) -c $(CFLAGS) -o $(OBJECT_DIR)/vfdiscrete.o $(SOURCE_DIR)/vfdiscrete.cpp

$(OBJECT_DIR)/vfrandom.o: $(SOURCE_DIR)/vfrandom.cpp $(SOURCE_DIR)/vfrandom.h
	$(CXX) -c $(CFLAGS) -o $(OBJECT_DIR)/vfrandom.o $(SOURCE_DIR)/vfrandom.cpp

$(OBJECT_DIR)/inout.o: $(SOURCE_DIR)/inout.cpp $(SOURCE_DIR)/simulation.h \
    $(SOURCE_DIR)/cad.h \
    $(SOURCE_DIR)/neuron.h \
    $(SOURCE_DIR)/vfdistributions.h \
    $(SOURCE_DIR)/synapse.h \
    $(SOURCE_DIR)/vfdiscrete.h \
    $(SOURCE_DIR)/topology.h \
    $(SOURCE_DIR)/vfrandom.h
	$(CXX) -c $(CFLAGS) -o $(OBJECT_DIR)/inout.o $(SOURCE_DIR)/inout.cpp
