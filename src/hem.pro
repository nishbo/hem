TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp \
    neuron.cpp \
    synapse.cpp \
    topology.cpp \
    simulation.cpp \
    vfdistributions.cpp \
    vfdiscrete.cpp \
    vfrandom.cpp \
    inout.cpp

HEADERS += \
    neuron.h \
    synapse.h \
    topology.h \
    var_functions.h \
    cad.h \
    simulation.h \
    vfdistributions.h \
    vfdiscrete.h \
    vfrandom.h

OTHER_FILES += \
    Makefile.a

