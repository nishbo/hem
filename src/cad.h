#ifndef CAD_H
#define CAD_H

/// Constants and definitions and libraries.

// Files import:
#define FILE_IMP_PARAMETERS         "./init/simulation.ini"
#define FILE_IMP_NEURONS            "./import/neurons.txt"
#define FILE_IMP_SYNAPSES           "./import/synapses.txt"

//Files export:
#define FILE_EXP_MAIN_OUTPUT        "./data/output.txt"
#define FILE_EXP_SPIKES             "./data/spikes.txt"
#define FILE_EXP_PARAMETERS         "./data/parameters.txt"
#define FILE_EXP_CONN_MATR          "./data/conn_matr.txt"
#define FILE_EXP_SYNAPSE            "./data/synapse.txt"
#define FILE_EXP_NEURON_PARAMS      "./data/neuron_parameters.txt"
#define FILE_EXP_WEIGHTS            "./data/weights.txt"
#define FILE_EXP_NEURONS            "./export/neurons.txt"
#define FILE_EXP_SYNAPSES           "./export/synapses.txt"
#define FILE_EXP_CURRENT            "./data/current.txt"

//0 for terminal, 1 for built-in, 2 for file
#define USE_DEFAULT_PARAMETERS      2
// Disables asking perm for overwriting,
#define ENABLE_TEST                 1
// below in mm

#define Malloc(n,t) (t*)std::malloc((n)*sizeof(t))

#define MIN(a,b) (((a)>(b))?(b):(a))
#define MAX(a,b) (((a)>(b))?(a):(b))

#endif // CAD_H
