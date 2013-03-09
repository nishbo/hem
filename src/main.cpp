#include "cad.h"
#include "neuron.h"
#include "synapse.h"
#include "topology.h"
#include "simulation.h"
#include "vfdistributions.h"
#include "vfdiscrete.h"

using namespace std;

int main(){
    srand(time(NULL));
    cout<<"Welcome to HEM simulator.\n\n";

    SimulationSingleton *storage = SimulationSingleton::instance();
    storage->loadParameters();
    storage->errorReport();

    storage->createNetwork();
    storage->errorReport();

    storage->outputParametersInFile();
    storage->outputConnectivityMatrixInFile();
    storage->setOutputFile();
    storage->errorReport();

    storage->exportNeurons();
    storage->exportSynapses();

    storage->length_of_simulation += storage->dt;
    cout.precision(5);
    cout<<endl;
    for(storage->time_now = 0.0; \
        storage->time_now < storage->length_of_simulation ;\
        storage->time_now += storage->dt){

        storage->sendNeuralNoise();
        storage->sendSynapseNoise();
        storage->evolveAllNeurons();
        storage->evolveAllSynapses();
        storage->outputChangingData();

        cout<<"Finished "<<fixed<<(storage->time_now)/ \
                                  (storage->length_of_simulation)<<"\r";

        storage->errorReport();
    }

    storage->closeOutputFile();

    cout<<"Finished 1.00000"<<endl<<"Press any key..."<<endl;
    getch();
    return 0;
}

