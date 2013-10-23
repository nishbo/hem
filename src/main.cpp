#include "cad.h"
#include "neuron.h"
#include "synapse.h"
#include "topology.h"
#include "simulation.h"
#include "vfdistributions.h"
#include "vfdiscrete.h"

#include <iomanip>
#include <locale>
#include <sstream>
#include <string>

using namespace std;

int main(){
    srand(time(NULL));
    cout.precision(5);

    cout<<"\tWelcome to HEM simulator.\n";

    SimulationSingleton *storage = SimulationSingleton::instance();
    storage->createNetwork();

    cout<<"\nSaving network parameters to file...\r";
    storage->outputParametersInFile();
    storage->outputConnectivityMatrixInFile();
    cout<<"Network parameters saved to file.   ";

    storage->setOutputFile();

    storage->exportNeurons();
    storage->exportSynapses();
    storage->exportStimulation();

    storage->length_of_simulation += storage->dt;

    //************
    int stamp = time(NULL);
    srand(stamp);
    std::cout<<"\nSEED: "<<stamp<<std::endl;
    std::cout<<"\nFirst rand: "<<rand()<<std::endl;
    storage->createSynapseStimulation();
    std::cout<<"\nSecond rand: "<<rand()<<std::endl;

    cout<<"\n\tStarting simulation.\n";

    for(storage->time_now = 0.0; \
        storage->time_now < storage->length_of_simulation ;\
        storage->time_now += storage->dt){

        storage->evolveAllNeurons();
        storage->sendNeuralNoise();
        storage->sendSynapseNoise();
        storage->evolveAllSynapses();
        storage->outputChangingData();

        cout<<"Finished "<<fixed<<(storage->time_now)/ \
                                  (storage->length_of_simulation)<<"\r";
    }
    std::cout<<"\nThird rand: "<<rand()<<std::endl;

    storage->closeOutputFile();
    storage->saveWeights("data/weight_target.txt");
    /// LOL**********************************************************************
    cout<<"\n\n\tStarting simulation. AGAIN\n";
    //************
    storage->resetNetwork();
    srand(stamp);
    std::cout<<"\nSEED: "<<stamp<<std::endl;
    std::cout<<"\nFirst rand: "<<rand()<<std::endl;
    storage->createSynapseStimulation();
    std::cout<<"\nSecond rand: "<<rand()<<std::endl;
    //************
    storage->saveWeights("data/weight_2sim.txt");

    for(storage->time_now = 0.0; \
        storage->time_now < storage->length_of_simulation ;\
        storage->time_now += storage->dt){

        storage->evolveAllNeurons();
        storage->sendNeuralNoise();
        storage->sendSynapseNoise();
        storage->evolveAllSynapses();
        storage->outputChangingData();

        if(vf_discrete::inBetween(storage->time_now, 100, storage->dt)){
            storage->saveWeights("data/we/weight_"+\
                static_cast<ostringstream*>( &(ostringstream() << storage->time_now) )->str() \
                +".txt");
        }

        cout<<"Finished "<<fixed<<(storage->time_now)/ \
                                  (storage->length_of_simulation)<<"\r";
    }
    std::cout<<"\nThird rand: "<<rand()<<std::endl;

    storage->saveWeights("data/weight_final.txt");
    storage->closeOutputFile();

    cout<<"Finished 1.00000"<<endl<<"Press RETURN to exit..."<<endl;
    // storage->test();
    cin.ignore();
    return 0;
}

