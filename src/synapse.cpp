#include "synapse.h"
#include <random>
#include <ctime>

double* Synapse::innerDataArr = NULL;

void Synapse::incSpike(double timen){
    last_spiked = timen;
}
void Synapse::incAfterSpike(double timen){
    last_spiked_post = timen;
}
void Synapse::setDeliveries(double dt){
    number_of_possible_transfers = delay / dt;
    delivery = new double[number_of_possible_transfers];
    for(int i=0; i<number_of_possible_transfers; i++)
        delivery[i] = 0;
}
double Synapse::moveDeliveries(){
    buf = delivery[0];
    for(int i=0; i<number_of_possible_transfers-1; i++)
        delivery[i] = delivery[i+1];
    delivery[number_of_possible_transfers-1] = out_current;
    return buf;
}
double* Synapse::getInnerData(){
    return innerDataArr;
}

std::string Synapse::getName(){
    return "Name_not_set";
}

int Synapse::importData(double *arr){
    return 0;
}

double* Synapse::exportData(){
    double* arr = new double[1];
    arr[0] = 0;
    return arr;
}

int Synapse::initSynapses(){
    using namespace vf_file;
    FILE* fid = fopen (getFilenameFromIni(DATAFILES, FILE_EXP_INIT_SYNAPSES).c_str(), "w");
    if(!fid){
        exit(15);
    }
    fprintf(fid, "No data specified.\n");
    fclose(fid);
    return 1;
}

int Synapse::numEssentialVariables(){
    return 0;
}

int Synapse::to(){
    return postsynaptic;
}

int Synapse::from(){
    return presynaptic;
}

double Synapse::test(){
    return -1;
}

/// Static synapse

std::string SynapseStatic::synapsetype = "Static_synapse";
double SynapseStatic::init_weight = 0.5;

int SynapseStatic::initSynapsesLocal(){
    using namespace vf_file;
    std::string buf30 = loadFileToString("./init/SynapseStatic.ini");
    init_weight = getParameterIni("Weight", buf30);

    FILE* fid = fopen (getFilenameFromIni(DATAFILES, FILE_EXP_INIT_SYNAPSES).c_str(), "w");
    if(!fid){
        exit(15);
    }
    fprintf(fid, "weight = %f;\n", init_weight);
    fclose(fid);

    return 0;
}

double* SynapseStatic::exportData(){
    double* arr = new double[numEssentialVariables()+1];
    arr[0] = numEssentialVariables()+1;
    arr[1] = presynaptic;
    arr[2] = postsynaptic;
    arr[3] = weight;
    return arr;
}

int SynapseStatic::importData(double *arr){
    presynaptic = arr[0];
    postsynaptic = arr[1];
    weight = arr[3];

    last_spiked = -100;
    last_spiked_post = -100;
    return 0;
}

int SynapseStatic::initSynapses(){
    return initSynapsesLocal();
}

std::string SynapseStatic::getName(){
    return synapsetype;
}

double SynapseStatic::evolve(double dt, double time, double Vpre, double Vpost){
    out_current = weight;

    return moveDeliveries();
}

void SynapseStatic::setData (int pre, int pos, int preex, int posex, double dt){
    // sets data like pre- and postsyn neuron, etc.
    presynaptic = pre;
    postsynaptic = pos;
    weight = init_weight;
    last_spiked = -100;
    last_spiked_post = -100;
}

int SynapseStatic::numEssentialVariables(){
    return 3;
}

/// Tsodyks-Markram model rk4

std::string SynapseTsodyksMarkramRK::synapsetype = \
        "Tsodyks-Markram_synapse_(rk4)";

std::string SynapseTsodyksMarkramRK::getName(){
    return synapsetype;
}

double SynapseTsodyksMarkramRK::init_tau_one = 3;
double SynapseTsodyksMarkramRK::init_x = 0.98;
double SynapseTsodyksMarkramRK::init_y = 0.01;
double SynapseTsodyksMarkramRK::init_z = 0.01;
double SynapseTsodyksMarkramRK::init_Aee = 38;
double SynapseTsodyksMarkramRK::init_Uee = 0.5;
double SynapseTsodyksMarkramRK::init_tau_recee = 800;
double SynapseTsodyksMarkramRK::init_tau_facilee = 0;
double SynapseTsodyksMarkramRK::init_Aei = 54;
double SynapseTsodyksMarkramRK::init_Uei = 0.5;
double SynapseTsodyksMarkramRK::init_tau_recei = 800;
double SynapseTsodyksMarkramRK::init_tau_facilei = 0;
double SynapseTsodyksMarkramRK::init_Aie = -72;
double SynapseTsodyksMarkramRK::init_Uie = 0.04;
double SynapseTsodyksMarkramRK::init_tau_recie = 100;
double SynapseTsodyksMarkramRK::init_tau_facilie = 1000;
double SynapseTsodyksMarkramRK::init_Aii = -72;
double SynapseTsodyksMarkramRK::init_Uii = 0.04;
double SynapseTsodyksMarkramRK::init_tau_recii = 100;
double SynapseTsodyksMarkramRK::init_tau_facilii = 1000;
int SynapseTsodyksMarkramRK::init_distribute_params = 0;
int SynapseTsodyksMarkramRK::synamount = 0;
double SynapseTsodyksMarkramRK::xav = 0;
double SynapseTsodyksMarkramRK::yav = 0;
double SynapseTsodyksMarkramRK::zav = 0;
double SynapseTsodyksMarkramRK::uav = 0;
double* SynapseTsodyksMarkramRK::innerDataArr = NULL;

void SynapseTsodyksMarkramRK::setData(int pre, int pos, \
                                    int preex, int posex, double dt){
    tau_one = init_tau_one;

    presynaptic = pre;
    postsynaptic = pos;
    weight = 1;
    last_spiked = -100;
    last_spiked_post = -100;

    if(preex && posex){
        //exc exc
        A = init_Aee;
        U = init_Uee;
        tau_rec = init_tau_recee;
        tau_facil = init_tau_facilee;
        exc = 1;
    } else if(preex){
        //exc inh
        A = init_Aei;
        U = init_Uei;
        tau_rec = init_tau_recei;
        tau_facil = init_tau_facilei;
        exc = 1;
    } else if(posex){
        // inh exc
        A = init_Aie;
        U = init_Uie;
        tau_rec = init_tau_recie;
        tau_facil = init_tau_facilie;
        exc = 0;
    } else {
        // inh inh
        A = init_Aii;
        U = init_Uii;
        tau_rec = init_tau_recii;
        tau_facil = init_tau_facilii;
        exc = 0;
    }
    if(init_distribute_params){
        if(A>0)
            A = vf_distributions::normal(A, A/2, 0, 4*A);
        else
            A = vf_distributions::normal(A, -A/2, 4*A, 0);
        U = vf_distributions::normal(U, U/2, 0, MIN(1, 4*U));
        tau_rec = vf_distributions::normal(tau_rec, tau_rec/2, dt, tau_rec*4);
        tau_facil = vf_distributions::normal(tau_facil, tau_facil/2, \
                                            dt, tau_facil*4);
    }

    x = init_x;
    y = init_y;
    z = init_z;
    u = U;

    synamount++;
}

double SynapseTsodyksMarkramRK::Xr(double x1, double y1, double z1, double u1, \
                                  double t1, double dt){
    return z1/tau_rec - vf_discrete::diracDelta(t1 - last_spiked, dt) * u1 * x1;
}

double SynapseTsodyksMarkramRK::Yr(double x1, double y1, double z1, double u1, \
                                  double t1, double dt){
    return -y1/tau_one + vf_discrete::diracDelta(t1 - last_spiked, dt) * u1 * x1;
}

double SynapseTsodyksMarkramRK::Zr(double x1, double y1, double z1, double u1, \
                                  double t1, double dt){
    return -z1/tau_rec + y1/tau_one;
}

double SynapseTsodyksMarkramRK::Ur(double x1, double y1, double z1, double u1, \
                                  double t1, double dt){
    if(exc)
        return U;
    else
        return -u1/tau_facil + vf_discrete::diracDelta(t1 - last_spiked, dt) \
                * U * (1 - u1);
}

double SynapseTsodyksMarkramRK::evolve(double dt, double time, \
                                       double Vpre, double Vpost){

    k_1_x = dt * Xr(x, y, z, u, time, dt);
    k_1_y = dt * Yr(x, y, z, u, time, dt);
    k_1_z = dt * Zr(x, y, z, u, time, dt);
    if(exc)
        k_1_u = dt * U;
    else
        k_1_u = dt * Ur(x, y, z, u, time, dt);

    k_2_x = dt * Xr(x + k_1_x/2.0, y + k_1_y/2.0, z + k_1_z/2.0, \
                    u + k_1_u/2.0, time + dt/2.0, dt);
    k_2_y = dt * Yr(x + k_1_x/2.0, y + k_1_y/2.0, z + k_1_z/2.0, \
                    u + k_1_u/2.0, time + dt/2.0, dt);
    k_2_z = dt * Zr(x + k_1_x/2.0, y + k_1_y/2.0, z + k_1_z/2.0, \
                    u + k_1_u/2.0, time + dt/2.0, dt);
    if(exc)
        k_2_u = dt * k_1_u;
    else
        k_2_u = dt * Ur(x + k_1_x/2.0, y + k_1_y/2.0, z + k_1_z/2.0, \
                    u + k_1_u/2.0, time + dt/2.0, dt);

    k_3_x = dt * Xr(x + k_2_x/2.0, y + k_2_y/2.0, z + k_2_z/2.0, \
                    u + k_2_u/2.0, time + dt/2.0, dt);
    k_3_y = dt * Yr(x + k_2_x/2.0, y + k_2_y/2.0, z + k_2_z/2.0, \
                    u + k_2_u/2.0, time + dt/2.0, dt);
    k_3_z = dt * Zr(x + k_2_x/2.0, y + k_2_y/2.0, z + k_2_z/2.0, \
                    u + k_2_u/2.0, time + dt/2.0, dt);
    if(exc)
        k_3_u = dt * k_2_u;
    else
        k_3_u = dt * Ur(x + k_2_x/2.0, y + k_2_y/2.0, z + k_2_z/2.0, \
                    u + k_2_u/2.0, time + dt/2.0, dt);

    k_4_x = dt * Xr(x + k_3_x, y + k_3_y, z + k_3_z, u + k_3_u, time + dt, dt);
    k_4_y = dt * Yr(x + k_3_x, y + k_3_y, z + k_3_z, u + k_3_u, time + dt, dt);
    k_4_z = dt * Zr(x + k_3_x, y + k_3_y, z + k_3_z, u + k_3_u, time + dt, dt);
    if(exc)
        k_4_u = dt * U;
    else
        k_4_u = dt * Ur(x + k_3_x, y + k_3_y, z + k_3_z, u + k_3_u, time + dt, dt);

    x += (k_1_x + 2.0*k_2_x + 2.0*k_3_x + k_4_x)/6.0;
    y += (k_1_y + 2.0*k_2_y + 2.0*k_3_y + k_4_y)/6.0;
    z += (k_1_z + 2.0*k_2_z + 2.0*k_3_z + k_4_z)/6.0;
    if(!exc)
        u += (k_1_u + 2.0*k_2_u + 2.0*k_3_u + k_4_u)/6.0;

    //mistakes:
    if(x >= 1) x = 0.9999999; if(x<=0) x = 0.0000001;
    if(y >= 1) y = 0.9999999; if(y<=0) y = 0.0000001;
    if(z >= 1) z = 0.9999999; if(z<=0) z = 0.0000001;
    if(u >= 1) u = 0.9999999; if(u<=0) u = 0.0000001;

    out_current = A * y;

    xav += x;
    yav += y;
    zav += z;
    uav += u;

    return moveDeliveries();
}

double* SynapseTsodyksMarkramRK::exportData(){
    double* arr = new double[numEssentialVariables()+1];
    arr[0] = numEssentialVariables()+1;
    arr[1] = presynaptic;
    arr[2] = postsynaptic;
    arr[3] = weight;
    arr[4] = delay;
    arr[5] = x;
    arr[6] = y;
    arr[7] = z;
    arr[8] = u;
    arr[9] = tau_one;
    arr[10] = tau_rec;
    arr[11] = tau_facil;
    arr[12] = U;
    arr[13] = A;
    arr[14] = exc;

    return arr;
}

int SynapseTsodyksMarkramRK::importData(double *arr){
    presynaptic = arr[0];
    postsynaptic = arr[1];
    weight = arr[2];
    delay = arr[3];
    x = arr[4];
    y = arr[5];
    z = arr[6];
    u = arr[7];
    tau_one = arr[8];
    tau_rec = arr[9];
    tau_facil = arr[10];
    U = arr[11];
    A = arr[12];
    exc = arr[13];

    synamount++;
    xav += x;
    yav += y;
    zav += z;
    uav += u;
    last_spiked = -100;
    last_spiked_post = -100;
    return 0;
}

int SynapseTsodyksMarkramRK::numEssentialVariables(){
    return 14;
}

int SynapseTsodyksMarkramRK::initSynapsesLocal(){
    using namespace vf_file;
    std::string buf30 = loadFileToString("./init/SynapseTsodyksMarkram.ini");

    init_tau_one = getParameterIni("tau_one", buf30);
    init_x = getParameterIni("INIT_x", buf30);
    init_y = getParameterIni("INIT_y", buf30);
    init_z = getParameterIni("INIT_z", buf30);

    init_Aee = getParameterIni("Aee", buf30);
    init_Uee = getParameterIni("Uee", buf30);
    init_tau_recee = getParameterIni("tau_recee", buf30);
    init_tau_facilee = getParameterIni("tau_facilee", buf30);

    init_Aei = getParameterIni("Aei", buf30);
    init_Uei = getParameterIni("Uei", buf30);
    init_tau_recei = getParameterIni("tau_recei", buf30);
    init_tau_facilei = getParameterIni("tau_facilei", buf30);

    init_Aie = getParameterIni("Aie", buf30);
    init_Uie = getParameterIni("Uie", buf30);
    init_tau_recie = getParameterIni("tau_recie", buf30);
    init_tau_facilie = getParameterIni("tau_facilie", buf30);

    init_Aii = getParameterIni("Aii", buf30);
    init_Uii = getParameterIni("Uii", buf30);
    init_tau_recii = getParameterIni("tau_recii", buf30);
    init_tau_facilii = getParameterIni("tau_facilii", buf30);

    init_distribute_params = getParameterIni("DISTRIBUTE_PARAMETERS", buf30);

    FILE* fid = fopen (getFilenameFromIni(DATAFILES, FILE_EXP_INIT_SYNAPSES).c_str(), "w");
    if(!fid){
        exit(15);
    }
    fprintf(fid, "tau_one = %.4f;\n", init_tau_one);
    fprintf(fid, "x = %.4f;\n", init_x);
    fprintf(fid, "y = %.4f;\n", init_y);
    fprintf(fid, "z = %.4f;\n\n", init_z);

    fprintf(fid, "Aee = %.4f;\n", init_Aee);
    fprintf(fid, "Uee = %.4f;\n", init_Uee);
    fprintf(fid, "tau_recee = %.4f;\n", init_tau_recee);
    fprintf(fid, "tau_facilee = %.4f;\n\n", init_tau_facilee);

    fprintf(fid, "Aei = %.4f;\n", init_Aei);
    fprintf(fid, "Uei = %.4f;\n", init_Uei);
    fprintf(fid, "tau_recei = %.4f;\n", init_tau_recei);
    fprintf(fid, "tau_facilei = %.4f;\n\n", init_tau_facilei);

    fprintf(fid, "Aie = %.4f;\n", init_Aie);
    fprintf(fid, "Uie = %.4f;\n", init_Uie);
    fprintf(fid, "tau_recie = %.4f;\n", init_tau_recie);
    fprintf(fid, "tau_facilie = %.4f;\n\n", init_tau_facilie);

    fprintf(fid, "Aii = %.4f;\n", init_Aii);
    fprintf(fid, "Uii = %.4f;\n", init_Uii);
    fprintf(fid, "tau_recii = %.4f;\n", init_tau_recii);
    fprintf(fid, "tau_facilii = %.4f;\n\n", init_tau_facilii);

    fprintf(fid, "DISTRIBUTE_PARAMETERS = %d;\n", init_distribute_params);
    fclose(fid);

    return 0;
}

int SynapseTsodyksMarkramRK::initSynapses(){
    synamount = 0;
    xav = 0;
    yav = 0;
    zav = 0;
    uav = 0;
    return initSynapsesLocal();
}

double* SynapseTsodyksMarkramRK::getInnerData(){
    innerDataArr = new double [5];
    innerDataArr[0] = 5;
    innerDataArr[1] = xav / synamount;
    innerDataArr[2] = yav / synamount;
    innerDataArr[3] = zav / synamount;
    innerDataArr[4] = uav / synamount;
    xav = 0;
    yav = 0;
    zav = 0;
    uav = 0;
    return innerDataArr;
}

double SynapseTsodyksMarkramRK::test(){
    return y;
}

/// Tsodyks-Markram rk nest style

std::string SynapseTsodyksMarkramRKNest::synapsetype = \
        "Tsodyks-Markram_synapse_(rk4_nest)";

std::string SynapseTsodyksMarkramRKNest::getName(){
    return synapsetype;
}

double SynapseTsodyksMarkramRKNest::init_tau_one = 3;
double SynapseTsodyksMarkramRKNest::init_x = 0.98;
double SynapseTsodyksMarkramRKNest::init_y = 0.01;
double SynapseTsodyksMarkramRKNest::init_z = 0.01;
double SynapseTsodyksMarkramRKNest::init_Aee = 38;
double SynapseTsodyksMarkramRKNest::init_Uee = 0.5;
double SynapseTsodyksMarkramRKNest::init_tau_recee = 800;
double SynapseTsodyksMarkramRKNest::init_tau_facilee = 0;
double SynapseTsodyksMarkramRKNest::init_Aei = 54;
double SynapseTsodyksMarkramRKNest::init_Uei = 0.5;
double SynapseTsodyksMarkramRKNest::init_tau_recei = 800;
double SynapseTsodyksMarkramRKNest::init_tau_facilei = 0;
double SynapseTsodyksMarkramRKNest::init_Aie = -72;
double SynapseTsodyksMarkramRKNest::init_Uie = 0.04;
double SynapseTsodyksMarkramRKNest::init_tau_recie = 100;
double SynapseTsodyksMarkramRKNest::init_tau_facilie = 1000;
double SynapseTsodyksMarkramRKNest::init_Aii = -72;
double SynapseTsodyksMarkramRKNest::init_Uii = 0.04;
double SynapseTsodyksMarkramRKNest::init_tau_recii = 100;
double SynapseTsodyksMarkramRKNest::init_tau_facilii = 1000;
double SynapseTsodyksMarkramRKNest::init_tau_ge = 3;
double SynapseTsodyksMarkramRKNest::init_tau_gi = 7;
double SynapseTsodyksMarkramRKNest::init_g = 0;
double SynapseTsodyksMarkramRKNest::init_Ee = 0;
double SynapseTsodyksMarkramRKNest::init_Ei = -80;

void SynapseTsodyksMarkramRKNest::setData(int pre, int pos, \
                                    int preex, int posex, double dt){
    tau_one = init_tau_one;

    presynaptic = pre;
    postsynaptic = pos;
    weight = 1;
    last_spiked = -100;
    lsp2 = 0;
    last_spiked_post = -100;

    if(preex && posex){
        //exc exc
        A = init_Aee;
        U = init_Uee;
        tau_rec = init_tau_recee;
        tau_facil = init_tau_facilee;
        exc = 1;
    } else if(preex){
        //exc inh
        A = init_Aei;
        U = init_Uei;
        tau_rec = init_tau_recei;
        tau_facil = init_tau_facilei;
        exc = 1;
    } else if(posex){
        // inh exc
        A = init_Aie;
        U = init_Uie;
        tau_rec = init_tau_recie;
        tau_facil = init_tau_facilie;
        exc = 0;
    } else {
        // inh inh
        A = init_Aii;
        U = init_Uii;
        tau_rec = init_tau_recii;
        tau_facil = init_tau_facilii;
        exc = 0;
    }
    if(A>0)
        A = vf_distributions::normal(A, A/2, 0, 4*A);
    else
        A = vf_distributions::normal(A, -A/2, 4*A, 0);
    U = vf_distributions::normal(U, U/2, 0, MIN(1, 4*U));
    tau_rec = vf_distributions::normal(tau_rec, tau_rec/2, dt, tau_rec*4);
    tau_facil = vf_distributions::normal(tau_facil, tau_facil/2, \
                                        dt, tau_facil*4);

    xo = init_x;
    yo = init_y;
    zo = init_z;
    uo = U;

    if(exc){
        E = init_Ee;
        tau_g = init_tau_ge;
    } else {
        E = init_Ei;
        tau_g = init_tau_gi;
    }
    g = init_g;
}

double SynapseTsodyksMarkramRKNest::evolve(double dt, double time, \
                                       double Vpre, double Vpost){

    if(vf_discrete::diracDelta(time - last_spiked, dt)){
        h = last_spiked - lsp2;
        lsp2 = last_spiked;

        zo = 1 - xo - yo;
        if(exc)
            u1 = 0;
        else
            u1 = uo * exp(- h / tau_facil);
        x1 = xo + ((exp(-h/tau_rec)-1)*tau_rec - (exp(-h/tau_one)-1)*tau_one) *\
                yo / (tau_one - tau_rec) + (1 - exp(-h/tau_rec))*zo;
        y1 = yo * exp(-h/tau_one);
        u2 = u1 + U * (1-u1);

        g += A * u2 * x1;

        xo = x1 - u2*x1;
        yo = y1 + u2*x1;
    }

    g -= dt * (g / tau_g);

    out_current = g;

    return moveDeliveries();
}

double* SynapseTsodyksMarkramRKNest::exportData(){
    double* arr = new double[numEssentialVariables()+1];
    arr[0] = numEssentialVariables()+1;
    arr[1] = presynaptic;
    arr[2] = postsynaptic;
    arr[3] = weight;
    arr[4] = delay;
    arr[5] = xo;
    arr[6] = yo;
    arr[7] = zo;
    arr[8] = uo;
    arr[9] = tau_one;
    arr[10] = tau_rec;
    arr[11] = tau_facil;
    arr[12] = U;
    arr[13] = A;
    arr[14] = exc;
    arr[15] = g;
    arr[16] = tau_g;
    arr[17] = E;
    return arr;
}

int SynapseTsodyksMarkramRKNest::importData(double *arr){
    presynaptic = arr[0];
    postsynaptic = arr[1];
    weight = arr[2];
    delay = arr[3];
    xo = arr[4];
    yo = arr[5];
    zo = arr[6];
    uo = arr[7];
    tau_one = arr[8];
    tau_rec = arr[9];
    tau_facil = arr[10];
    U = arr[11];
    A = arr[12];
    exc = arr[13];
    g = arr[14];
    tau_g = arr[15];
    E = arr[16];

    if(tau_rec < 1e-5)
        tau_rec = 1e-5;
    if(tau_facil < 1e-5)
        tau_facil = 1e-5;

    last_spiked = -100;
    lsp2 = -100;
    last_spiked_post = -100;

    return 0;
}

int SynapseTsodyksMarkramRKNest::numEssentialVariables(){
    return 17;
}

int SynapseTsodyksMarkramRKNest::initSynapsesLocal(){
    using namespace vf_file;
    std::string buf30 = loadFileToString("./init/SynapseTsodyksMarkramRKNest.ini");

    init_tau_one = getParameterIni("tau_one", buf30);
    init_x = getParameterIni("x", buf30);
    init_y = getParameterIni("y", buf30);
    init_z = getParameterIni("z", buf30);

    init_Aee = getParameterIni("Aee", buf30);
    init_Uee = getParameterIni("Uee", buf30);
    init_tau_recee = getParameterIni("tau_recee", buf30);
    init_tau_facilee = getParameterIni("tau_facilee", buf30);

    init_Aei = getParameterIni("Aei", buf30);
    init_Uei = getParameterIni("Uei", buf30);
    init_tau_recei = getParameterIni("tau_recei", buf30);
    init_tau_facilei = getParameterIni("tau_facilei", buf30);

    init_Aie = getParameterIni("Aie", buf30);
    init_Uie = getParameterIni("Uie", buf30);
    init_tau_recie = getParameterIni("tau_recie", buf30);
    init_tau_facilie = getParameterIni("tau_facilie", buf30);

    init_Aii = getParameterIni("Aii", buf30);
    init_Uii = getParameterIni("Uii", buf30);
    init_tau_recii = getParameterIni("tau_recii", buf30);
    init_tau_facilii = getParameterIni("tau_facilii", buf30);

    init_g = getParameterIni("g", buf30);
    init_Ee = getParameterIni("Ee", buf30);
    init_Ei = getParameterIni("Ei", buf30);
    init_tau_ge = getParameterIni("tau_ge", buf30);
    init_tau_gi = getParameterIni("tau_gi", buf30);

    FILE* fid = fopen (getFilenameFromIni(DATAFILES, FILE_EXP_INIT_SYNAPSES).c_str(), "w");
    if(!fid){
        exit(15);
    }
    fprintf(fid, "tau_one = %.4f;\n", init_tau_one);
    fprintf(fid, "x = %.4f;\n", init_x);
    fprintf(fid, "y = %.4f;\n", init_y);
    fprintf(fid, "z = %.4f;\n\n", init_z);

    fprintf(fid, "Aee = %.4f;\n", init_Aee);
    fprintf(fid, "Uee = %.4f;\n", init_Uee);
    fprintf(fid, "tau_recee = %.4f;\n", init_tau_recee);
    fprintf(fid, "tau_facilee = %.4f;\n\n", init_tau_facilee);

    fprintf(fid, "Aei = %.4f;\n", init_Aei);
    fprintf(fid, "Uei = %.4f;\n", init_Uei);
    fprintf(fid, "tau_recei = %.4f;\n", init_tau_recei);
    fprintf(fid, "tau_facilei = %.4f;\n\n", init_tau_facilei);

    fprintf(fid, "Aie = %.4f;\n", init_Aie);
    fprintf(fid, "Uie = %.4f;\n", init_Uie);
    fprintf(fid, "tau_recie = %.4f;\n", init_tau_recie);
    fprintf(fid, "tau_facilie = %.4f;\n\n", init_tau_facilie);

    fprintf(fid, "Aii = %.4f;\n", init_Aii);
    fprintf(fid, "Uii = %.4f;\n", init_Uii);
    fprintf(fid, "tau_recii = %.4f;\n", init_tau_recii);
    fprintf(fid, "tau_facilii = %.4f;\n\n", init_tau_facilii);

    fprintf(fid, "g = %.4f;\n", init_g);
    fprintf(fid, "Ee = %.4f;\n", init_Ee);
    fprintf(fid, "Ei = %.4f;\n", init_Ei);
    fprintf(fid, "tau_ge = %.4f;\n", init_tau_ge);
    fprintf(fid, "tau_gi = %.4f;\n", init_tau_gi);

    fclose(fid);

    return 0;
}

int SynapseTsodyksMarkramRKNest::initSynapses(){
    return initSynapsesLocal();
}

/// Tsodyks-Markram for kostya style

std::string SynapseKostya::synapsetype = \
        "Kostya_synapse";

std::string SynapseKostya::getName(){
    return synapsetype;
}

SynapseKostya::SynapseKostya(){
    toggler = 1;
}

void SynapseKostya::setData(int pre, int pos, int preex, int posex, double dt){
    last_spiked = -100;
    last_spiked2 = -100;
    last_spiked_post = -100;

    if(toggler){
        WP = 0.3;//0.3;
        WM = 0.03;//0.3105;
        w_min = 0;
        taup = 20;
        taum = 20;
        
        presynaptic = pre;
        postsynaptic = pos;

        // STDP:
        if(rand()%2){
            w = vf_distributions::normal(54, 10.8, 21.6, 86.4);
            w_max = w;
        } else {
            w = 0;
            w_max = vf_distributions::normal(54, 10.8, 21.6, 86.4);
        }

        // tm
        if(!preex){
            //exc
            U = 0.5;
            D = 1100;
            F = 50;
            exc = 1;
        } else {
            // inh
            U = 0.25;
            D = 700;
            F = 20;
            exc = 0;
        }
        U = vf_distributions::normal(U, U/10, 0, 1);
        D = vf_distributions::normal(D, D/10, 0, 10*D);
        F = vf_distributions::normal(F, F/10, 0, 10*F);

        tau_one = 1;
        toggler = 0;
    } else {
        if(exc){
            std::default_random_engine generator(rand());
            std::gamma_distribution<double> distribution (90.0/7.0,0.7);
            w = MIN(distribution(generator), w_max);
        }
        toggler = 1;
    }

    u = U;
    R = 1;
    g = 0;
    weight = w;
}

double SynapseKostya::evolve(double dt, double time, double Vpre, double Vpost){

    if(vf_discrete::diracDelta(time - last_spiked, dt)){
        // buffer
        w1 = w; u1 = u; R1 = R;

        // t-m
        h = last_spiked - last_spiked2;
        last_spiked2 = last_spiked;
        u = U + u1 * ( 1 - U ) * exp( - h / F );
        R = 1 + ( R1 - u1 * R1 - 1) * exp( - h / D);

        // stdp
        if(toggler && exc){
            h = last_spiked_post - last_spiked;
            if(h > 0){
                w = MIN( w_max , w1 + WP * exp( - h / taup ));
            } else {
                w = MAX( w_min , w1 - WM * exp( h / taup ));
            }
        }

        A = w * u * R;
        if(exc)
            g += A;
        else
            g -= A;
    }

    g -= dt * (g / tau_one);
    out_current = g;

    weight = w;

    return moveDeliveries();
}

double* SynapseKostya::exportData(){
    double* arr = new double[numEssentialVariables()+1];
    arr[0] = numEssentialVariables()+1;
    arr[1] = 1;
    return arr;
}

int SynapseKostya::importData(double *arr){
    return 0;
}

int SynapseKostya::numEssentialVariables(){
    return 1;
}
int SynapseKostya::initSynapses(){
    return 0;
}


/// Spike-Timing Dependent Plasticity Model with transmission function (g)

std::string SynapseSTDPG::synapsetype = "STDP_with_g-function";

double SynapseSTDPG::init_g = 0.0;
double SynapseSTDPG::init_tau_s = 5;
double SynapseSTDPG::init_lambda = 0.1;
double SynapseSTDPG::init_alpha = 1;
double SynapseSTDPG::init_tau_corr = 20;
double SynapseSTDPG::init_weight = 0.5;

std::string SynapseSTDPG::getName(){
    return synapsetype;
}

void SynapseSTDPG::setData(int pre, int post, int preex, int posex, double dt){
    presynaptic = pre;
    postsynaptic = post;
    last_spiked = -100;
    last_spiked_post = -100;

    g = init_g;
    tau_s = init_tau_s;

    // STDP rule:
    // Dt = t_post - t_pre,
    // Dw = lambda*(1-w)^mu*exp(-|Dt|/tau), if Dt > 0,
    // Dw = -lambda*alpha*w^mu*exp(-|Dt|/tau), if Dt <= 0,
    // where 0 < lambda << 1, 0 <= mu <= 1, alpha ~ 1.
    // tau_s * dg/dt = -g + sum_total (delta(t-t_j^pre),t_i<t),
    //    mu = 0.5;

    lambda = init_lambda;
    alpha = init_alpha;
    tau_corr = init_tau_corr; //ms

    weight = init_weight;
}

double SynapseSTDPG::evolve(double dt, double time, double Vpre, double Vpost){

    if(last_spiked>0 && last_spiked_post>0 && (\
            vf_discrete::diracDelta(last_spiked_post - time, dt) || \
            vf_discrete::diracDelta(last_spiked - time, dt)\
            ) ) {
        if(last_spiked_post > last_spiked){
            weight += lambda * sqrt(1-weight)*\
                    exp(-abs(last_spiked_post - last_spiked)/tau_corr);
        } else {
            weight -= lambda * alpha * sqrt(weight) * \
                    exp(-abs(last_spiked_post - last_spiked)/tau_corr);
        }
    }
    if(weight < 1e-8) weight = 1e-8; //bug? shows NaN else
    if(weight > 1) weight = 1; //calc mistakes
    g += dt*(-g/tau_s) + vf_discrete::diracDelta(last_spiked - time, dt);
    out_current = weight * g;

    return moveDeliveries();
}

double* SynapseSTDPG::exportData(){
    double* arr = new double[numEssentialVariables()+1];
    arr[0] = numEssentialVariables()+1;
    arr[1] = presynaptic;
    arr[2] = postsynaptic;
    arr[3] = weight;
    arr[4] = delay;
    arr[5] = g;
    arr[6] = tau_s;
    arr[7] = lambda;
    arr[8] = alpha;
    arr[9] = tau_corr;
    return arr;
}

int SynapseSTDPG::importData(double *arr){
    presynaptic = arr[0];
    postsynaptic = arr[1];
    weight = arr[2];
    delay = arr[3];
    g = arr[4];
    tau_s = arr[5];
    lambda = arr[6];
    alpha = arr[7];
    tau_corr = arr[8];

    if(tau_corr < 1e-5)
        tau_corr = 1e-5;

    last_spiked = -100;
    last_spiked_post = -100;
    return 0;
}

int SynapseSTDPG::numEssentialVariables(){
    return 9;
}

int SynapseSTDPG::initSynapsesLocal(){
    using namespace vf_file;
    std::string buf30 = loadFileToString("./init/SynapseSTDPG.ini");

    init_g = getParameterIni("g", buf30);
    init_tau_s = getParameterIni("tau_s", buf30);
    init_lambda = getParameterIni("lambda", buf30);
    init_alpha = getParameterIni("alpha", buf30);
    init_tau_corr = getParameterIni("tau_corr", buf30);
    init_weight = getParameterIni("weight", buf30);

    FILE* fid = fopen (getFilenameFromIni(DATAFILES, FILE_EXP_INIT_SYNAPSES).c_str(), "w");
    if(!fid){
        exit(15);
    }
    fprintf(fid, "g = %.4f\n", init_g);
    fprintf(fid, "tau_s = %.4f\n", init_tau_s);
    fprintf(fid, "lambda = %.4f\n", init_lambda);
    fprintf(fid, "alpha = %.4f\n", init_alpha);
    fprintf(fid, "tau_corr = %.4f\n", init_tau_corr);
    fprintf(fid, "weight = %.4f\n", init_weight);

    fclose(fid);

    return 0;
}

int SynapseSTDPG::initSynapses(){
    return initSynapsesLocal();
}

/// STDP on exc synapses AND T-M

std::string SynapseTMexcSTDP::synapsetype = "Tsodyks-Markram_and_exc_syn_STDP_model_(no2mult)";

double SynapseTMexcSTDP::init_lambda = 0.1;
double SynapseTMexcSTDP::init_alpha = 1;
double SynapseTMexcSTDP::init_tau_corr = 20;
double SynapseTMexcSTDP::init_weight = 0.5;
double SynapseTMexcSTDP::init_tau_one = 3;
double SynapseTMexcSTDP::init_x = 0.98;
double SynapseTMexcSTDP::init_y = 0.01;
double SynapseTMexcSTDP::init_z = 0.01;
double SynapseTMexcSTDP::init_Aee = 38;
double SynapseTMexcSTDP::init_Uee = 0.5;
double SynapseTMexcSTDP::init_tau_recee = 800;
double SynapseTMexcSTDP::init_tau_facilee = 0;
double SynapseTMexcSTDP::init_Aei = 54;
double SynapseTMexcSTDP::init_Uei = 0.5;
double SynapseTMexcSTDP::init_tau_recei = 800;
double SynapseTMexcSTDP::init_tau_facilei = 0;
double SynapseTMexcSTDP::init_Aie = -72;
double SynapseTMexcSTDP::init_Uie = 0.04;
double SynapseTMexcSTDP::init_tau_recie = 100;
double SynapseTMexcSTDP::init_tau_facilie = 1000;
double SynapseTMexcSTDP::init_Aii = -72;
double SynapseTMexcSTDP::init_Uii = 0.04;
double SynapseTMexcSTDP::init_tau_recii = 100;
double SynapseTMexcSTDP::init_tau_facilii = 1000;
double SynapseTMexcSTDP::init_t_start = 5000;
int SynapseTMexcSTDP::init_type_of_weight = 0;

std::string SynapseTMexcSTDP::getName(){
    return synapsetype;
}

void SynapseTMexcSTDP::setData(int pre, int pos, \
                                    int preex, int posex, double dt){
    tau_one = init_tau_one;
    t_start = init_t_start;

    presynaptic = pre;
    postsynaptic = pos;
    weight = 1;
    last_spiked = -100;
    last_spiked_post = -100;

    if(preex && posex){
        //exc exc
        A = init_Aee;
        U = init_Uee;
        tau_rec = init_tau_recee;
        tau_facil = init_tau_facilee;
        exc = 1;
    } else if(preex){
        //exc inh
        A = init_Aei;
        U = init_Uei;
        tau_rec = init_tau_recei;
        tau_facil = init_tau_facilei;
        exc = 1;
    } else if(posex){
        // inh exc
        A = init_Aie;
        U = init_Uie;
        tau_rec = init_tau_recie;
        tau_facil = init_tau_facilie;
        exc = 0;
    } else {
        // inh inh
        A = init_Aii;
        U = init_Uii;
        tau_rec = init_tau_recii;
        tau_facil = init_tau_facilii;
        exc = 0;
    }
    if(A>0)
        A = vf_distributions::normal(A, A/2, 0, 4*A);
    else
        A = vf_distributions::normal(A, -A/2, 4*A, 0);
    U = vf_distributions::normal(U, U/2, 0, MIN(1, 4*U));
    tau_rec = vf_distributions::normal(tau_rec, tau_rec/2, dt, tau_rec*4);
    tau_facil = vf_distributions::normal(tau_facil, tau_facil/2, \
                                        dt, tau_facil*4);

    x = init_x;
    y = init_y;
    z = init_z;
    u = U;

    // STDP rule:
    // Dt = t_post - t_pre,
    // Dw = lambda*(1-w)^mu*exp(-|Dt|/tau), if Dt > 0,
    // Dw = -lambda*alpha*w^mu*exp(-|Dt|/tau), if Dt <= 0,
    // where 0 < lambda << 1, 0 <= mu <= 1, alpha ~ 1.
    // tau_s * dg/dt = -g + sum_total (delta(t-t_j^pre),t_i<t),
    //    mu = 1;

    lambda = init_lambda;
    alpha = init_alpha;
    tau_corr = init_tau_corr; //ms
    switch(init_type_of_weight){
    case 0:
        weight = init_weight;
        break;
    case 1:
        weight = vf_distributions::uniform(0, 1);
        break;
    case 2:
        weight = vf_distributions::normal(init_weight, init_weight*0.1, 0, 1);
        break;
    default:
        weight = init_weight;
    }

}

double SynapseTMexcSTDP::Xr(double x1, double y1, double z1, double u1, \
                                  double t1, double dt){
    return z1/tau_rec - vf_discrete::diracDelta(t1 - last_spiked, dt) * u1 * x1;
}

double SynapseTMexcSTDP::Yr(double x1, double y1, double z1, double u1, \
                                  double t1, double dt){
    return -y1/tau_one + vf_discrete::diracDelta(t1 - last_spiked, dt) * u1 * x1;
}

double SynapseTMexcSTDP::Zr(double x1, double y1, double z1, double u1, \
                                  double t1, double dt){
    return -z1/tau_rec + y1/tau_one;
}

double SynapseTMexcSTDP::Ur(double x1, double y1, double z1, double u1, \
                                  double t1, double dt){
    if(exc)
        return U;
    else
        return -u1/tau_facil + vf_discrete::diracDelta(t1 - last_spiked, dt) \
                * U * (1 - u1);
}

double SynapseTMexcSTDP::evolve(double dt, double time, \
                                       double Vpre, double Vpost){

    k_1_x = dt * Xr(x, y, z, u, time, dt);
    k_1_y = dt * Yr(x, y, z, u, time, dt);
    k_1_z = dt * Zr(x, y, z, u, time, dt);
    if(exc)
        k_1_u = dt * U;
    else
        k_1_u = dt * Ur(x, y, z, u, time, dt);

    k_2_x = dt * Xr(x + k_1_x/2.0, y + k_1_y/2.0, z + k_1_z/2.0, \
                    u + k_1_u/2.0, time + dt/2.0, dt);
    k_2_y = dt * Yr(x + k_1_x/2.0, y + k_1_y/2.0, z + k_1_z/2.0, \
                    u + k_1_u/2.0, time + dt/2.0, dt);
    k_2_z = dt * Zr(x + k_1_x/2.0, y + k_1_y/2.0, z + k_1_z/2.0, \
                    u + k_1_u/2.0, time + dt/2.0, dt);
    if(exc)
        k_2_u = dt * k_1_u;
    else
        k_2_u = dt * Ur(x + k_1_x/2.0, y + k_1_y/2.0, z + k_1_z/2.0, \
                    u + k_1_u/2.0, time + dt/2.0, dt);

    k_3_x = dt * Xr(x + k_2_x/2.0, y + k_2_y/2.0, z + k_2_z/2.0, \
                    u + k_2_u/2.0, time + dt/2.0, dt);
    k_3_y = dt * Yr(x + k_2_x/2.0, y + k_2_y/2.0, z + k_2_z/2.0, \
                    u + k_2_u/2.0, time + dt/2.0, dt);
    k_3_z = dt * Zr(x + k_2_x/2.0, y + k_2_y/2.0, z + k_2_z/2.0, \
                    u + k_2_u/2.0, time + dt/2.0, dt);
    if(exc)
        k_3_u = dt * k_2_u;
    else
        k_3_u = dt * Ur(x + k_2_x/2.0, y + k_2_y/2.0, z + k_2_z/2.0, \
                    u + k_2_u/2.0, time + dt/2.0, dt);

    x += (k_1_x + 2.0*k_2_x + 2.0*k_3_x + k_4_x)/6.0;
    y += (k_1_y + 2.0*k_2_y + 2.0*k_3_y + k_4_y)/6.0;
    z += (k_1_z + 2.0*k_2_z + 2.0*k_3_z + k_4_z)/6.0;
    if(!exc)
        u += (k_1_u + 2.0*k_2_u + 2.0*k_3_u + k_4_u)/6.0;

    if(x >= 1) x = 0.9999999; if(x<=0) x = 0.0000001;
    if(y >= 1) y = 0.9999999; if(y<=0) y = 0.0000001;
    if(z >= 1) z = 0.9999999; if(z<=0) z = 0.0000001;
    if(u >= 1) u = 0.9999999; if(u<=0) u = 0.0000001;

    if(time>t_start && exc){
        if(last_spiked>0 && last_spiked_post>0 && (\
                vf_discrete::diracDelta(last_spiked_post - time, dt) || \
                vf_discrete::diracDelta(last_spiked - time, dt)
                ) ) {
            if(last_spiked_post > last_spiked){
                weight += lambda * (1-weight)*\
                        exp(-abs(last_spiked_post - last_spiked)/tau_corr);
            } else {
                weight -= lambda * alpha * (weight) * \
                        exp(-abs(last_spiked_post - last_spiked)/tau_corr);
            }
        }
        if(weight < 1e-8) weight = 0; //bug? shows NaN else
        if(weight > 1) weight = 1; //calc mistakes
    }

    out_current = 2 * A * y * weight;

    return moveDeliveries();
}

double* SynapseTMexcSTDP::exportData(){
    double* arr = new double[numEssentialVariables()+1];
    arr[0] = numEssentialVariables()+1;
    arr[1] = presynaptic;
    arr[2] = postsynaptic;
    arr[3] = weight;
    arr[4] = delay;
    arr[5] = x;
    arr[6] = y;
    arr[7] = z;
    arr[8] = u;
    arr[9] = tau_one;
    arr[10] = tau_rec;
    arr[11] = tau_facil;
    arr[12] = U;
    arr[13] = A;
    arr[14] = exc;
    return arr;
}

int SynapseTMexcSTDP::importData(double *arr){
    presynaptic = arr[0];
    postsynaptic = arr[1];
    weight = arr[2];
    delay = arr[3];
    x = arr[4];
    y = arr[5];
    z = arr[6];
    u = arr[7];
    tau_one = arr[8];
    tau_rec = arr[9];
    tau_facil = arr[10];
    U = arr[11];
    A = arr[12];
    exc = arr[13];

    if(tau_rec < 1e-5)
        tau_rec = 1e-5;
    if(tau_facil < 1e-5)
        tau_facil = 1e-5;

    lambda = init_lambda;
    alpha = init_alpha;
    tau_corr = init_tau_corr;
    t_start = init_t_start;
    switch(init_type_of_weight){
    case 0:
        weight = init_weight;
        break;
    case 1:
        weight = vf_distributions::uniform(0, 1);
        break;
    case 2:
        weight = vf_distributions::normal(init_weight, init_weight*0.1, 0, 1);
        break;
    default:
        weight = init_weight;
    }

    last_spiked = -100;
    last_spiked_post = -100;
    return 0;
}

int SynapseTMexcSTDP::numEssentialVariables(){
    return 14;
}

int SynapseTMexcSTDP::initSynapsesLocal(){
    using namespace vf_file;
    std::string buf30 = loadFileToString("./init/SynapseTMSTDP.ini");

    init_tau_one = getParameterIni("tau_one", buf30);
    init_x = getParameterIni("x", buf30);
    init_y = getParameterIni("y", buf30);
    init_z = getParameterIni("z", buf30);

    init_Aee = getParameterIni("Aee", buf30);
    init_Uee = getParameterIni("Uee", buf30);
    init_tau_recee = getParameterIni("tau_recee", buf30);
    init_tau_facilee = getParameterIni("tau_facilee", buf30);

    init_Aei = getParameterIni("Aei", buf30);
    init_Uei = getParameterIni("Uei", buf30);
    init_tau_recei = getParameterIni("tau_recei", buf30);
    init_tau_facilei = getParameterIni("tau_facilei", buf30);

    init_Aie = getParameterIni("Aie", buf30);
    init_Uie = getParameterIni("Uie", buf30);
    init_tau_recie = getParameterIni("tau_recie", buf30);
    init_tau_facilie = getParameterIni("tau_facilie", buf30);

    init_Aii = getParameterIni("Aii", buf30);
    init_Uii = getParameterIni("Uii", buf30);
    init_tau_recii = getParameterIni("tau_recii", buf30);
    init_tau_facilii = getParameterIni("tau_facilii", buf30);

    init_t_start = getParameterIni("t_start", buf30);
    init_type_of_weight = getParameterIni("type_of_weight", buf30);
    init_lambda = getParameterIni("lambda", buf30);
    init_alpha = getParameterIni("alpha", buf30);
    init_tau_corr = getParameterIni("tau_corr", buf30);
    init_weight = getParameterIni("weight", buf30);

    FILE* fid = fopen (getFilenameFromIni(DATAFILES, FILE_EXP_INIT_SYNAPSES).c_str(), "w");
    if(!fid){
        exit(15);
    }
    fprintf(fid, "tau_one = %.4f;\n", init_tau_one);
    fprintf(fid, "x = %.4f;\n", init_x);
    fprintf(fid, "y = %.4f;\n", init_y);
    fprintf(fid, "z = %.4f;\n\n", init_z);

    fprintf(fid, "Aee = %.4f;\n", init_Aee);
    fprintf(fid, "Uee = %.4f;\n", init_Uee);
    fprintf(fid, "tau_recee = %.4f;\n", init_tau_recee);
    fprintf(fid, "tau_facilee = %.4f;\n\n", init_tau_facilee);

    fprintf(fid, "Aei = %.4f;\n", init_Aei);
    fprintf(fid, "Uei = %.4f;\n", init_Uei);
    fprintf(fid, "tau_recei = %.4f;\n", init_tau_recei);
    fprintf(fid, "tau_facilei = %.4f;\n\n", init_tau_facilei);

    fprintf(fid, "Aie = %.4f;\n", init_Aie);
    fprintf(fid, "Uie = %.4f;\n", init_Uie);
    fprintf(fid, "tau_recie = %.4f;\n", init_tau_recie);
    fprintf(fid, "tau_facilie = %.4f;\n\n", init_tau_facilie);

    fprintf(fid, "Aii = %.4f;\n", init_Aii);
    fprintf(fid, "Uii = %.4f;\n", init_Uii);
    fprintf(fid, "tau_recii = %.4f;\n", init_tau_recii);
    fprintf(fid, "tau_facilii = %.4f;\n\n", init_tau_facilii);

    fprintf(fid, "t_start = %.4f;\n", init_t_start);
    fprintf(fid, "type_of_weight = %d;\n", init_type_of_weight);
    fprintf(fid, "lambda = %.4f;\n", init_lambda);
    fprintf(fid, "alpha = %.4f;\n", init_alpha);
    fprintf(fid, "tau_corr = %.4f;\n", init_tau_corr);
    fprintf(fid, "weight = %.4f;\n", init_weight);
    fclose(fid);

    return 0;
}

int SynapseTMexcSTDP::initSynapses(){
    return initSynapsesLocal();
}

/// STDP AND T-M

std::string SynapseTMSTDP::synapsetype = \
        "Tsodyks-Markram_and_STDP_model_(full)";

double SynapseTMSTDP::init_A_plus = 0.8;
double SynapseTMSTDP::init_A_minus = 0.3;
double SynapseTMSTDP::init_tau_corr_plus = 20;
double SynapseTMSTDP::init_tau_corr_minus = 30;
double SynapseTMSTDP::init_weight = 0.5;
double SynapseTMSTDP::init_tau_one = 3;
double SynapseTMSTDP::init_x = 0.98;
double SynapseTMSTDP::init_y = 0.01;
double SynapseTMSTDP::init_z = 0.01;
double SynapseTMSTDP::init_Aee = 38;
double SynapseTMSTDP::init_Uee = 0.5;
double SynapseTMSTDP::init_tau_recee = 800;
double SynapseTMSTDP::init_tau_facilee = 0;
double SynapseTMSTDP::init_Aei = 54;
double SynapseTMSTDP::init_Uei = 0.5;
double SynapseTMSTDP::init_tau_recei = 800;
double SynapseTMSTDP::init_tau_facilei = 0;
double SynapseTMSTDP::init_Aie = -72;
double SynapseTMSTDP::init_Uie = 0.04;
double SynapseTMSTDP::init_tau_recie = 100;
double SynapseTMSTDP::init_tau_facilie = 1000;
double SynapseTMSTDP::init_Aii = -72;
double SynapseTMSTDP::init_Uii = 0.04;
double SynapseTMSTDP::init_tau_recii = 100;
double SynapseTMSTDP::init_tau_facilii = 1000;
double SynapseTMSTDP::init_t_start = 5000;
int SynapseTMSTDP::init_type_of_weight = 0;
int SynapseTMSTDP::init_distribute_params = 0;

std::string SynapseTMSTDP::getName(){
    return synapsetype;
}

void SynapseTMSTDP::setData(int pre, int pos, \
                                    int preex, int posex, double dt){
    tau_one = init_tau_one;
    t_start = init_t_start;

    presynaptic = pre;
    postsynaptic = pos;
    weight = 1;
    last_spiked = -100;
    last_spiked_post = -100;

    if(preex && posex){
        //exc exc
        A = init_Aee;
        U = init_Uee;
        tau_rec = init_tau_recee;
        tau_facil = init_tau_facilee;
        exc = 1;
    } else if(preex){
        //exc inh
        A = init_Aei;
        U = init_Uei;
        tau_rec = init_tau_recei;
        tau_facil = init_tau_facilei;
        exc = 1;
    } else if(posex){
        // inh exc
        A = init_Aie;
        U = init_Uie;
        tau_rec = init_tau_recie;
        tau_facil = init_tau_facilie;
        exc = 0;
    } else {
        // inh inh
        A = init_Aii;
        U = init_Uii;
        tau_rec = init_tau_recii;
        tau_facil = init_tau_facilii;
        exc = 0;
    }

    if(init_distribute_params){
        if(A>0)
            A = vf_distributions::normal(A, A/2, 0, 4*A);
        else
            A = vf_distributions::normal(A, -A/2, 4*A, 0);
        U = vf_distributions::normal(U, U/2, 0, MIN(1, 4*U));
        tau_rec = vf_distributions::normal(tau_rec, tau_rec/2, dt, tau_rec*4);
        tau_facil = vf_distributions::normal(tau_facil, tau_facil/2, \
                                            dt, tau_facil*4);
    }

    x = init_x;
    y = init_y;
    z = init_z;
    u = U;

    // STDP rule:
    // Dt = t_post - t_pre,
    // Dw = A_plus*(1-w)^mu*exp(-|Dt|/tau_plus), if Dt > 0,
    // Dw = -A_minus*w^mu*exp(-|Dt|/tau_minus), if Dt <= 0,
    // where 0 < lambda << 1, 0 <= mu <= 1, alpha ~ 1.
    //    mu = 1;

    A_plus = init_A_plus;
    A_minus = init_A_minus;
    tau_corr_plus = init_tau_corr_plus;
    tau_corr_minus = init_tau_corr_minus;
    switch(init_type_of_weight){
    case 0:
        weight = init_weight;
        break;
    case 1:
        weight = vf_distributions::uniform(0, 1);
        break;
    case 2:
        weight = vf_distributions::normal(init_weight, init_weight*0.1, 0, 1);
        break;
    default:
        weight = init_weight;
    }

}

double SynapseTMSTDP::Xr(double x1, double y1, double z1, \
                                     double u1, double t1, double dt){
    return z1/tau_rec - vf_discrete::diracDelta(t1 - last_spiked, dt) * u1 * x1;
}

double SynapseTMSTDP::Yr(double x1, double y1, double z1, \
                                     double u1, double t1, double dt){
    return -y1/tau_one + vf_discrete::diracDelta(t1 - last_spiked, dt) * u1 * x1;
}

double SynapseTMSTDP::Zr(double x1, double y1, double z1, \
                                     double u1, double t1, double dt){
    return -z1/tau_rec + y1/tau_one;
}

double SynapseTMSTDP::Ur(double x1, double y1, double z1, \
                                     double u1, double t1, double dt){
    if(exc)
        return U;
    else
        return -u1/tau_facil + vf_discrete::diracDelta(t1 - last_spiked, dt) \
                * U * (1 - u1);
}

double SynapseTMSTDP::evolve(double dt, double time, \
                                         double Vpre, double Vpost){

    k_1_x = dt * Xr(x, y, z, u, time, dt);
    k_1_y = dt * Yr(x, y, z, u, time, dt);
    k_1_z = dt * Zr(x, y, z, u, time, dt);
    if(exc)
        k_1_u = dt * U;
    else
        k_1_u = dt * Ur(x, y, z, u, time, dt);

    k_2_x = dt * Xr(x + k_1_x/2.0, y + k_1_y/2.0, z + k_1_z/2.0, \
                    u + k_1_u/2.0, time + dt/2.0, dt);
    k_2_y = dt * Yr(x + k_1_x/2.0, y + k_1_y/2.0, z + k_1_z/2.0, \
                    u + k_1_u/2.0, time + dt/2.0, dt);
    k_2_z = dt * Zr(x + k_1_x/2.0, y + k_1_y/2.0, z + k_1_z/2.0, \
                    u + k_1_u/2.0, time + dt/2.0, dt);
    if(exc)
        k_2_u = dt * k_1_u;
    else
        k_2_u = dt * Ur(x + k_1_x/2.0, y + k_1_y/2.0, z + k_1_z/2.0, \
                    u + k_1_u/2.0, time + dt/2.0, dt);

    k_3_x = dt * Xr(x + k_2_x/2.0, y + k_2_y/2.0, z + k_2_z/2.0, \
                    u + k_2_u/2.0, time + dt/2.0, dt);
    k_3_y = dt * Yr(x + k_2_x/2.0, y + k_2_y/2.0, z + k_2_z/2.0, \
                    u + k_2_u/2.0, time + dt/2.0, dt);
    k_3_z = dt * Zr(x + k_2_x/2.0, y + k_2_y/2.0, z + k_2_z/2.0, \
                    u + k_2_u/2.0, time + dt/2.0, dt);
    if(exc)
        k_3_u = dt * k_2_u;
    else
        k_3_u = dt * Ur(x + k_2_x/2.0, y + k_2_y/2.0, z + k_2_z/2.0, \
                    u + k_2_u/2.0, time + dt/2.0, dt);

    k_4_x = dt * Xr(x + k_3_x, y + k_3_y, z + k_3_z, u + k_3_u, time + dt, dt);
    k_4_y = dt * Yr(x + k_3_x, y + k_3_y, z + k_3_z, u + k_3_u, time + dt, dt);
    k_4_z = dt * Zr(x + k_3_x, y + k_3_y, z + k_3_z, u + k_3_u, time + dt, dt);
    if(exc)
        k_4_u = dt * U;
    else
        k_4_u = dt * Ur(x + k_3_x, y + k_3_y, z + k_3_z, u + k_3_u, time + dt, dt);


    x += (k_1_x + 2.0*k_2_x + 2.0*k_3_x + k_4_x)/6.0;
    y += (k_1_y + 2.0*k_2_y + 2.0*k_3_y + k_4_y)/6.0;
    z += (k_1_z + 2.0*k_2_z + 2.0*k_3_z + k_4_z)/6.0;
    if(!exc)
        u += (k_1_u + 2.0*k_2_u + 2.0*k_3_u + k_4_u)/6.0;

    x += (k_1_x + 2.0*k_2_x + 2.0*k_3_x + k_4_x)/6.0;
    y += (k_1_y + 2.0*k_2_y + 2.0*k_3_y + k_4_y)/6.0;
    z += (k_1_z + 2.0*k_2_z + 2.0*k_3_z + k_4_z)/6.0;
    u += (k_1_u + 2.0*k_2_u + 2.0*k_3_u + k_4_u)/6.0;

    if(x >= 1) x = 0.9999999; if(x<=0) x = 0.0000001;
    if(y >= 1) y = 0.9999999; if(y<=0) y = 0.0000001;
    if(z >= 1) z = 0.9999999; if(z<=0) z = 0.0000001;
    if(u >= 1) u = 0.9999999; if(u<=0) u = 0.0000001;

    if(time>t_start && exc){
        if(last_spiked>0 && last_spiked_post>0 && (\
                vf_discrete::diracDelta(last_spiked_post - time, dt) || \
                vf_discrete::diracDelta(last_spiked - time, dt)
                ) ) {
            if(last_spiked_post > last_spiked){
                weight += A_plus * (1-weight)*\
                       exp(-abs(last_spiked_post - last_spiked)/tau_corr_plus);
            } else {
                weight -= A_minus * (weight) * \
                       exp(-abs(last_spiked_post - last_spiked)/tau_corr_minus);
            }
        }
        if(weight < 1e-8) weight = 0; //bug? shows NaN else
        if(weight > 1) weight = 1; //calc mistakes
    }

    out_current = 2 * A * y * weight;

    return moveDeliveries();
}

double* SynapseTMSTDP::exportData(){
    double* arr = new double[numEssentialVariables()+1];
    arr[0] = numEssentialVariables()+1;
    arr[1] = presynaptic;
    arr[2] = postsynaptic;
    arr[3] = weight;
    arr[4] = delay;
    arr[5] = x;
    arr[6] = y;
    arr[7] = z;
    arr[8] = u;
    arr[9] = tau_one;
    arr[10] = tau_rec;
    arr[11] = tau_facil;
    arr[12] = U;
    arr[13] = A;
    arr[14] = exc;
    return arr;
}

int SynapseTMSTDP::importData(double *arr){
    presynaptic = arr[0];
    postsynaptic = arr[1];
    weight = arr[2];
    delay = arr[3];
    x = arr[4];
    y = arr[5];
    z = arr[6];
    u = arr[7];
    tau_one = arr[8];
    tau_rec = arr[9];
    tau_facil = arr[10];
    U = arr[11];
    A = arr[12];
    exc = arr[13];

    if(tau_rec < 1e-5)
        tau_rec = 1e-5;
    if(tau_facil < 1e-5)
        tau_facil = 1e-5;

    A_plus = init_A_plus;
    A_minus = init_A_minus;
    tau_corr_plus = init_tau_corr_plus;
    tau_corr_minus = init_tau_corr_minus;
    t_start = init_t_start;
    switch(init_type_of_weight){
    case 0:
        weight = init_weight;
        break;
    case 1:
        weight = vf_distributions::uniform(0, 1);
        break;
    case 2:
        weight = vf_distributions::normal(init_weight, init_weight*0.1, 0, 1);
        break;
    default:
        weight = init_weight;
    }

    last_spiked = -100;
    last_spiked_post = -100;
    return 0;
}

int SynapseTMSTDP::numEssentialVariables(){
    return 14;
}

int SynapseTMSTDP::initSynapsesLocal(){
    using namespace vf_file;
    std::string buf30 = loadFileToString("./init/SynapseTMSTDP.ini");

    init_tau_one = getParameterIni("tau_one", buf30);
    init_x = getParameterIni("INIT_x", buf30);
    init_y = getParameterIni("INIT_y", buf30);
    init_z = getParameterIni("INIT_z", buf30);

    init_Aee = getParameterIni("Aee", buf30);
    init_Uee = getParameterIni("Uee", buf30);
    init_tau_recee = getParameterIni("tau_recee", buf30);
    init_tau_facilee = getParameterIni("tau_facilee", buf30);

    init_Aei = getParameterIni("Aei", buf30);
    init_Uei = getParameterIni("Uei", buf30);
    init_tau_recei = getParameterIni("tau_recei", buf30);
    init_tau_facilei = getParameterIni("tau_facilei", buf30);

    init_Aie = getParameterIni("Aie", buf30);
    init_Uie = getParameterIni("Uie", buf30);
    init_tau_recie = getParameterIni("tau_recie", buf30);
    init_tau_facilie = getParameterIni("tau_facilie", buf30);

    init_Aii = getParameterIni("Aii", buf30);
    init_Uii = getParameterIni("Uii", buf30);
    init_tau_recii = getParameterIni("tau_recii", buf30);
    init_tau_facilii = getParameterIni("tau_facilii", buf30);

    init_t_start = getParameterIni("t_start", buf30);
    init_type_of_weight = getParameterIni("type_of_weight", buf30);
    init_weight = getParameterIni("Weight", buf30);

    init_A_plus = getParameterIni("A_plus", buf30);
    init_A_minus = getParameterIni("A_minus", buf30);
    init_tau_corr_plus = getParameterIni("tau_corr_plus", buf30);
    init_tau_corr_minus = getParameterIni("tau_corr_minus", buf30);

    init_distribute_params = getParameterIni("DISTRIBUTE_PARAMETERS", buf30);

    FILE* fid = fopen (getFilenameFromIni(DATAFILES, FILE_EXP_INIT_SYNAPSES).c_str(), "w");
    if(!fid){
        exit(15);
    }
    fprintf(fid, "tau_one = %.4f;\n", init_tau_one);
    fprintf(fid, "INIT_x = %.4f;\n", init_x);
    fprintf(fid, "INIT_y = %.4f;\n", init_y);
    fprintf(fid, "INIT_z = %.4f;\n\n", init_z);

    fprintf(fid, "Aee = %.4f;\n", init_Aee);
    fprintf(fid, "Uee = %.4f;\n", init_Uee);
    fprintf(fid, "tau_recee = %.4f;\n", init_tau_recee);
    fprintf(fid, "tau_facilee = %.4f;\n\n", init_tau_facilee);

    fprintf(fid, "Aei = %.4f;\n", init_Aei);
    fprintf(fid, "Uei = %.4f;\n", init_Uei);
    fprintf(fid, "tau_recei = %.4f;\n", init_tau_recei);
    fprintf(fid, "tau_facilei = %.4f;\n\n", init_tau_facilei);

    fprintf(fid, "Aie = %.4f;\n", init_Aie);
    fprintf(fid, "Uie = %.4f;\n", init_Uie);
    fprintf(fid, "tau_recie = %.4f;\n", init_tau_recie);
    fprintf(fid, "tau_facilie = %.4f;\n\n", init_tau_facilie);

    fprintf(fid, "Aii = %.4f;\n", init_Aii);
    fprintf(fid, "Uii = %.4f;\n", init_Uii);
    fprintf(fid, "tau_recii = %.4f;\n", init_tau_recii);
    fprintf(fid, "tau_facilii = %.4f;\n\n", init_tau_facilii);

    fprintf(fid, "t_start = %.4f;\n", init_t_start);
    fprintf(fid, "type_of_weight = %d;\n", init_type_of_weight);
    fprintf(fid, "Weight = %.4f;\n\n", init_weight);

    fprintf(fid, "A_plus = %.4f;\n", init_A_plus);
    fprintf(fid, "A_minus = %.4f;\n", init_A_minus);
    fprintf(fid, "tau_corr_plus = %.4f;\n", init_tau_corr_plus);
    fprintf(fid, "tau_corr_minus = %.4f;\n", init_tau_corr_minus);

    fprintf(fid, "DISTRIBUTE_PARAMETERS = %d;\n", init_distribute_params);
    fclose(fid);

    return 0;
}

int SynapseTMSTDP::initSynapses(){
    return initSynapsesLocal();
}

double SynapseTMSTDP::test(){
    return init_type_of_weight;
}


//Simple synapse with no short- or long-term plasticity.

std::string SynapseSimpleG::synapsetype = "Simple_g-synapse";

double SynapseSimpleG::init_g = 0.0;
double SynapseSimpleG::init_tau_se = 3;
double SynapseSimpleG::init_tau_si = 7;
double SynapseSimpleG::init_g_jump = 0.2;
double SynapseSimpleG::init_A = 30;
int SynapseSimpleG::TYPE = 0;

double* SynapseSimpleG::exportData(){
    double* arr = new double[numEssentialVariables()+1];
    arr[0] = numEssentialVariables()+1;
    arr[1] = presynaptic;
    arr[2] = postsynaptic;
    arr[3] = weight;
    arr[4] = delay;
    arr[5] = g;
    arr[6] = tau_s;
    arr[7] = g_jump;
    arr[8] = A;
    arr[9] = max_spikes;
    return arr;
}

int SynapseSimpleG::importData(double *arr){
    presynaptic =   arr[0];
    postsynaptic =  arr[1];
    weight =        arr[2];
    delay =         arr[3];
    g =             arr[4];
    tau_s =         arr[5];
    g_jump =        arr[6];
    A =             arr[7];
    max_spikes =    arr[8];
    if(tau_s < 1e-5)
        tau_s = 1e-5;
    last_spiked = -100;
    last_spiked_post = -100;
    return 0;
}

int SynapseSimpleG::numEssentialVariables(){
    return 9;
}

int SynapseSimpleG::initSynapsesLocal(){
    using namespace vf_file;
    std::string buf30 = loadFileToString("./init/SynapseGSimple.ini");

    TYPE = getParameterIni("TYPE", buf30);
    init_g = getParameterIni("INIT_g", buf30);
    init_A = getParameterIni("INIT_A", buf30);
    init_tau_se = getParameterIni("tau_se", buf30);
    init_tau_si = getParameterIni("tau_si", buf30);
    init_g_jump = getParameterIni("g_jump", buf30);

    FILE* fid = fopen (getFilenameFromIni(DATAFILES, FILE_EXP_INIT_SYNAPSES).c_str(), "w");
    if(!fid){
        exit(15);
    }
    fprintf(fid, "TYPE = %d;\n", TYPE);
    fprintf(fid, "INIT_g = %.4f;\n", init_g);
    fprintf(fid, "INIT_A = %.4f;\n", init_A);
    fprintf(fid, "tau_se = %.4f;\n", init_tau_se);
    fprintf(fid, "tau_si = %.4f;\n", init_tau_si);
    fprintf(fid, "g_jump = %.4f;\n", init_g_jump);
    fclose(fid);
    
    return 0;
}

int SynapseSimpleG::initSynapses(){
    return initSynapsesLocal();
}

std::string SynapseSimpleG::getName(){
    return synapsetype;
}

void SynapseSimpleG::setData(int pre, int pos, int preex, int posex, \
                                 double dt){
    presynaptic = pre;
    postsynaptic = pos;
    weight = 1;
    last_spiked = -100;
    last_spiked_post = -100;
    g_jump = init_g_jump;
    g = init_g;
    A = init_A;

    if(preex){
        tau_s = init_tau_se;
    } else {
        tau_s = init_tau_si;
    }

    now_spikes = 0;
    max_spikes = (tau_s*6.6383 + 1)/3; //accuracy ~0.01
    if(TYPE == 1){
        spikerow = Malloc(max_spikes, double);
        for(int i=0; i<max_spikes; i++){
            spikerow[i] = 0;
        }
    } else 
        max_spikes = 0;

}

double SynapseSimpleG::alphaFunction(double x){
    return x/tau_s * exp(-x/tau_s);
}

double SynapseSimpleG::evolve(double dt, double timen, double Vpre, \
                                  double Vpost){

    if(TYPE==1){
        if(timen - last_spiked < dt/2){
            spikerow[now_spikes%max_spikes] = timen;
            now_spikes++;
        }

        g = 0.0;
        for(int i=0; i<now_spikes && i<max_spikes; i++)
            g += alphaFunction(timen - spikerow[i]);
        g = g * g_jump;
    } else {
        g += -g / tau_s + g_jump*vf_discrete::diracDelta(timen - last_spiked, dt);
    }

    out_current = g * A;

    return moveDeliveries();
}


/// Models with g
// First type. No cut on post, differential equation.

std::string SynapseGFirstType::synapsetype = "G-model_1st_type";

double SynapseGFirstType::init_g = 0.0;
double SynapseGFirstType::init_tau_se = 3;
double SynapseGFirstType::init_tau_si = 7;
double SynapseGFirstType::init_gs = 0.4;
double SynapseGFirstType::init_Ee = 0;
double SynapseGFirstType::init_Ei = -80;

double* SynapseGFirstType::exportData(){
    double* arr = new double[9];
    arr[0] = 9;
    arr[1] = presynaptic;
    arr[2] = postsynaptic;
    arr[3] = weight;
    arr[4] = delay;
    arr[5] = g;
    arr[6] = tau_s;
    arr[7] = gs;
    arr[8] = E;
    return arr;
}

int SynapseGFirstType::importData(double *arr){
    presynaptic =   arr[0];
    postsynaptic =  arr[1];
    weight =        arr[2];
    delay =         arr[3];
    g =             arr[4];
    tau_s =         arr[5];
    gs =            arr[6];
    E =             arr[7];

    if(tau_s < 1e-5)
        tau_s = 1e-5;

    last_spiked = -100;
    last_spiked_post = -100;
    return 0;
}

int SynapseGFirstType::numEssentialVariables(){
    return 8;
}

int SynapseGFirstType::initSynapsesLocal(){
    using namespace vf_file;
    std::string buf30 = loadFileToString("./init/SynapseG.ini");

    init_g = getParameterIni("INIT_g", buf30);
    init_tau_se = getParameterIni("tau_se", buf30);
    init_tau_si = getParameterIni("tau_si", buf30);
    init_gs = getParameterIni("g_s", buf30);
    init_Ee = getParameterIni("Ee", buf30);
    init_Ei = getParameterIni("Ei", buf30);

    FILE* fid = fopen (getFilenameFromIni(DATAFILES, FILE_EXP_INIT_SYNAPSES).c_str(), "w");
    if(!fid){
        exit(15);
    }
    fprintf(fid, "g = %.4f;\n", init_g);
    fprintf(fid, "tau_se = %.4f;\n", init_tau_se);
    fprintf(fid, "tau_si = %.4f;\n", init_tau_si);
    fprintf(fid, "g_s = %.4f;\n", init_gs);
    fprintf(fid, "Ee = %.4f;\n", init_Ee);
    fprintf(fid, "Ei = %.4f;\n", init_Ei);
    fclose(fid);

    return 0;
}

int SynapseGFirstType::initSynapses(){
    return initSynapsesLocal();
}

std::string SynapseGFirstType::getName(){
    return synapsetype;
}

void SynapseGFirstType::setData(int pre, int pos, int preex, int posex, \
                                double dt){
    presynaptic = pre;
    postsynaptic = pos;
    weight = vf_distributions::uniform(0, 1);
    last_spiked = -100;
    last_spiked_post = -100;
    gs = init_gs;
    g = init_g;

    if(preex){
        tau_s = init_tau_se;
        E = init_Ee;
    } else {
        tau_s = init_tau_si;
        E = init_Ei;
    }
}

double SynapseGFirstType::evolve(double dt, double timen, \
                                 double Vpre, double Vpost){
    g += dt * (-g / tau_s) + gs * vf_discrete::diracDelta(timen-last_spiked, dt);
    out_current = g * (E - Vpost);

    return moveDeliveries();
}
// First type. Cut on post, differential equation.

std::string SynapseGFirstTypeWCUT::synapsetype = "G-model_1st_type_with_cut";

double SynapseGFirstTypeWCUT::init_g = 0.0;
double SynapseGFirstTypeWCUT::init_tau_se = 3;
double SynapseGFirstTypeWCUT::init_tau_si = 7;
double SynapseGFirstTypeWCUT::init_gs = 0.4;
double SynapseGFirstTypeWCUT::init_Ee = 0;
double SynapseGFirstTypeWCUT::init_Ei = -80;

double* SynapseGFirstTypeWCUT::exportData(){
    double* arr = new double[numEssentialVariables()+1];
    arr[0] = numEssentialVariables()+1;
    arr[1] = presynaptic;
    arr[2] = postsynaptic;
    arr[3] = weight;
    arr[4] = delay;
    arr[5] = g;
    arr[6] = tau_s;
    arr[7] = gs;
    arr[8] = E;
    return arr;
}

int SynapseGFirstTypeWCUT::importData(double *arr){
    presynaptic =   arr[0];
    postsynaptic =  arr[1];
    weight =        arr[2];
    delay =         arr[3];
    g =             arr[4];
    tau_s =         arr[5];
    gs =            arr[6];
    E =             arr[7];
    if(tau_s < 1e-5)
        tau_s = 1e-5;
    last_spiked = -100;
    last_spiked_post = -100;
    return 0;
}

int SynapseGFirstTypeWCUT::numEssentialVariables(){
    return 8;
}

int SynapseGFirstTypeWCUT::initSynapsesLocal(){
    using namespace vf_file;
    std::string buf30 = loadFileToString("./init/SynapseG.ini");

    init_g = getParameterIni("INIT_g", buf30);
    init_tau_se = getParameterIni("tau_se", buf30);
    init_tau_si = getParameterIni("tau_si", buf30);
    init_gs = getParameterIni("g_s", buf30);
    init_Ee = getParameterIni("Ee", buf30);
    init_Ei = getParameterIni("Ei", buf30);

    FILE* fid = fopen (getFilenameFromIni(DATAFILES, FILE_EXP_INIT_SYNAPSES).c_str(), "w");
    if(!fid){
        exit(15);
    }
    fprintf(fid, "g = %.4f;\n", init_g);
    fprintf(fid, "tau_se = %.4f;\n", init_tau_se);
    fprintf(fid, "tau_si = %.4f;\n", init_tau_si);
    fprintf(fid, "g_s = %.4f;\n", init_gs);
    fprintf(fid, "Ee = %.4f;\n", init_Ee);
    fprintf(fid, "Ei = %.4f;\n", init_Ei);
    fclose(fid);

    return 0;
}

int SynapseGFirstTypeWCUT::initSynapses(){
    return initSynapsesLocal();
}

std::string SynapseGFirstTypeWCUT::getName(){
    return synapsetype;
}

void SynapseGFirstTypeWCUT::setData(int pre, int pos, int preex, int posex, \
                                    double dt){
    presynaptic = pre;
    postsynaptic = pos;
    weight = vf_distributions::uniform(0, 1);
    last_spiked = -100;
    last_spiked_post = -100;
    gs = init_gs;
    g = init_g;

    if(preex){
        tau_s = init_tau_se;
        E = init_Ee;
    } else {
        tau_s = init_tau_si;
        E = init_Ei;
    }
}

double SynapseGFirstTypeWCUT::evolve(double dt, double timen, double Vpre, \
                                     double Vpost){
    if(timen - last_spiked_post<dt/2){
        g = 0;
    }
    g += dt * (-g / tau_s) + gs * vf_discrete::diracDelta(timen-last_spiked, dt);
    out_current = g * (E - Vpost);

    return moveDeliveries();
}

//Second type. No cut on post, normal equation.

std::string SynapseGSecondType::synapsetype = "G-model_2nd_type";

double SynapseGSecondType::init_g = 0.0;
double SynapseGSecondType::init_tau_se = 3;
double SynapseGSecondType::init_tau_si = 7;
double SynapseGSecondType::init_gs = 0.4;
double SynapseGSecondType::init_Ee = 0;
double SynapseGSecondType::init_Ei = -80;

double* SynapseGSecondType::exportData(){
    double* arr = new double[numEssentialVariables()+1];
    arr[0] = numEssentialVariables()+1;
    arr[1] = presynaptic;
    arr[2] = postsynaptic;
    arr[3] = weight;
    arr[4] = delay;
    arr[5] = g;
    arr[6] = tau_s;
    arr[7] = gs;
    arr[8] = E;
    return arr;
}

int SynapseGSecondType::importData(double *arr){
    presynaptic =   arr[0];
    postsynaptic =  arr[1];
    weight =        arr[2];
    delay =         arr[3];
    g =             arr[4];
    tau_s =         arr[5];
    gs =            arr[6];
    E =             arr[7];
    if(tau_s < 1e-5)
        tau_s = 1e-5;
    last_spiked = -100;
    last_spiked_post = -100;
    return 0;
}

int SynapseGSecondType::numEssentialVariables(){
    return 8;
}

int SynapseGSecondType::initSynapsesLocal(){
    using namespace vf_file;
    std::string buf30 = loadFileToString("./init/SynapseG.ini");

    init_g = getParameterIni("INIT_g", buf30);
    init_tau_se = getParameterIni("tau_se", buf30);
    init_tau_si = getParameterIni("tau_si", buf30);
    init_gs = getParameterIni("g_s", buf30);
    init_Ee = getParameterIni("Ee", buf30);
    init_Ei = getParameterIni("Ei", buf30);

    FILE* fid = fopen (getFilenameFromIni(DATAFILES, FILE_EXP_INIT_SYNAPSES).c_str(), "w");
    if(!fid){
        exit(15);
    }
    fprintf(fid, "g = %.4f;\n", init_g);
    fprintf(fid, "tau_se = %.4f;\n", init_tau_se);
    fprintf(fid, "tau_si = %.4f;\n", init_tau_si);
    fprintf(fid, "g_s = %.4f;\n", init_gs);
    fprintf(fid, "Ee = %.4f;\n", init_Ee);
    fprintf(fid, "Ei = %.4f;\n", init_Ei);
    fclose(fid);
    
    return 0;
}

int SynapseGSecondType::initSynapses(){
    return initSynapsesLocal();
}


std::string SynapseGSecondType::getName(){
    return synapsetype;
}

void SynapseGSecondType::setData(int pre, int pos, int preex, int posex, \
                                 double dt){
    presynaptic = pre;
    postsynaptic = pos;
    weight = vf_distributions::uniform(0, 1);
    last_spiked = -100;
    last_spiked_post = -100;
    gs = init_gs;
    g = init_g;

    if(preex){
        tau_s = init_tau_se;
        E = init_Ee;
    } else {
        tau_s = init_tau_si;
        E = init_Ei;
    }

    now_spikes = 0;
    max_spikes = (tau_s*6.6383 + 1)/dt; //accuracy 0.01
    spikerow = Malloc(max_spikes, double);//
    for(int i=0; i<max_spikes; i++)
        spikerow[i] = 0;

}

double SynapseGSecondType::alphaFunction(double x){
    return x/tau_s * exp(-x/tau_s);
}

double SynapseGSecondType::evolve(double dt, double timen, double Vpre, \
                                  double Vpost){

    if(timen - last_spiked < dt/2){
        spikerow[now_spikes%max_spikes] = timen;
        now_spikes++;
    }

    g = 0.0;
    for(int i=0; i<now_spikes && i<max_spikes; i++)
        g += alphaFunction(timen - spikerow[i]);

    out_current = g * (E - Vpost) * gs;

    return moveDeliveries();
}

//Second type. Cut on post, normal equation.

std::string SynapseGSecondTypeWCUT::synapsetype = "G-model_2nd_type_with_cut";

double SynapseGSecondTypeWCUT::init_g = 0.0;
double SynapseGSecondTypeWCUT::init_tau_se = 3;
double SynapseGSecondTypeWCUT::init_tau_si = 7;
double SynapseGSecondTypeWCUT::init_gs = 0.4;
double SynapseGSecondTypeWCUT::init_Ee = 0;
double SynapseGSecondTypeWCUT::init_Ei = -80;

double* SynapseGSecondTypeWCUT::exportData(){
    double* arr = new double[numEssentialVariables()+1];
    arr[0] = numEssentialVariables()+1;
    arr[1] = presynaptic;
    arr[2] = postsynaptic;
    arr[3] = weight;
    arr[4] = delay;
    arr[5] = g;
    arr[6] = tau_s;
    arr[7] = gs;
    arr[8] = E;
    return arr;
}

int SynapseGSecondTypeWCUT::importData(double *arr){
    presynaptic =   arr[0];
    postsynaptic =  arr[1];
    weight =        arr[2];
    delay =         arr[3];
    g =             arr[4];
    tau_s =         arr[5];
    gs =            arr[6];
    E =             arr[7];
    if(tau_s < 1e-5)
        tau_s = 1e-5;

    last_spiked = -100;
    last_spiked_post = -100;
    return 0;
}

int SynapseGSecondTypeWCUT::numEssentialVariables(){
    return 8;
}

int SynapseGSecondTypeWCUT::initSynapsesLocal(){
    using namespace vf_file;
    std::string buf30 = loadFileToString("./init/SynapseG.ini");

    init_g = getParameterIni("INIT_g", buf30);
    init_tau_se = getParameterIni("tau_se", buf30);
    init_tau_si = getParameterIni("tau_si", buf30);
    init_gs = getParameterIni("g_s", buf30);
    init_Ee = getParameterIni("Ee", buf30);
    init_Ei = getParameterIni("Ei", buf30);

    FILE* fid = fopen (getFilenameFromIni(DATAFILES, FILE_EXP_INIT_SYNAPSES).c_str(), "w");
    if(!fid){
        exit(15);
    }
    fprintf(fid, "g = %.4f;\n", init_g);
    fprintf(fid, "tau_se = %.4f;\n", init_tau_se);
    fprintf(fid, "tau_si = %.4f;\n", init_tau_si);
    fprintf(fid, "g_s = %.4f;\n", init_gs);
    fprintf(fid, "Ee = %.4f;\n", init_Ee);
    fprintf(fid, "Ei = %.4f;\n", init_Ei);
    fclose(fid);
    
    return 0;
}

int SynapseGSecondTypeWCUT::initSynapses(){
    return initSynapsesLocal();
}


std::string SynapseGSecondTypeWCUT::getName(){
    return synapsetype;
}

void SynapseGSecondTypeWCUT::setData(int pre, int pos, int preex, int posex, \
                                     double dt){
    presynaptic = pre;
    postsynaptic = pos;
    weight = vf_distributions::uniform(0, 1);
    last_spiked = -100;
    last_spiked_post = -100;
    gs = init_gs;
    g = init_g;

    if(preex){
        tau_s = init_tau_se;
        E = init_Ee;
    } else {
        tau_s = init_tau_si;
        E = init_Ei;
    }

    now_spikes = 0;
    max_spikes = (tau_s*6.6383 + 1)/dt; //accuracy 0.01
    spikerow = Malloc(max_spikes, double);//new double(max_spikes);//
    for(int i=0; i<max_spikes; i++)
        spikerow[i] = 0;

}

double SynapseGSecondTypeWCUT::alphaFunction(double x){
    return x/tau_s * exp(-x/tau_s);
}

double SynapseGSecondTypeWCUT::evolve(double dt, double timen, double Vpre, \
                                      double Vpost){
    if(timen - last_spiked_post<dt/2){
        now_spikes = 0;
    }

    if(timen - last_spiked < dt/2){
        spikerow[now_spikes%max_spikes] = timen;
        now_spikes++;
    }

    g = 0.0;
    for(int i=0; i<now_spikes && i<max_spikes; i++)
        g += alphaFunction(timen - spikerow[i]);

    out_current = g * (E - Vpost) * gs;

    return moveDeliveries();
}
