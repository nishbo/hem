#include "synapse.h"

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
    exit(11);
    return 0;
}

double* Synapse::exportData(){
    return NULL;
}

int Synapse::initSynapses(){
    return 1;
}

int Synapse::numEssentialVariables(){
    return 0;
}

/// Static synapse

std::string SynapseStatic::synapsetype = "Static_synapse";
double SynapseStatic::init_weight = 0.5;

int SynapseStatic::initSynapsesLocal(){
    FILE* fid = fopen("./init/SynapseStatic.ini", "r");
    if(!fid){
        exit(29);
    }
    float buf01;
    fscanf(fid, "Weight = %f\n", &buf01);
    init_weight = buf01;

    fclose(fid);
    return 0;
}

double* SynapseStatic::exportData(){
    if(!working)
        return new double;
    double* arr = new double[4];
    arr[0] = 4;
    arr[1] = presynaptic;
    arr[2] = postsynaptic;
    arr[3] = weight;
    return arr;
}


int SynapseStatic::initSynapses(){
    return initSynapsesLocal();
}

std::string SynapseStatic::getName(){
    return synapsetype;
}

SynapseStatic::SynapseStatic (){
    working = 0;
}

double SynapseStatic::evolve(double dt, double time, double Vpre, double Vpost){
    out_current = weight;

    return moveDeliveries();
}

void SynapseStatic::setData (int pre, int pos, int preex, int posex, double dt){
  // sets data like pre- and postsyn neuron, etc.
  working = 1;
  presynaptic = pre;
  postsynaptic = pos;
  weight = init_weight;
  last_spiked = -100;
  last_spiked_post = -100;
}

/// Tsodyks-Markram model

std::string SynapseTsodyksMarkram::synapsetype = \
        "Tsodyks-Markram synapse (euler)";

std::string SynapseTsodyksMarkram::getName(){
    return synapsetype;
}

double SynapseTsodyksMarkram::xav = 0;
double SynapseTsodyksMarkram::yav = 0;
double SynapseTsodyksMarkram::zav = 0;
double SynapseTsodyksMarkram::uav = 0;
int SynapseTsodyksMarkram::number_of_synapses = 0;

double* SynapseTsodyksMarkram::getInnerData(){
    innerDataArr[0] = 5;
    innerDataArr[1] = xav;
    innerDataArr[2] = yav;
    innerDataArr[3] = zav;
    innerDataArr[4] = uav;
    xav = 0;
    yav = 0;
    zav = 0;
    uav = 0;
    return innerDataArr;
}

SynapseTsodyksMarkram::SynapseTsodyksMarkram (){
    working = 0;
}
void SynapseTsodyksMarkram::setData(int pre, int pos, \
                                    int preex, int posex, double dt){
  tau_one = 3;

  working = 1;
  presynaptic = pre;
  postsynaptic = pos;
  weight = VFDistributions::drand();
  last_spiked = -100;
  last_spiked_post = -100;

  if(preex && posex){
    //exc exc
    A = 38;
    U = 0.5;
    tau_rec = 800;
    tau_facil = 0;
    exc = 1;
    A = VFDistributions::normal(A, A/2);
  } else if(preex){
    //exc inh
    A = 54;
    U = 0.5;
    tau_rec = 800;
    tau_facil = 0;
    exc = 1;
    A = VFDistributions::normal(A, A/2);
  } else if(posex){
    // inh exc
    A = -72;
    U = 0.04;
    tau_rec = 100;
    tau_facil = 1000;
    exc = 0;
    A = VFDistributions::normal(A, -A/2);
  } else {
    // inh inh
    A = -72;
    U = 0.04;
    tau_rec = 100;
    tau_facil = 1000;
    exc = 0;
    A = VFDistributions::normal(A, -A/2);
  }
  U = VFDistributions::normal(U, U/2);
  tau_rec = dt + abs(VFDistributions::normal(tau_rec, tau_rec/2));
  tau_facil = dt + abs(VFDistributions::normal(tau_facil, tau_facil/2));

  x = 0.28;
  y = 0.01;
  z = 0.71;
  u = U;

  number_of_synapses++;
  if(number_of_synapses==1)
      innerDataArr = new double[5];
}

double SynapseTsodyksMarkram::evolve(double dt, double time, \
                                     double Vpre, double Vpost){
  xo = x;
  yo = y;
  zo = z;
  uo = u;
  dd = VFDiscrete::diracDelta(time-last_spiked, dt);

  x += dt*(zo/tau_rec - uo * xo * dd);
  y += dt*(-yo/tau_one + uo * xo * dd);
  z += dt*(-zo/tau_rec + yo/tau_one);
  if(!exc)
    u += dt*(-uo/tau_facil + U * (1 - uo) * dd);

  //mistakes:
  if(x >= 1) x = 0.9999999; if(x<=0) x = 0.0000001;
  if(y >= 1) y = 0.9999999; if(y<=0) y = 0.0000001;
  if(z >= 1) z = 0.9999999; if(z<=0) z = 0.0000001;
  if(u >= 1) u = 0.9999999; if(u<=0) u = 0.0000001;

  out_current = A * y;
  weight = y;

  xav += x/number_of_synapses;
  yav += y/number_of_synapses;
  zav += z/number_of_synapses;
  uav += u/number_of_synapses;

  return moveDeliveries();
}

/// Tsodyks-Markram model rk4

std::string SynapseTsodyksMarkramRK::synapsetype = \
        "Tsodyks-Markram_synapse_(rk4)";

std::string SynapseTsodyksMarkramRK::getName(){
    return synapsetype;
}

double SynapseTsodyksMarkramRK::xav = 0;
double SynapseTsodyksMarkramRK::yav = 0;
double SynapseTsodyksMarkramRK::zav = 0;
double SynapseTsodyksMarkramRK::uav = 0;
int SynapseTsodyksMarkramRK::number_of_synapses = 0;

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

double* SynapseTsodyksMarkramRK::getInnerData(){
    innerDataArr[0] = 5;
    innerDataArr[1] = xav;
    innerDataArr[2] = yav;
    innerDataArr[3] = zav;
    innerDataArr[4] = uav;
    xav = 0;
    yav = 0;
    zav = 0;
    uav = 0;
    return innerDataArr;
}

SynapseTsodyksMarkramRK::SynapseTsodyksMarkramRK (){
    working = 0;
}
void SynapseTsodyksMarkramRK::setData(int pre, int pos, \
                                    int preex, int posex, double dt){
    tau_one = init_tau_one;

    working = 1;
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
    A = VFDistributions::normal(A, abs(A/2));
    U = VFDistributions::normal(U, U/2);
    tau_rec = dt + abs(VFDistributions::normal(tau_rec, tau_rec/2));
    tau_facil = dt + abs(VFDistributions::normal(tau_facil, tau_facil/2));

    x = init_x;
    y = init_y;
    z = init_z;
    u = U;

    number_of_synapses++;
    if(number_of_synapses==1)
      innerDataArr = new double[5];
}

double SynapseTsodyksMarkramRK::Xr(double x1, double y1, double z1, double u1, \
                                  double t1, double dt){
    return z1/tau_rec - VFDiscrete::diracDelta(t1 - last_spiked, dt) * u1 * x1;
}

double SynapseTsodyksMarkramRK::Yr(double x1, double y1, double z1, double u1, \
                                  double t1, double dt){
    return -y1/tau_one + VFDiscrete::diracDelta(t1 - last_spiked, dt) * u1 * x1;
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
        return -u1/tau_facil + VFDiscrete::diracDelta(t1 - last_spiked, dt) \
                * U * (1 - u1);
}

double SynapseTsodyksMarkramRK::evolve(double dt, double time, \
                                       double Vpre, double Vpost){


    k_1_x = dt * Xr(x, y, z, u, time, dt);
    k_1_y = dt * Yr(x, y, z, u, time, dt);
    k_1_z = dt * Zr(x, y, z, u, time, dt);
    k_1_u = dt * Ur(x, y, z, u, time, dt);

    k_2_x = dt * Xr(x + k_1_x/2.0, y + k_1_y/2.0, z + k_1_z/2.0, \
                    u + k_1_u/2.0, time + dt/2.0, dt);
    k_2_y = dt * Yr(x + k_1_x/2.0, y + k_1_y/2.0, z + k_1_z/2.0, \
                    u + k_1_u/2.0, time + dt/2.0, dt);
    k_2_z = dt * Zr(x + k_1_x/2.0, y + k_1_y/2.0, z + k_1_z/2.0, \
                    u + k_1_u/2.0, time + dt/2.0, dt);
    k_2_u = dt * Ur(x + k_1_x/2.0, y + k_1_y/2.0, z + k_1_z/2.0, \
                    u + k_1_u/2.0, time + dt/2.0, dt);

    k_3_x = dt * Xr(x + k_2_x/2.0, y + k_2_y/2.0, z + k_2_z/2.0, \
                    u + k_2_u/2.0, time + dt/2.0, dt);
    k_3_y = dt * Yr(x + k_2_x/2.0, y + k_2_y/2.0, z + k_2_z/2.0, \
                    u + k_2_u/2.0, time + dt/2.0, dt);
    k_3_z = dt * Zr(x + k_2_x/2.0, y + k_2_y/2.0, z + k_2_z/2.0, \
                    u + k_2_u/2.0, time + dt/2.0, dt);
    k_3_u = dt * Ur(x + k_2_x/2.0, y + k_2_y/2.0, z + k_2_z/2.0, \
                    u + k_2_u/2.0, time + dt/2.0, dt);

    x += (k_1_x + 2.0*k_2_x + 2.0*k_3_x + k_4_x)/6.0;
    y += (k_1_y + 2.0*k_2_y + 2.0*k_3_y + k_4_y)/6.0;
    z += (k_1_z + 2.0*k_2_z + 2.0*k_3_z + k_4_z)/6.0;
    u += (k_1_u + 2.0*k_2_u + 2.0*k_3_u + k_4_u)/6.0;

    //mistakes:
    if(x >= 1) x = 0.9999999; if(x<=0) x = 0.0000001;
    if(y >= 1) y = 0.9999999; if(y<=0) y = 0.0000001;
    if(z >= 1) z = 0.9999999; if(z<=0) z = 0.0000001;
    if(u >= 1) u = 0.9999999; if(u<=0) u = 0.0000001;

    out_current = A * y;

    xav += x/number_of_synapses;
    yav += y/number_of_synapses;
    zav += z/number_of_synapses;
    uav += u/number_of_synapses;

    return moveDeliveries();
}

double* SynapseTsodyksMarkramRK::exportData(){
    if(!working)
        return new double;
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
    working = 1;
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

    if(tau_rec < 1e-3)
        tau_rec = 1e-3;
    if(tau_facil < 1e-3)
        tau_facil = 1e-3;
    return 0;
}

int SynapseTsodyksMarkramRK::numEssentialVariables(){
    return 14;
}

int SynapseTsodyksMarkramRK::initSynapsesLocal(){
    FILE* fid = fopen("./init/SynapseTsodyksMarkram.ini", "r");
    if(!fid){
        exit(29);
    }
    float buf01;
    fscanf(fid, "tau_one = %f\n", &buf01);
    init_tau_one = buf01;
    fscanf(fid, "x = %f\n", &buf01);
    init_x = buf01;
    fscanf(fid, "y = %f\n", &buf01);
    init_y = buf01;
    fscanf(fid, "z = %f\n\n", &buf01);
    init_z = buf01;

    fscanf(fid, "Aee = %f\n", &buf01);
    init_Aee = buf01;
    fscanf(fid, "Uee = %f\n", &buf01);
    init_Uee = buf01;
    fscanf(fid, "tau_recee = %f\n", &buf01);
    init_tau_recee = buf01;
    fscanf(fid, "tau_facilee = %f\n\n", &buf01);
    init_tau_facilee = buf01;

    fscanf(fid, "Aei = %f\n", &buf01);
    init_Aei = buf01;
    fscanf(fid, "Uei = %f\n", &buf01);
    init_Uei = buf01;
    fscanf(fid, "tau_recei = %f\n", &buf01);
    init_tau_recei = buf01;
    fscanf(fid, "tau_facilei = %f\n\n", &buf01);
    init_tau_facilei = buf01;

    fscanf(fid, "Aie = %f\n", &buf01);
    init_Aie = buf01;
    fscanf(fid, "Uie = %f\n", &buf01);
    init_Uie = buf01;
    fscanf(fid, "tau_recie = %f\n", &buf01);
    init_tau_recie = buf01;
    fscanf(fid, "tau_facilie = %f\n\n", &buf01);
    init_tau_facilie = buf01;

    fscanf(fid, "Aii = %f\n", &buf01);
    init_Aii = buf01;
    fscanf(fid, "Uii = %f\n", &buf01);
    init_Uii = buf01;
    fscanf(fid, "tau_recii = %f\n", &buf01);
    init_tau_recii = buf01;
    fscanf(fid, "tau_facilii = %f\n", &buf01);
    init_tau_facilii = buf01;

    fclose(fid);
    return 0;
}

int SynapseTsodyksMarkramRK::initSynapses(){
    return initSynapsesLocal();
}

/// Tsodyks-Markram analytical

std::string SynapseTsodyksMarkramAn::synapsetype = \
        "Tsodyks-Markram_synapse_(an)";

std::string SynapseTsodyksMarkramAn::getName(){
    return synapsetype;
}

double SynapseTsodyksMarkramAn::init_tau_one = 3;
double SynapseTsodyksMarkramAn::init_x = 0.98;
double SynapseTsodyksMarkramAn::init_y = 0.01;
double SynapseTsodyksMarkramAn::init_z = 0.01;
double SynapseTsodyksMarkramAn::init_Aee = 38;
double SynapseTsodyksMarkramAn::init_Uee = 0.5;
double SynapseTsodyksMarkramAn::init_tau_recee = 800;
double SynapseTsodyksMarkramAn::init_tau_facilee = 0;
double SynapseTsodyksMarkramAn::init_Aei = 54;
double SynapseTsodyksMarkramAn::init_Uei = 0.5;
double SynapseTsodyksMarkramAn::init_tau_recei = 800;
double SynapseTsodyksMarkramAn::init_tau_facilei = 0;
double SynapseTsodyksMarkramAn::init_Aie = -72;
double SynapseTsodyksMarkramAn::init_Uie = 0.04;
double SynapseTsodyksMarkramAn::init_tau_recie = 100;
double SynapseTsodyksMarkramAn::init_tau_facilie = 1000;
double SynapseTsodyksMarkramAn::init_Aii = -72;
double SynapseTsodyksMarkramAn::init_Uii = 0.04;
double SynapseTsodyksMarkramAn::init_tau_recii = 100;
double SynapseTsodyksMarkramAn::init_tau_facilii = 1000;

SynapseTsodyksMarkramAn::SynapseTsodyksMarkramAn (){
    working = 0;
}
void SynapseTsodyksMarkramAn::setData(int pre, int pos, \
                                    int preex, int posex, double dt){
    tau_one = init_tau_one;

    working = 1;
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
    A = VFDistributions::normal(A, abs(A/2));
    U = VFDistributions::normal(U, U/2);
    tau_rec = dt + abs(VFDistributions::normal(tau_rec, tau_rec/2));
    tau_facil = dt + abs(VFDistributions::normal(tau_facil, tau_facil/2));

    x = init_x;
    y = init_y;
    z = init_z;
    u = U;
}

double SynapseTsodyksMarkramAn::evolve(double dt, double time, \
                                       double Vpre, double Vpost){
    buf00 = VFDiscrete::diracDelta(time - last_spiked, dt);

    if(y > 0.00001){
        eone    = exp(- (time - last_spiked) / tau_one);
        erec    = exp(- (time - last_spiked) / tau_rec);

        if(!exc){
            efacil  = exp(- (time - last_spiked) / tau_facil);
            u = usp * efacil;
        }
        y = ysp * eone;
        z = erec / (tau_rec - tau_one) * (zsp * (tau_rec - tau_one) + \
                                          ysp * tau_rec * (1 - eone/erec));
        x = (xsp + ysp + zsp) + 1/(tau_one - tau_rec)*(-ysp * tau_one * eone + \
                                   (zsp*(tau_rec - tau_one)+ysp*tau_rec)*erec);
    }
    if(buf00){
        eone    = exp(- (time - last_spiked - dt) / tau_one);
        erec    = exp(- (time - last_spiked - dt) / tau_rec);

        if(!exc){
            efacil  = exp(- (time - last_spiked - dt) / tau_facil);
            u = usp * efacil;
        }
        y = ysp * eone;
        z = erec / (tau_rec - tau_one) * (zsp * (tau_rec - tau_one) + \
                                          ysp * tau_rec * (1 - eone/erec));
        x = (xsp + ysp + zsp) + 1/(tau_one - tau_rec)*(-ysp * tau_one * eone + \
                                   (zsp*(tau_rec - tau_one)+ysp*tau_rec)*erec);

        x -= u * x;
        y += u * x;
        u += U * (1 - u);
    }


    //mistakes:
    if(x >= 1) x = 0.9999999; else if(x<=0) x = 0.0000001;
    if(y >= 1) y = 0.9999999; else if(y<=0) y = 0.0000001;
    if(z >= 1) z = 0.9999999; else if(z<=0) z = 0.0000001;
    if(u >= 1) u = 0.9999999; else if(u<=0) u = 0.0000001;

    out_current = A * y;

    return moveDeliveries();
}

double* SynapseTsodyksMarkramAn::exportData(){
    if(!working)
        return new double;
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

int SynapseTsodyksMarkramAn::importData(double *arr){
    working = 1;
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

    if(tau_rec < 1e-3)
        tau_rec = 1e-3;
    if(tau_facil < 1e-3)
        tau_facil = 1e-3;
    return 0;
}

int SynapseTsodyksMarkramAn::numEssentialVariables(){
    return 14;
}

int SynapseTsodyksMarkramAn::initSynapsesLocal(){
    FILE* fid = fopen("./init/SynapseTsodyksMarkram.ini", "r");
    if(!fid){
        exit(29);
    }
    float buf01;
    fscanf(fid, "tau_one = %f\n", &buf01);
    init_tau_one = buf01;
    fscanf(fid, "x = %f\n", &buf01);
    init_x = buf01;
    fscanf(fid, "y = %f\n", &buf01);
    init_y = buf01;
    fscanf(fid, "z = %f\n\n", &buf01);
    init_z = buf01;

    fscanf(fid, "Aee = %f\n", &buf01);
    init_Aee = buf01;
    fscanf(fid, "Uee = %f\n", &buf01);
    init_Uee = buf01;
    fscanf(fid, "tau_recee = %f\n", &buf01);
    init_tau_recee = buf01;
    fscanf(fid, "tau_facilee = %f\n\n", &buf01);
    init_tau_facilee = buf01;

    fscanf(fid, "Aei = %f\n", &buf01);
    init_Aei = buf01;
    fscanf(fid, "Uei = %f\n", &buf01);
    init_Uei = buf01;
    fscanf(fid, "tau_recei = %f\n", &buf01);
    init_tau_recei = buf01;
    fscanf(fid, "tau_facilei = %f\n\n", &buf01);
    init_tau_facilei = buf01;

    fscanf(fid, "Aie = %f\n", &buf01);
    init_Aie = buf01;
    fscanf(fid, "Uie = %f\n", &buf01);
    init_Uie = buf01;
    fscanf(fid, "tau_recie = %f\n", &buf01);
    init_tau_recie = buf01;
    fscanf(fid, "tau_facilie = %f\n\n", &buf01);
    init_tau_facilie = buf01;

    fscanf(fid, "Aii = %f\n", &buf01);
    init_Aii = buf01;
    fscanf(fid, "Uii = %f\n", &buf01);
    init_Uii = buf01;
    fscanf(fid, "tau_recii = %f\n", &buf01);
    init_tau_recii = buf01;
    fscanf(fid, "tau_facilii = %f\n", &buf01);
    init_tau_facilii = buf01;

    fclose(fid);
    return 0;
}

int SynapseTsodyksMarkramAn::initSynapses(){
    return initSynapsesLocal();
}

/// Spike-Timing Dependent Plasticity Model

std::string SynapseSTDP::synapsetype = "STDP";

std::string SynapseSTDP::getName(){
    return synapsetype;
}

SynapseSTDP::SynapseSTDP (){
    working = 0;
}

void SynapseSTDP::setData(int pre, int post, int preex, int posex, double dt){
    working = 1;
    presynaptic = pre;
    postsynaptic = post;
    last_spiked = -100;
    last_spiked_post = -100;

    // STDP rule:
    // Dt = t_post - t_pre,
    // Dw = lambda*(1-w)^mu*exp(-|Dt|/tau), if Dt > 0,
    // Dw = -lambda*alpha*w^mu*exp(-|Dt|/tau), if Dt <= 0,
    // where 0 < lambda << 1, 0 <= mu <= 1, alpha ~ 1.
    lambda = 0.1;
//    mu = 0.5;
    alpha = 1;
    tau_corr = 20; //ms

    weight = VFDistributions::drand();
}

double SynapseSTDP::evolve(double dt, double time, double Vpre, double Vpost){

    if(last_spiked>0 && last_spiked_post>0 && \
            VFDiscrete::diracDelta(last_spiked_post - time, dt) && \
            VFDiscrete::diracDelta(last_spiked - time, dt)){
        if(last_spiked_post > last_spiked){
            weight += lambda * sqrt(1-weight)*\
                    exp(-abs(last_spiked_post - last_spiked)/tau_corr);
        } else {
            weight -= lambda * alpha * sqrt(weight) * \
                    exp(-abs(last_spiked_post - last_spiked)/tau_corr);
        }
        out_current = weight;
    } else {
        out_current = 0;
    }
    if(weight < 1e-8) weight = 1e-8; //bug? shows NaN else
    if(weight > 1) weight = 1; //calc mistakes

    return moveDeliveries();
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

SynapseSTDPG::SynapseSTDPG (){
    working = 0;
}

void SynapseSTDPG::setData(int pre, int post, int preex, int posex, double dt){
    working = 1;
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

    if(last_spiked>0 && last_spiked_post>0 && \
            VFDiscrete::diracDelta(last_spiked_post - time, dt) && \
            VFDiscrete::diracDelta(last_spiked - time, dt)){
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
    g += dt*(-g/tau_s) + VFDiscrete::diracDelta(last_spiked - time, dt);
    out_current = weight * g;

    return moveDeliveries();
}

double* SynapseSTDPG::exportData(){
    if(!working)
        return new double;
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
    working = 1;
    presynaptic = arr[0];
    postsynaptic = arr[1];
    weight = arr[2];
    delay = arr[3];
    g = arr[4];
    tau_s = arr[5];
    lambda = arr[6];
    alpha = arr[7];
    tau_corr = arr[8];

    if(tau_corr < 1e-3)
        tau_corr = 1e-3;
    return 0;
}

int SynapseSTDPG::numEssentialVariables(){
    return 9;
}

int SynapseSTDPG::initSynapsesLocal(){
    FILE* fid = fopen("./init/SynapseSTDPG.ini", "r");
    if(!fid){
        exit(29);
    }
    float buf01;
    fscanf(fid, "g = %f\n", &buf01);
    init_g = buf01;
    fscanf(fid, "tau_s = %f\n", &buf01);
    init_tau_s = buf01;
    fscanf(fid, "lambda = %f\n", &buf01);
    init_lambda = buf01;
    fscanf(fid, "alpha = %f\n", &buf01);
    init_alpha = buf01;
    fscanf(fid, "tau_corr = %f\n", &buf01);
    init_tau_corr = buf01;
    fscanf(fid, "weight = %f\n", &buf01);
    init_weight = buf01;

    fclose(fid);
    return 0;
}

int SynapseSTDPG::initSynapses(){
    return initSynapsesLocal();
}

/// STDP AND T-M

std::string SynapseTMSTDP::synapsetype = "Tsodyks-Markram_and_STDP_model";

double SynapseTMSTDP::init_lambda = 0.1;
double SynapseTMSTDP::init_alpha = 1;
double SynapseTMSTDP::init_tau_corr = 20;
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

std::string SynapseTMSTDP::getName(){
    return synapsetype;
}

SynapseTMSTDP::SynapseTMSTDP (){
  working = 0;
}
void SynapseTMSTDP::setData(int pre, int pos, \
                                    int preex, int posex, double dt){
    tau_one = init_tau_one;
    t_start = init_t_start;

    working = 1;
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
    A = VFDistributions::normal(A, abs(A/2));
    U = VFDistributions::normal(U, U/2);
    tau_rec = dt + abs(VFDistributions::normal(tau_rec, tau_rec/2));
    tau_facil = dt + abs(VFDistributions::normal(tau_facil, tau_facil/2));

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
        weight = VFDistributions::uniform(0, 1);
        break;
    case 2:
        weight = VFDistributions::normal(init_weight, init_weight*0.1, 0, 1);
        break;
    default:
        weight = init_weight;
    }

}

double SynapseTMSTDP::Xr(double x1, double y1, double z1, double u1, \
                                  double t1, double dt){
    return z1/tau_rec - VFDiscrete::diracDelta(t1 - last_spiked, dt) * u1 * x1;
}

double SynapseTMSTDP::Yr(double x1, double y1, double z1, double u1, \
                                  double t1, double dt){
    return -y1/tau_one + VFDiscrete::diracDelta(t1 - last_spiked, dt) * u1 * x1;
}

double SynapseTMSTDP::Zr(double x1, double y1, double z1, double u1, \
                                  double t1, double dt){
    return -z1/tau_rec + y1/tau_one;
}

double SynapseTMSTDP::Ur(double x1, double y1, double z1, double u1, \
                                  double t1, double dt){
    if(exc)
        return U;
    else
        return -u1/tau_facil + VFDiscrete::diracDelta(t1 - last_spiked, dt) \
                * U * (1 - u1);
}

double SynapseTMSTDP::evolve(double dt, double time, \
                                       double Vpre, double Vpost){

    k_1_x = dt * Xr(x, y, z, u, time, dt);
    k_1_y = dt * Yr(x, y, z, u, time, dt);
    k_1_z = dt * Zr(x, y, z, u, time, dt);
    k_1_u = dt * Ur(x, y, z, u, time, dt);

    k_2_x = dt * Xr(x + k_1_x/2.0, y + k_1_y/2.0, z + k_1_z/2.0, \
                    u + k_1_u/2.0, time + dt/2.0, dt);
    k_2_y = dt * Yr(x + k_1_x/2.0, y + k_1_y/2.0, z + k_1_z/2.0, \
                    u + k_1_u/2.0, time + dt/2.0, dt);
    k_2_z = dt * Zr(x + k_1_x/2.0, y + k_1_y/2.0, z + k_1_z/2.0, \
                    u + k_1_u/2.0, time + dt/2.0, dt);
    k_2_u = dt * Ur(x + k_1_x/2.0, y + k_1_y/2.0, z + k_1_z/2.0, \
                    u + k_1_u/2.0, time + dt/2.0, dt);

    k_3_x = dt * Xr(x + k_2_x/2.0, y + k_2_y/2.0, z + k_2_z/2.0, \
                    u + k_2_u/2.0, time + dt/2.0, dt);
    k_3_y = dt * Yr(x + k_2_x/2.0, y + k_2_y/2.0, z + k_2_z/2.0, \
                    u + k_2_u/2.0, time + dt/2.0, dt);
    k_3_z = dt * Zr(x + k_2_x/2.0, y + k_2_y/2.0, z + k_2_z/2.0, \
                    u + k_2_u/2.0, time + dt/2.0, dt);
    k_3_u = dt * Ur(x + k_2_x/2.0, y + k_2_y/2.0, z + k_2_z/2.0, \
                    u + k_2_u/2.0, time + dt/2.0, dt);

    k_4_x = dt * Xr(x + k_3_x, y + k_3_y, z + k_3_z, u + k_3_u, time + dt, dt);
    k_4_y = dt * Yr(x + k_3_x, y + k_3_y, z + k_3_z, u + k_3_u, time + dt, dt);
    k_4_z = dt * Zr(x + k_3_x, y + k_3_y, z + k_3_z, u + k_3_u, time + dt, dt);
    k_4_u = dt * Ur(x + k_3_x, y + k_3_y, z + k_3_z, u + k_3_u, time + dt, dt);

    x += (k_1_x + 2.0*k_2_x + 2.0*k_3_x + k_4_x)/6.0;
    y += (k_1_y + 2.0*k_2_y + 2.0*k_3_y + k_4_y)/6.0;
    z += (k_1_z + 2.0*k_2_z + 2.0*k_3_z + k_4_z)/6.0;
    u += (k_1_u + 2.0*k_2_u + 2.0*k_3_u + k_4_u)/6.0;

    if(x >= 1) x = 0.9999999; if(x<=0) x = 0.0000001;
    if(y >= 1) y = 0.9999999; if(y<=0) y = 0.0000001;
    if(z >= 1) z = 0.9999999; if(z<=0) z = 0.0000001;
    if(u >= 1) u = 0.9999999; if(u<=0) u = 0.0000001;

    if(time>t_start){
        if(last_spiked>0 && last_spiked_post>0 && \
                VFDiscrete::diracDelta(last_spiked_post - time, dt) && \
                VFDiscrete::diracDelta(last_spiked - time, dt)){
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

    out_current = A * y * weight * 2;

    return moveDeliveries();
}

double* SynapseTMSTDP::exportData(){
    if(!working)
        return new double;
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
    working = 1;
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

    if(tau_rec < 1e-3)
        tau_rec = 1e-3;
    if(tau_facil < 1e-3)
        tau_facil = 1e-3;

    lambda = init_lambda;
    alpha = init_alpha;
    tau_corr = init_tau_corr;
    t_start = init_t_start;
    return 0;
}

int SynapseTMSTDP::numEssentialVariables(){
    return 14;
}

int SynapseTMSTDP::initSynapsesLocal(){
    FILE* fid = fopen("./init/SynapseTMSTDP.ini", "r");
    if(!fid){
        exit(29);
    }
    float buf01;
    fscanf(fid, "t_start = %f\n", &buf01);
    init_t_start = buf01;
    fscanf(fid, "type_of_weight = %d\n", &init_type_of_weight);
    fscanf(fid, "lambda = %f\n", &buf01);
    init_lambda = buf01;
    fscanf(fid, "alpha = %f\n", &buf01);
    init_alpha = buf01;
    fscanf(fid, "tau_corr = %f\n", &buf01);
    init_tau_corr = buf01;
    fscanf(fid, "weight = %f\n", &buf01);
    init_weight = buf01;
    fscanf(fid, "tau_one = %f\n", &buf01);
    init_tau_one = buf01;
    fscanf(fid, "x = %f\n", &buf01);
    init_x = buf01;
    fscanf(fid, "y = %f\n", &buf01);
    init_y = buf01;
    fscanf(fid, "z = %f\n\n", &buf01);
    init_z = buf01;

    fscanf(fid, "Aee = %f\n", &buf01);
    init_Aee = buf01;
    fscanf(fid, "Uee = %f\n", &buf01);
    init_Uee = buf01;
    fscanf(fid, "tau_recee = %f\n", &buf01);
    init_tau_recee = buf01;
    fscanf(fid, "tau_facilee = %f\n\n", &buf01);
    init_tau_facilee = buf01;

    fscanf(fid, "Aei = %f\n", &buf01);
    init_Aei = buf01;
    fscanf(fid, "Uei = %f\n", &buf01);
    init_Uei = buf01;
    fscanf(fid, "tau_recei = %f\n", &buf01);
    init_tau_recei = buf01;
    fscanf(fid, "tau_facilei = %f\n\n", &buf01);
    init_tau_facilei = buf01;

    fscanf(fid, "Aie = %f\n", &buf01);
    init_Aie = buf01;
    fscanf(fid, "Uie = %f\n", &buf01);
    init_Uie = buf01;
    fscanf(fid, "tau_recie = %f\n", &buf01);
    init_tau_recie = buf01;
    fscanf(fid, "tau_facilie = %f\n\n", &buf01);
    init_tau_facilie = buf01;

    fscanf(fid, "Aii = %f\n", &buf01);
    init_Aii = buf01;
    fscanf(fid, "Uii = %f\n", &buf01);
    init_Uii = buf01;
    fscanf(fid, "tau_recii = %f\n", &buf01);
    init_tau_recii = buf01;
    fscanf(fid, "tau_facilii = %f\n", &buf01);
    init_tau_facilii = buf01;

    fclose(fid);
    return 0;
}

int SynapseTMSTDP::initSynapses(){
    return initSynapsesLocal();
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
    if(!working)
        return new double;
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
    working = 1;
    presynaptic =   arr[0];
    postsynaptic =  arr[1];
    weight =        arr[2];
    delay =         arr[3];
    g =             arr[4];
    tau_s =         arr[5];
    gs =            arr[6];
    E =             arr[7];
    if(tau_s < 1e-3)
        tau_s = 1e-3;
    return 0;
}

int SynapseGFirstType::numEssentialVariables(){
    return 8;
}

int SynapseGFirstType::initSynapsesLocal(){
    FILE* fid = fopen("./init/SynapseG.ini", "r");
    if(!fid){
        exit(29);
    }
    float buf01;
    fscanf(fid, "g = %f\n", &buf01);
    init_g = buf01;
    fscanf(fid, "tau_se = %f\n", &buf01);
    init_tau_se = buf01;
    fscanf(fid, "tau_si = %f\n", &buf01);
    init_tau_si = buf01;
    fscanf(fid, "g_s = %f\n", &buf01);
    init_gs = buf01;
    fscanf(fid, "Ee = %f\n", &buf01);
    init_Ee = buf01;
    fscanf(fid, "Ei = %f\n", &buf01);
    init_Ei = buf01;

    fclose(fid);
    return 0;
}

int SynapseGFirstType::initSynapses(){
    return initSynapsesLocal();
}

std::string SynapseGFirstType::getName(){
    return synapsetype;
}

SynapseGFirstType::SynapseGFirstType(){
    working = 0;
}
void SynapseGFirstType::setData(int pre, int pos, int preex, int posex, \
                                double dt){
    working = 1;
    presynaptic = pre;
    postsynaptic = pos;
    weight = VFDistributions::drand();
    last_spiked = -100;
    last_spiked_post = -100;
    gs = init_gs;

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
    g += dt * (-g / tau_s) + gs * VFDiscrete::diracDelta(timen-last_spiked, dt);
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
    if(!working)
        return new double;
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
    working = 1;
    presynaptic =   arr[0];
    postsynaptic =  arr[1];
    weight =        arr[2];
    delay =         arr[3];
    g =             arr[4];
    tau_s =         arr[5];
    gs =            arr[6];
    E =             arr[7];
    if(tau_s < 1e-3)
        tau_s = 1e-3;
    return 0;
}

int SynapseGFirstTypeWCUT::numEssentialVariables(){
    return 8;
}

int SynapseGFirstTypeWCUT::initSynapsesLocal(){
    FILE* fid = fopen("./init/SynapseG.ini", "r");
    if(!fid){
        exit(29);
    }
    float buf01;
    fscanf(fid, "g = %f\n", &buf01);
    init_g = buf01;
    fscanf(fid, "tau_se = %f\n", &buf01);
    init_tau_se = buf01;
    fscanf(fid, "tau_si = %f\n", &buf01);
    init_tau_si = buf01;
    fscanf(fid, "g_s = %f\n", &buf01);
    init_gs = buf01;
    fscanf(fid, "Ee = %f\n", &buf01);
    init_Ee = buf01;
    fscanf(fid, "Ei = %f\n", &buf01);
    init_Ei = buf01;

    fclose(fid);
    return 0;
}

int SynapseGFirstTypeWCUT::initSynapses(){
    return initSynapsesLocal();
}

std::string SynapseGFirstTypeWCUT::getName(){
    return synapsetype;
}

SynapseGFirstTypeWCUT::SynapseGFirstTypeWCUT(){
    working = 0;
}

void SynapseGFirstTypeWCUT::setData(int pre, int pos, int preex, int posex, \
                                    double dt){
    working = 1;
    presynaptic = pre;
    postsynaptic = pos;
    weight = VFDistributions::drand();
    last_spiked = -100;
    last_spiked_post = -100;
    gs = init_gs;

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
    g += dt * (-g / tau_s) + gs * VFDiscrete::diracDelta(timen-last_spiked, dt);
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
    if(!working)
        return new double;
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
    working = 1;
    presynaptic =   arr[0];
    postsynaptic =  arr[1];
    weight =        arr[2];
    delay =         arr[3];
    g =             arr[4];
    tau_s =         arr[5];
    gs =            arr[6];
    E =             arr[7];
    if(tau_s < 1e-3)
        tau_s = 1e-3;
    return 0;
}

int SynapseGSecondType::numEssentialVariables(){
    return 8;
}

int SynapseGSecondType::initSynapsesLocal(){
    FILE* fid = fopen("./init/SynapseG.ini", "r");
    if(!fid){
        exit(29);
    }
    float buf01;
    fscanf(fid, "g = %f\n", &buf01);
    init_g = buf01;
    fscanf(fid, "tau_se = %f\n", &buf01);
    init_tau_se = buf01;
    fscanf(fid, "tau_si = %f\n", &buf01);
    init_tau_si = buf01;
    fscanf(fid, "g_s = %f\n", &buf01);
    init_gs = buf01;
    fscanf(fid, "Ee = %f\n", &buf01);
    init_Ee = buf01;
    fscanf(fid, "Ei = %f\n", &buf01);
    init_Ei = buf01;

    fclose(fid);
    return 0;
}

int SynapseGSecondType::initSynapses(){
    return initSynapsesLocal();
}


std::string SynapseGSecondType::getName(){
    return synapsetype;
}

SynapseGSecondType::SynapseGSecondType(){
    working = 0;
}

void SynapseGSecondType::setData(int pre, int pos, int preex, int posex, \
                                 double dt){
    working = 1;
    presynaptic = pre;
    postsynaptic = pos;
    weight = VFDistributions::drand();
    last_spiked = -100;
    last_spiked_post = -100;
    gs = init_gs;

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
    if(!working)
        return new double;
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
    working = 1;
    presynaptic =   arr[0];
    postsynaptic =  arr[1];
    weight =        arr[2];
    delay =         arr[3];
    g =             arr[4];
    tau_s =         arr[5];
    gs =            arr[6];
    E =             arr[7];
    if(tau_s < 1e-3)
        tau_s = 1e-3;
    return 0;
}

int SynapseGSecondTypeWCUT::numEssentialVariables(){
    return 8;
}

int SynapseGSecondTypeWCUT::initSynapsesLocal(){
    FILE* fid = fopen("./init/SynapseG.ini", "r");
    if(!fid){
        exit(29);
    }
    float buf01;
    fscanf(fid, "g = %f\n", &buf01);
    init_g = buf01;
    fscanf(fid, "tau_se = %f\n", &buf01);
    init_tau_se = buf01;
    fscanf(fid, "tau_si = %f\n", &buf01);
    init_tau_si = buf01;
    fscanf(fid, "g_s = %f\n", &buf01);
    init_gs = buf01;
    fscanf(fid, "Ee = %f\n", &buf01);
    init_Ee = buf01;
    fscanf(fid, "Ei = %f\n", &buf01);
    init_Ei = buf01;

    fclose(fid);
    return 0;
}

int SynapseGSecondTypeWCUT::initSynapses(){
    return initSynapsesLocal();
}


std::string SynapseGSecondTypeWCUT::getName(){
    return synapsetype;
}

SynapseGSecondTypeWCUT::SynapseGSecondTypeWCUT(){
    working = 0;
}

void SynapseGSecondTypeWCUT::setData(int pre, int pos, int preex, int posex, \
                                     double dt){
    working = 1;
    presynaptic = pre;
    postsynaptic = pos;
    weight = VFDistributions::drand();
    last_spiked = -100;
    last_spiked_post = -100;
    gs = init_gs;

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

/// Prototype of a synapse
// Use this as a prototype. Do not change, copy, then alter.
std::string SynapsePrototype::synapsetype = "Namespace";

std::string SynapsePrototype::getName(){
    return synapsetype;
}

SynapsePrototype::SynapsePrototype (){
    // Do not save anything here.
    working = 0;
}

double SynapsePrototype::evolve(double dt, double time, \
                                double Vpre, double Vpost){
    // Change the state of synapse. Set to out_current current you want
    // to send to neuron. In the end MUST be return moveDeliveries();

    // For example:
    out_current = weight;

    return moveDeliveries();
}

void SynapsePrototype::setData (int pre, int pos, int preex, int posex, \
                                double dt){
  working = 1;
  presynaptic = pre;
  postsynaptic = pos;
  last_spiked = -100;
  last_spiked_post = -100;

  // Set all necessary parameters:
  weight = VFDistributions::drand();
}

