#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include <ctime>
#include <cstdio>
#include <iostream>

#include "cad.h"
#include "vfdistributions.h"
#include "vfdiscrete.h"
#include "vffile.h"

class Topology{
public:
    // static int setTopology(const int N, double *x, double *y, int *con, double *delays);
    static int setTopology(const int N, int *Mfull, double *x, double *y, int **sout, double **delays);
private:
	static std::string fle;
	static int type, delay_type;
	static double r_p, smw_beta;
	static int smw_local;
	static double border, velocity, max_delay, min_delay;

	static int setCoordinates(const int N, double *x, double *y);

	// static int randomTopology(const int N, int *con);
	// static int smallWorldTopology(const int N, int *con);

	static int randomTopology(const int N, int *Mfull, int **sout);
	static int smallWorldTopology(const int N, int *Mfull, int **sout);

	// static int setDelaysdt(const int N, const int *con, double *delays);
	// static int setDelaysRandom(const int N, const int *con, double *delays);
	// static int setDelaysCoord(const int N, const int *con, const double *x, \
	// 						  const double *y, double *delays);

	static int setDelaysdt(const int N, const int Mfull, double *delays);
	static int setDelaysRandom(const int N, const int Mfull, double *delays);
	static int setDelaysCoord(const int N, int const* const* sout, const double *x, \
                            	const double *y, double *delays);

    Topology();
};

#endif // TOPOLOGY_H
