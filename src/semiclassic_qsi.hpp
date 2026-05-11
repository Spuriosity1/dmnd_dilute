#pragma once
#include "struct/cell_geometry.hpp"
#include "struct/chain.hpp"
#include "struct/preset_cellspecs.hpp"
#include <complex>
#include <vector>
#include <random>


using namespace CellGeometry;

typedef std::complex<double> cplx;

inline double abs2(const std::complex<double> c){
    return c.real()*c.real() + c.imag()*c.imag();
}

struct Tetra : public Cell<0> {
//    double Q;
};


struct PyroSite : public Cell<1> {
    double state[3]; 
    double Sz() const { 
        return state[2];
    }
    std::complex<double> xy() const {
        return *reinterpret_cast<const cplx*>(state);
    }
};

struct Plaquette : public Cell<2> {
    double g =1.0;
    double g_prime = 0;
};


typedef PeriodicPlaqLattice<Tetra, PyroSite, Plaquette> sc_QSI;

/*
// An object responsible for precomputing which geometric objects are capable 
// of receiving what kind of update
struct MC_sampler{
    MC_sampler(const sc_QSI& qsi){
    }

    std::vector<const Tetra*> full_tetras;
    std::vector<const Plaquette*> full_plaquettes;

    private:
    // Determines which tetrahedra can be gauges WLOG
    void calc_tetras(const sc_QSI& qsi) {
        for (const Tetra& t : qsi.points){
            // checks for all 4-membered tetrahedra
            if (t.coboundary.size() == 4) {
                full_tetras.push_back(&t);
            }
        }

        for (const Plaquette& p : qsi.plaqs){
            if (p.boundary.size() == 6){
                full_plaquettes.push_back(&p);
            }
        }
    }
};
*/

// MONTE CARLO ROUTINES
namespace MC {
    // Rotates all spins of sc_QSI in sucha  way as to not break emergent gauge symmetry
    void apply_tetra_gauge(sc_QSI& lat, std::mt19937& rng); 

    // Attempts to rotate each spin
	unsigned apply_angle(sc_QSI& lat, double beta, std::mt19937& rng);

    unsigned apply_ring_Ising(sc_QSI& lat, double beta, std::mt19937&rng);
 
	unsigned apply_step(sc_QSI& lat, double beta, std::mt19937& rng);
}; // end namespace MC

