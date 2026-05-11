#include "semiclassic_qsi.hpp"
#include "random.hpp"
#include <algorithm>
#include <chrono>
#include <complex>
#include <random>




	// ENERGY AND FIELD IMPLEMENTATION
	//
	
    inline cplx ring(const Plaquette* plaq) {
        double re_ring=0;
        double im_ring=0;
        for ( auto& [site, m] : plaq->boundary ){
            auto dp = static_cast<PyroSite*>(site);
            auto exp_ia = dp->xy();
            assert(m*m == 1);
			// performs ring -> ring * e^+- ia
            re_ring = re_ring*exp_ia.real() - im_ring * m* exp_ia.imag();
            im_ring = im_ring*exp_ia.real() + re_ring * m*exp_ia.imag();
        }
        return cplx(re_ring, im_ring);
    }


    inline cplx partial_ring(const Plaquette* plaq, const PyroSite* self_site) {
        double re_ring=0;
        double im_ring=0;

        for ( auto& [site, m] : plaq->boundary ){
			if (site == self_site) continue;
            auto dp = static_cast<PyroSite*>(site);
            auto exp_ia = dp->xy();
            assert(m*m == 1);
			// performs ring -> ring * e^+- ia
            re_ring = re_ring*exp_ia.real() - im_ring * m* exp_ia.imag();
            im_ring = im_ring*exp_ia.real() + re_ring * m*exp_ia.imag();
        }
        return cplx(re_ring, im_ring);
    }

	// Calculates dH/dS+
	// H = sum 
	cplx local_field(const PyroSite& p){	
		cplx h = 0;

		// iterates over plaquettes sharing this link
		for (auto& [plaq_, m] : p.coboundary){
			auto plaq = static_cast<Plaquette*>(plaq_);
			// IMPORTANT --> Note m * g' to account for the sign of the plaquette in the direction it is traced
			h += partial_ring(plaq, &p)*cplx(plaq->g, m * plaq->g_prime);
		}

		return h;

	}

	double ring_energy(const Plaquette& plaq){
		// Returns the one-plauette part of the ring energy
        double re_ring=0;
        double im_ring=0;
        for ( auto& [site, m] : plaq.boundary ){
            auto dp = static_cast<PyroSite*>(site);
            auto exp_ia = dp->xy();
            assert(m*m == 1);
			// performs ring -> ring * e^+- ia
            re_ring = re_ring*exp_ia.real() - im_ring * m* exp_ia.imag();
            im_ring = im_ring*exp_ia.real() + re_ring * m*exp_ia.imag();
        }
		return 2*(plaq.g * re_ring - plaq.g_prime * im_ring);
	}



namespace MC {	


	inline double propose_Ising_tilt(double beta, std::mt19937& rng){
		double stdev = 0.25/sqrt(beta);
		// The Ising component can't be larger than 1 even if T is really large...
		if (stdev > 0.5){
			stdev = 0.5;
		}
		std::normal_distribution<double> normal(0, stdev);
		return normal(rng);
	}


	unsigned apply_step(sc_QSI& lat, double beta, std::mt19937& rng){
		unsigned success=0;
		apply_tetra_gauge(lat, rng);
		success += apply_angle(lat, beta, rng);
		success += apply_ring_Ising(lat, beta, rng);
		return success;
	}


	// Rotates all spins of sc_QSI in a gauge respecting way
	void apply_tetra_gauge(sc_QSI& lat, std::mt19937& rng) {
		std::uniform_real_distribution<> unif(-1,1);

		for (Tetra& t : lat.points){
			cplx theta(unif(rng), unif(rng));
			theta /= std::norm(theta);
			// theta is now a complex number of norm one
			for (auto& [spin, m] : t.coboundary){
				static_cast<PyroSite*>(spin)->xy()*=theta;
			}
		}
	}

	// Applies an angle update using Von Mises distribution
	unsigned apply_angle(sc_QSI& lat, double beta, std::mt19937& rng){	
		std::uniform_real_distribution<> unif(0,1);

		unsigned success = 0;

		for (PyroSite& p : lat.links){
			cplx theta(unif(rng)*2-1, unif(rng)*2-1);
			theta /= std::norm(theta);
			// theta is now a complex number of norm one
			
			cplx h = local_field(p);
			const cplx xy = p.xy();

			double rot = von_mises(rng, arg(h), abs(h * xy)*beta);

			// The distribution sampled is ~exp(beta|h||S| cos(arg(S)-arg(h)))
			double dE = -std::real( std::conj(local_field(p)) * p.xy() *
					(std::polar(1.0,rot) - 1.0) );
			if ((dE < 0) || (exp(-dE*beta) > unif(rng)) ) {
				p.xy() *= std::polar(1.0,rot);
				++success;
			}	
		}
		
		return success;
	}

	unsigned plaq_apply_MC_Ising(Plaquette& plaq, double old_plaq_state[], double beta, std::mt19937& rng){
		double E_old = ring_energy(plaq);
		std::uniform_real_distribution<>uniform;

		double tilt = propose_Ising_tilt(beta, rng);	

		double* state_ptr;

		int n=0;
		bool reject_move = false;

		// copy old state
		for ( const auto& [site, m] : plaq.boundary ){
			state_ptr = static_cast<PyroSite*>(site)->state;
			memcpy(state_ptr, old_plaq_state+3*n, 3*sizeof(double));
		}

		for (const auto& [site, m] : plaq.boundary ){

			state_ptr[3] += tilt * m;
			if (fabs(state_ptr[3]) > 1) {
				reject_move = true;
				break;
			}

			
			++n;
		}

		double dE;
		if (!reject_move){
			dE = ring_energy(plaq) - E_old; 
			reject_move = (dE>0) && (exp(-dE * beta) < uniform(rng));
		}

		if (reject_move){
			// restore the state of all the affected spins
			n=0;
			for ( const auto& [site, m] : plaq.boundary ){	
				memcpy(old_plaq_state+3*n, static_cast<PyroSite*>(site)->state, 3*sizeof(double));
				++n;
			}
			return 0;
		}
		return 1;
	}


	// The goal -- propose an MC Ising tilt that preserves the value of arg(ring())
	unsigned apply_ring_Ising(sc_QSI& lat, double beta, std::mt19937&rng){
		double old_plaq_state[3*6];

		unsigned success=0;
		for (Plaquette& plaq : lat.plaqs){
			assert(plaq.boundary.size() <= 6);
			success += plaq_apply_MC_Ising(plaq, old_plaq_state, beta, rng);
		}
		return success;
	}

};
