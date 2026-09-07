#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <unordered_set>
#include <argparse.hpp>

// LatticeLab
#include <chain.hpp>
#include <cell_geometry.hpp>
#include <preset_cellspecs.hpp>
#include <UnitCellSpecifier.hpp>

#include "format_bits.hpp"
#include "geom_traverser.hpp"
/**
 * Adds link disorder to a diaomnd lattice and removes any 
 * even length intermediaries.
 */


using namespace CellGeometry;
using namespace nlohmann;
using namespace std;


struct Tetra : public Cell<0> {
};

struct Spin : public Cell<1> {
    unsigned depth=0;
};

struct Plaq : public Cell<2> {
    // 
};

struct Vol : public Cell<3> {
    //
};

typedef PeriodicVolLattice<Tetra, Spin, Plaq, Vol> Lattice;



inline std::map<int, int> count_link_neighbours(Lattice& lat, Spin* origin, unsigned len) {
    /** 
     * Counts the total number of Spin neighbours 
     * separated from origin by a length of <=len. 
     * len=0 -> 1 (only origin itself)
     * len=1 -> 6 (nearest neighbours in pyrochlore language)
     */
    std::map<int, int> hist;

    for (unsigned i=0; i<=len; i++){
        hist[i]=0;
    }

    // algo: colour all spins by depth from origin.
    // Initialise all depths to be infinity
    for (auto& [sl, l]: lat.links){
        l->depth = std::numeric_limits<unsigned>::max();
    }

    std::queue<Spin*> to_visit;
    std::unordered_set<Tetra*> seen_tetras;

    to_visit.push(origin);
    origin->depth=0;
    hist[0]=1; // origin itself sits at depth 0
    while(!to_visit.empty()){
        auto curr_l = to_visit.front();
        to_visit.pop();
        for (auto& [tt, m] : curr_l->boundary){
            auto t = static_cast<Tetra*>(tt);
            if (seen_tetras.contains(t)) continue;
            seen_tetras.insert(t);

            for (auto& [l, m2] : t->coboundary){
                auto new_l = static_cast<Spin*>(l);
                // skip anything already coloured (including origin/back-tracking);
                // a spin can be reached via either of its two boundary tetras.
                if (new_l->depth != std::numeric_limits<unsigned>::max()) continue;
                new_l->depth = curr_l->depth+1;

                if (new_l->depth <= len)
                    hist[new_l->depth]++;
                if (new_l->depth < len)
                    to_visit.push(new_l);
            }

        }
    }
    return hist;
}



int main (int argc, const char *argv[]) {

    argparse::ArgumentParser prog(argv[0]);
    prog.add_argument("L")
        .help("Linear dimension of supercell for search")
        .scan<'i',int>();
    prog.add_argument("maxN")
        .help("max length to search up to")
        .scan<'i',int>();

    prog.parse_args(argc, argv);
    int L = prog.get<int>("L");
    int maxN = prog.get<int>("maxN");


    imat33_t supercell_spec = imat33_t::from_cols({L,-L,-L},{-L,L,-L},{-L,-L,L});
    std::cout<<"Constructing supercell of dimensions \n"<<supercell_spec<<std::endl;

    const auto spec = PrimitiveSpecifiers::DiamondSpec();
    Lattice lat(spec, supercell_spec);

    auto counts = count_link_neighbours(lat, lat.links.begin()->second, maxN);

    printf("Count summary: \n");
    for (auto [n, count]: counts){
        printf("%4d\t%8d\n", n,count);
    }
    
    return 0;

}
