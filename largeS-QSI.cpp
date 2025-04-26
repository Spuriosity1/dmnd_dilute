#include <cstdlib>
#include <ostream>
#include <random>
#include <sstream>
#include "admin/basic_parser.hh"
#include "semiclassic_qsi.hpp"

using namespace std;


void parse_supercell_spec(imat33_t& supercell_spec,
        std::string L1_s, std::string L2_s, std::string L3_s){

    std::stringstream L1_ss(L1_s);
    std::stringstream L2_ss(L2_s);
    std::stringstream L3_ss(L3_s);
    for (int row=0; row<3; row++){
        L1_ss >> supercell_spec(row,0);
        L2_ss >> supercell_spec(row,1);
        L3_ss >> supercell_spec(row,2);
    }
}

int main (int argc, const char *argv[]) {
    std::string L1_s,L2_s,L3_s;
    unsigned num_anneal;
    double T_hot;
    std::filesystem::path outpath;

    basic_parser::Parser args(1,2);
    args.declare("Z1",&L1_s);
    args.declare("Z2",&L2_s);
    args.declare("Z3",&L3_s);

    args.declare("num_anneal", &num_anneal);
    args.declare("T_hot", &T_hot);

    if (argc < 3){
        cout<<"USAGE: "<<argv[0]<<" <input file> <output directory> [overrides...]";
    }

    args.from_file(argv[1]);
    outpath = argv[2];
    args.from_argv(argc, argv, 3);

    args.assert_initialised();

    // parse L1 L2 L3
    imat33_t supercell_spec;
    parse_supercell_spec(supercell_spec, L1_s, L2_s, L3_s);

    // Import and construct the cell specifier (this can be done statically
    // in principle, but there is very little point optimising it)
    const PrimitiveSpecifers::Diamond pyrospec;

    // Construct the periodic lattice extension using these extensions
    sc_QSI lat(pyrospec, supercell_spec);

    // Initialise RNG
    static std::random_device dev;
    static auto rng = std::mt19937(dev());

    // Run MC steps 
    unsigned num_success = 0;
    for (unsigned n=0; n<num_anneal; n++) {
        num_success += MC::apply_step(lat, 1/T_hot, rng);
    }
    cout << "Success rate: " << num_success*100.0 / num_anneal <<std::endl;


    return 0;
}
