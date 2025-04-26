#include <cell_geometry.hpp>
#include <argparse.hpp>
#include <chain.hpp>
#include <iostream>
#include <lattice_IO.hpp>
#include <preset_cellspecs.hpp>
#include <UnitCellSpecifier.hpp>
#include <algorithm>
#include <random>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>
#include <XoshiroCpp.hpp>
/**
 * dmndlat, a utility for generating diamond lattices precisely
 */

using namespace CellGeometry;


void parse_supercell_spec(imat33_t& supercell_spec, std::string Z1_s, std::string Z2_s, std::string Z3_s){

    std::stringstream Z1_ss(Z1_s);
    std::stringstream Z2_ss(Z2_s);
    std::stringstream Z3_ss(Z3_s);
    for (int row=0; row<3; row++){
        Z1_ss >> supercell_spec(row,0);
        Z2_ss >> supercell_spec(row,1);
        Z3_ss >> supercell_spec(row,2);
    }
}



inline std::string hash_parameters(
        const std::string& Z1_s, 
        const std::string& Z2_s, 
        const std::string& Z3_s
        ){
    // hashing the arguments 
    std::ostringstream name;
    name << "Z1="+Z1_s+";Z2="+Z2_s+";Z3="+Z3_s+";";
    auto s = name.str();
    std::replace(s.begin(), s.end(), ' ', ',');  // replace space by commas
    return s;
}

inline std::string hash_deletions(
        const std::vector<int>& dlist
        ){
    std::stringstream name_ss;
    if (dlist.size() > 0){
        // hash the setup for naming
        name_ss << "d1=";
        auto sep = "";
        for (const auto& j : dlist){
            name_ss << sep << j;
            sep = ",";
        }
        name_ss << ";";
    }

    return name_ss.str();
}


struct Tetra : public Cell<0> {
};

struct Spin : public Cell<1> {
    bool visited = false;
};

struct Plaq : public Cell<2> {
};

struct Vol : public Cell<3> {
    bool visited = false;
};

typedef PeriodicVolLattice<Tetra, Spin, Plaq, Vol> Lattice;

using namespace std;

struct bfs_node {
    const Cell<0>* point;
    Chain<1> path;
};



template <typename T>
std::vector<T*> get_neighbours(T* x0){
    std::vector<T*> retval;
    for (auto& [pl, m] : x0->boundary){
        for (auto& [x1, _] : pl->coboundary){
            if (x1 != x0) retval.push_back(static_cast<T*>(x1));
        }
    }
    return retval;
}


template <typename T>
std::vector<T*> get_coneighbours(T* x0){
    std::vector<T*> retval;
    for (auto& [pl, m] : x0->coboundary){
        for (auto& [x1, _] : pl->boundary){
            if (x1 != x0) retval.push_back(static_cast<T*>(x1));
        }
    }
    return retval;
}


inline std::vector<Chain<1>> find_paths_neighbours(Lattice& lat, Tetra* origin, Tetra* finish, unsigned len){
    // finds all chains of length 'len' connecting p0 to p1

    // ensure nothing is visited
    for (auto& l : lat.get_links()){
        l.visited = false;
    }
    std::stack<bfs_node> to_visit;
    to_visit.push(bfs_node(origin, Chain<1>()));

    Chain<0> expected_boundary;
    expected_boundary[origin] = 1;
    expected_boundary[finish] = -1;


#ifdef DEBUG
    cout<<"[find_paths_neighbours] expect boundary " << expected_boundary <<"\n";
    cout<<"\nposition          \t| len(path)\t| d(path)\n";
#endif

    std::vector<Chain<1>> res;
    while (!to_visit.empty()){
        auto curr = to_visit.top();
        cleanup_chain(curr.path);

        to_visit.pop();
        if (curr.path.size() == len){
            if (curr.point == finish) {
                assert( d(curr.path) == expected_boundary );
                res.push_back(curr.path);
            }
            continue;
        } 

#ifdef DEBUG
        cout<<curr.point->position<<"\t| "<< curr.path.size() <<"     \t| "<< d(curr.path)<<"\n";
#endif

        for (const auto& [l, m] : curr.point->coboundary){
            Spin* ln = static_cast<Spin*>(l);
            if (ln->visited) continue;
            ln->visited = true;
            Chain<1> this_ln; this_ln[l]=m;
            for (auto [p2, n] : d(this_ln)){
                if (p2 != curr.point){ 
                    to_visit.push(bfs_node(p2,curr.path + this_ln));
                }
            }
        }
    }
    return res;
}


inline std::vector<std::set<Vol*>> find_connected_components(Lattice& lat){
    // greedy algorithm -- traverse via all boundaries recursively
    std::stack<Vol*> vols;

    for (auto& [_,v] : lat.vols){
        v->visited = false;
    }

    std::vector<std::set<Vol*>> retval;

    for (auto& [_,v] : lat.vols){
        if (!v->visited){
            retval.push_back({});
            // start a DFS
            vols.push(v);
            while(!vols.empty()){
                auto curr =vols.top();
                retval.back().insert(curr);
                vols.pop();
                curr->visited = true;
                for(auto v2 :  get_neighbours<Vol>(curr)){
                    if (! v2->visited ){
                        vols.push(v2);
                    }
                }
            }
        }
    }

    return retval;
}

void sort_and_remove_duplicates(std::vector<int>& v){
    // Sort all of these so we don't muck up any indexing
    std::sort(v.begin(), v.end());
    // silently remove any duplicates
    v.erase( unique( v.begin(), v.end() ), v.end() );
}


// Deletes specified spin ids from the lattice
// Returns a std::set of Point* of 3 or less-member tetras
void del_spins_get_dtetras(Lattice& lat, std::vector<Spin*>& spins_to_delete, std::set<Tetra*>& defect_pts){
    // make sure it's in order
    std::sort(spins_to_delete.begin(), spins_to_delete.end());
    // delete the spins
    for (auto i=spins_to_delete.rbegin(); i<spins_to_delete.rend(); i++){
        auto l = *i;
        if (!lat.has_link(l)){
            throw std::out_of_range("Tried to delete spins at nonexistent index");
        }
        for (const auto& [p, _] : l->boundary){
            defect_pts.insert(static_cast<Tetra*>(p));
        }
        lat.erase_link(l);
    }
}


int main (int argc, const char *argv[]) {

    argparse::ArgumentParser prog("dmndlat");

    prog.add_argument("Z1")
        .help("First lattice vector in primitive units (as space separated string)");
    prog.add_argument("Z2")
        .help("First lattice vector in primitive units (as space separated string)");
    prog.add_argument("Z3")
        .help("First lattice vector in primitive units (as space separated string)");

    std::string outdir;
    prog.add_argument("--output_dir", "-o")
        .help("Path to output")
        .required()
        .store_into(outdir);

    
    prog.add_argument("--verbosity")
        .scan<'i', int>()
        .default_value(0);

    prog.add_argument("--neighbour")
        .scan<'i', int>()
        .default_value(2);
    
    std::vector<int> spin_ids_to_delete;
    prog.add_argument("--delete_spins", "-d")
        .help("Specific spin indexes to delete.")
        .default_value<std::vector<int>>({})
        .nargs(argparse::nargs_pattern::at_least_one)
        .store_into(spin_ids_to_delete);

    double dilution_prob;
    prog.add_argument("--dilution_prob","-p")
        .help("Probability of deleting spin i")
        .store_into(dilution_prob);

    uint64_t seed;
    prog.add_argument("--seed")
        .help("64-bit int to seed the RNG")
        .scan<'x', uint64_t>()
        .store_into(seed);

    prog.add_argument("--save_lattice")
        .help("Flag to save the full lattice file")
        .default_value(false)
        .implicit_value(true);
        

    try {
        prog.parse_args(argc, argv);
    } catch (const std::exception& err){
        cerr << err.what() << endl;
        cerr << prog;
        std::exit(1);
    }

    //////////////////////////////////////////////////////// 
    /// End program argument definitions
    ///

    std::filesystem::path outpath(outdir);
    if (! filesystem::exists(outpath) ){
        throw std::runtime_error("Cannot open outdir");
    }

    string Z1_s = prog.get<string>("Z1");
    string Z2_s = prog.get<string>("Z2");
    string Z3_s = prog.get<string>("Z3");


    // parse L1 L2 L3
    imat33_t supercell_spec;
    parse_supercell_spec(supercell_spec, Z1_s, Z2_s, Z3_s);
    std::cout<<"Constructing supercell of dimensions \n"<<supercell_spec<<std::endl;

    const PrimitiveSpecifers::Diamond spec;
    std::stringstream name;
    name << hash_parameters(Z1_s, Z2_s, Z3_s);

    auto neighbour = prog.get<int>("--neighbour");
    name <<"nn="<< neighbour << ";";

    Lattice lat(spec, supercell_spec);

    /// Erase the specified spins
    std::vector<Spin*> spins_to_yeet;
    if (spin_ids_to_delete.size() > 0){
        std::cout << "Erasing SPECIFIC spins\n";
        sort_and_remove_duplicates(spin_ids_to_delete);
        for (auto& x : spin_ids_to_delete){
            spins_to_yeet.push_back(lat.links.at(x));
        }
        name << hash_deletions(spin_ids_to_delete);
    }

    // Erase the randomly chosen spins
    // std::mt19937 gen(seed);
    XoshiroCpp::Xoshiro256PlusPlus gen(seed);
 
    std::bernoulli_distribution d1(dilution_prob);

    if (dilution_prob > 0){
        std::cout << "Erasing RANDOM spins (p="<<dilution_prob<<")\n";
        for (const auto& [_, p] : lat.links) {
            if (d1(gen)) spins_to_yeet.push_back(p);
        }
        name<<"p="<<dilution_prob<<";";
    }

    std::set<Tetra*> defect_tetras;
    del_spins_get_dtetras(lat, spins_to_yeet, defect_tetras);


    auto verbosity = prog.get<int>("--verbosity");
    lat.print_state(verbosity);
      
    
    if (verbosity >= 3){
        // sanity check: ensure all tetras are really defective
        std::cout<<"Defect tetras:\n";
        for (const auto t : defect_tetras){
            assert(t->coboundary.size() < 4);
            std::cout<<"[Dtetra] "<<t->position<<"\n";
        }
    }

    if (verbosity >= 1){
        std::cout<<"\n Finding links...\n";
    }

    // note that link erasures are guaranteed not to delete any points.
    vector<Tetra*> defect_tetras_vec(defect_tetras.begin(), defect_tetras.end());

    unsigned total_n = defect_tetras_vec.size() * (defect_tetras_vec.size() -1)/2;
    unsigned print_counter = 0ul;

    std::vector<ipos_t> deleted_link_locs;
    for (unsigned i=0; i<defect_tetras_vec.size(); i++){
        for (unsigned j=0; j<i; j++){
            auto t1 = defect_tetras_vec[i];
            auto t2 = defect_tetras_vec[j];

            auto links = find_paths_neighbours(lat, t1, t2, neighbour);
#ifdef DEBUG
            std::cout<<"("<<t1->position<<", "<<t2->position<<")\n";
            std::cout<<" connection size: "<<links.size()<<std::endl;
#endif
            for (auto path : links){
                for (auto [l, _] : path){
                    auto s = static_cast<Spin*>(l);
                    if (lat.has_link(s)){
                        deleted_link_locs.push_back(s->position);
                        lat.erase_link(s);
                    }
                }
            }
            if (print_counter % 100 == 0){
                printf("%5d / %5d (%02d%%)\r", print_counter, total_n, print_counter * 100 / total_n);
            }
            print_counter++;
        }
    }
    printf("\n");

    auto path = outpath/(name.str()+".lat.json");
    if (prog.get<bool>("--save_lattice")){
        cout<<"Saving to "<<path<<std::endl;

        using namespace nlohmann;
        json j = {};
        write_data(lat, j); // write the lattice data
        std::ofstream of(path); 
        of << j;
        of.close();
        j["defect_link_locs"] = deleted_link_locs;
    }
    return 0;
}


