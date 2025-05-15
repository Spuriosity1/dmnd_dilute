#include <cassert>
#include <cell_geometry.hpp>
#include <argparse.hpp>
#include <chain.hpp>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <lattice_IO.hpp>
#include <ostream>
#include <preset_cellspecs.hpp>
#include <UnitCellSpecifier.hpp>
#include <algorithm>
#include <queue>
#include <random>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>
#include <XoshiroCpp.hpp>
#include "format_bits.hpp"
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
    //bool visited = false;
    Tetra* origin = nullptr;
};

struct Plaq : public Cell<2> {
    bool visited = false;
    const Plaq* root = nullptr;
};

struct Vol : public Cell<3> {
    bool visited = false;
    const Vol* root = nullptr;
};

typedef PeriodicVolLattice<Tetra, Spin, Plaq, Vol> Lattice;


struct search_node {
    const Cell<0>* point;
    Chain<1> path;
};



/*
inline std::vector<Chain<1>> find_paths_neighbours(Lattice& lat, Tetra* origin, Tetra* finish, unsigned len){
    // finds all chains of length 'len' connecting p0 to p1

    // ensure nothing is visited
    for (auto& l : lat.get_links()){
        l.visited = false;
    }
    std::stack<search_node> to_visit;
    to_visit.push(search_node(origin, Chain<1>()));

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
                    to_visit.push(search_node(p2,curr.path + this_ln));
                }
            }
        }
    }
    return res;
}
*/

void excise_path(Lattice& lat, Chain<1>& path, std::set<void*>& deleted_link_ptrs, std::vector<ipos_t>& deleted_link_locs){
    for (auto [l, _] : path){
        auto s = static_cast<Spin*>(l);
        //if (lat.has_link(s)){
        if (!deleted_link_ptrs.contains(s)){
//            assert(lat.has_link(s));
            deleted_link_locs.push_back(s->position);
            lat.erase_link(s);
            deleted_link_ptrs.insert(s);
        }
    }
}

inline std::vector<Chain<1>> find_defect_links(
        Tetra* origin, unsigned len ){
    /** 
     * Finds all paths of specified length(s) connecting origin to a 
     * defect node and removes them from the lattice
     * @param lat: the base lattice
     * @param origin: the starting point
     * @param lens_to_trim: the set of lattice-seps to delete 
     * measured as the number of Tetras that are part of the path
     * EXCLUDING start (len=1 correspnds to nearest-neighbour pyrochlore sites)
     */

    /*
    // ensure nothing is visited
    for (auto& l : lat.get_links()){
        //l.visited = false;
    }
    */

    std::queue<search_node> to_visit;
    to_visit.push(search_node(origin, Chain<1>()));

    std::vector<Chain<1>> res;
    while (!to_visit.empty()){
        auto& curr = to_visit.front();

        // Check if we are at another defect point
        if (curr.path.size() == len){
            if (curr.point->coboundary.size() < 4){
                // trimmable path: mark it for deletion
                cleanup_chain(curr.path);
                res.push_back(curr.path);
            }
            to_visit.pop(); // curr invalidated!
            continue;
        }

        assert (curr.path.size() < len);
#ifdef DEBUG
        cout<<curr.point->position<<"\t| "<< curr.path.size() <<"     \t| "<< d(curr.path)<<"\n";
#endif
        for (const auto& [l, m] : curr.point->coboundary){
            Spin* ln = static_cast<Spin*>(l);
            if (ln->origin == origin) continue;
            ln->origin = origin;
            Chain<1> this_ln; this_ln[l]=m;
            for (auto [p2, n] : d(this_ln)){
                if (p2 != curr.point){ 
                    to_visit.push(search_node(p2,curr.path + this_ln));
                }
            }
        }
        to_visit.pop();
    }

    return res;
}


template<typename T>
concept Visitable = requires(T t) {
    // { t.visited } -> std::same_as<bool&>;
    { t.root } -> std::same_as<const T*&>;
};

template<Visitable T>
struct conn_components {
    std::set<T*> elems;
    bool wraps = false;
};

template<Visitable T>
inline std::vector<conn_components<T>> find_connected(std::unordered_map<sl_t, T*>& elems){
    /**
     * Starting from all points, either start a new cluster or expand cluster
     * until unable.
     */
    std::stack<T*> stack;

    // mark all unvisited
    for (auto& [_, v] : elems){
         // v->visited = false;
        v->root = nullptr;
    }

    std::vector<conn_components<T>> retval;

    for (auto& [_, root_cell] : elems){
        if (root_cell->root != nullptr) continue;

        retval.push_back({});
        // start a DFS
        stack.push(root_cell);
        while(!stack.empty()){
            auto curr = stack.top();
            retval.back().elems.insert(curr);
            stack.pop();
            curr->root = root_cell;
            for(auto next : get_neighbours<T>(curr)){
                if (next->root == nullptr){
                    stack.push(next);
                } else if (next->root == root_cell) {
                    // we have linked with the hive! Check if it was around a boundary -
                    auto dx1 = curr->position - root_cell->position;
                    auto dx2 = next->position - root_cell->position;
                    rational::Rational r = dot(dx1, dx2);
                    if ( r.num*r.denom < 0  ) {
                        retval.back().wraps = true;
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
    std::unordered_set<Spin*> present_spins;
    for (auto& [_, l] : lat.links){
        present_spins.insert(l);
    }
    // make sure it's in order
    std::sort(spins_to_delete.begin(), spins_to_delete.end());
    // delete the spins
    for (auto i=spins_to_delete.rbegin(); i<spins_to_delete.rend(); i++){
        auto l = *i;
        //if (!lat.has_link(l)){
        if (!present_spins.contains(l)){
            throw std::out_of_range("Tried to delete spins at nonexistent index");
        }
        for (const auto& [p, _] : l->boundary){
            defect_pts.insert(static_cast<Tetra*>(p));
        }
        lat.erase_link(l);
        present_spins.erase(l);
    }
}

/////////////////////////////////
/// EXPORTING TO FILE ///////////

void export_lattice(
        const filesystem::path& path,
        const Lattice& lat,
        const std::vector<ipos_t>& deleted_link_locs
        ){


        cout<<"Saving lattice to \n"<<path<<std::endl;

        json j = {};
        write_data(lat, j); // write the lattice data to j

        j["defect_link_locs"] = deleted_link_locs;

        std::ofstream of(path); 
        of << j;
        of.close();
}


inline json latstats_to_json(const Lattice& lat){
    json counts = {};
    counts["points"] = lat.points.size();
    counts["links"] = lat.links.size();
    counts["plaqs"] = lat.plaqs.size();
    counts["vols"] = lat.vols.size();

    return counts;
}

template <typename T>
inline std::vector<size_t> get_sorted_sizes(const std::vector<conn_components<T>>& parts ){
    std::vector<size_t> size_set;
    for (const auto& p : parts){
        auto size = p.elems.size();
        auto it = std::lower_bound(size_set.begin(), size_set.end(), size);
        size_set.insert(it, size);
    }
    return size_set;
}

template <typename T>
inline bool test_wraps(const std::vector<conn_components<T>>& parts ){
    for (const auto& p : parts){
        if (p.wraps) { return true; }
    }
    return false;
}


inline json percolstats_to_json(
        const std::vector<conn_components<Plaq>> connected_plaqs,
        const std::vector<conn_components<Vol>> connected_vols
        ) {

    json percolstats = {};
    percolstats["n_plaq_parts"] = connected_plaqs.size();
    percolstats["plaq_part_nelem"] = get_sorted_sizes(connected_plaqs);
    percolstats["plaqs_wrap"] = test_wraps(connected_plaqs);

    percolstats["n_vol_parts"] = connected_vols.size();
    percolstats["vol_part_nelem"] = get_sorted_sizes(connected_vols);
    percolstats["vols_wrap"] = test_wraps(connected_vols);
    return percolstats;
}



void export_stats(
        const filesystem::path& path,
        const Lattice& lat,
        const std::vector<conn_components<Plaq>> connected_plaqs,
        const std::vector<conn_components<Vol>> connected_vols
        ){
    cout<<"Saving statistics to \n"<<path<<std::endl;


    json j = {};


    j["counts"] = latstats_to_json(lat);
    j["percolation"] = percolstats_to_json(connected_plaqs, connected_vols);

    std::ofstream of(path); 
    of << j;
    of.close();
}


int main (int argc, const char *argv[]) {

    argparse::ArgumentParser prog("dmndlat");
    prog.add_argument("Z1")
        .help("First lattice vector in primitive units (as space separated string)")
        .nargs(3)
        .scan<'i', int>();
    prog.add_argument("Z2")
        .help("Second lattice vector in primitive units (as space separated string)")
        .nargs(3)
        .scan<'i', int>();
    prog.add_argument("Z3")
        .help("Third lattice vector in primitive units (as space separated string)")
        .nargs(3)
        .scan<'i', int>();

    std::string outdir;
    prog.add_argument("--output_dir", "-o")
        .help("Path to output")
        .required()
        .store_into(outdir);
    
    prog.add_argument("--verbosity", "-v")
        .scan<'i', int>()
        .default_value(0);

    prog.add_argument("--force", "-f")
        .help("Overwrites output files")
        .default_value(false)
        .implicit_value(true);

    std::vector<int> neighbours;
    prog.add_argument("--neighbours", "-n")
        .scan<'i', int>()
        .nargs(argparse::nargs_pattern::at_least_one)
        .default_value<std::vector<int>>({2})
        .store_into(neighbours);
    
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

    std::string seed_s;
    prog.add_argument("--seed")
        .help("64-bit int to seed the RNG")
        .store_into(seed_s);

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


    uint64_t seed; // ugly hack
    std::stringstream ss;
    ss << std::hex << seed_s;
    ss >> seed;


    std::stringstream name; // accumulates hashed options

    // parse L1 L2 L3
    imat33_t supercell_spec;
    name << parse_supercell_spec(supercell_spec, prog);
    std::cout<<"Constructing supercell of dimensions \n"<<supercell_spec<<std::endl;

    const PrimitiveSpecifers::Diamond spec;
 
    name << comma_separate("nn", neighbours);


    Lattice lat(spec, supercell_spec);

    /// Erase the specified spins
    std::vector<Spin*> spins_to_yeet;
    if (spin_ids_to_delete.size() > 0){
        std::cout << "Erasing SPECIFIC spins\n";
        sort_and_remove_duplicates(spin_ids_to_delete);
        for (auto& x : spin_ids_to_delete){
            spins_to_yeet.push_back(lat.links.at(x));
        }
        name << comma_separate("d1", spin_ids_to_delete);
    }

    // Erase the randomly chosen spins
    // std::mt19937 gen(seed);
    XoshiroCpp::Xoshiro256PlusPlus gen(seed);
 
    std::bernoulli_distribution d1(dilution_prob);

    for (const auto& [_, p] : lat.links) {
        if (d1(gen)) spins_to_yeet.push_back(p);
    }

    char buf[1024];
    snprintf(buf, 1024, "p=%.04f;seed=%llx;", dilution_prob, seed);
    name<<buf;

    bool save_lattice = prog.get<bool>("--save_lattice");

    // Name now fully specified
    auto statpath = outpath/(name.str()+".stats.json");
    auto latpath = outpath/(name.str()+".lat.json");

    if (!prog.get<bool>("--force")){
        // Check if these files already exist, if so abort early
        if ( !save_lattice && filesystem::exists(statpath)){
            cerr << "Statfile " << statpath << "already exists" << std::endl;
            throw std::runtime_error("Statfile exists");
        }
        if ( save_lattice && filesystem::exists(latpath)){
            cerr << "latfile " << latpath << "already exists" << std::endl;
            throw std::runtime_error("Latfile exists");
        }
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

    unsigned total_n = defect_tetras_vec.size();

    std::vector<ipos_t> deleted_link_locs;

    // initialise the spins to have nullptr roots
    for (auto& [_, s] : lat.links){
        s->origin = nullptr;
    }
    
    for (auto len : neighbours){
        printf("[search] finding %d neighbours\n", len);
        unsigned print_counter = 0ul;
        for (auto t1 : defect_tetras_vec){
            auto links = find_defect_links(t1, len);
            std::set<void*> deleted_link_ptrs;
            
            for (auto path : links){
                excise_path(lat, path, deleted_link_ptrs, deleted_link_locs);
            }
            
            printf("%5d / %5d (%02d%%)\r", print_counter, total_n,
                    print_counter * 100 / total_n);
            fflush(stdout); 
            print_counter++;
            
        }

        printf("\n");
    }
    

    if (save_lattice){
        export_lattice(latpath, lat, deleted_link_locs);
    }

    // Counting complete, output observables
    // Simple ones: raw counts

    
    // more complex: connected components 
    auto connected_plaqs = find_connected(lat.plaqs);
    auto connected_vols = find_connected(lat.vols);


    export_stats(statpath, lat, connected_plaqs, connected_vols);

    return 0;
}

