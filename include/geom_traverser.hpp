#pragma once
#include <cell_geometry.hpp>
#include <stack>
#include <set>
#include <algorithm>
#include <UnitCellSpecifier.hpp>
#include <chain.hpp>

template<typename T>
concept Visitable = requires(T t) {
    // { t.visited } -> std::same_as<bool&>;
    { t.root } -> std::same_as<const T*&>;
    { t.position } -> std::convertible_to<const ipos_t>;
};

template<Visitable T>
struct conn_components {
    std::set<T*> elems;
    bool wraps = false;
};


inline double calc_Lmin2(const imat33_t& cell_vectors){
    int64_t l2[3] = {0,0,0};
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            l2[i] += cell_vectors(j,i)*cell_vectors(j,i);
        }
    }
    return std::min(l2[0], std::min(l2[1],l2[2]));
}

template<Visitable T>
inline std::vector<conn_components<T>> find_connected(
        CellGeometry::PeriodicAbstractLattice lat,
		CellGeometry::SparseMap<CellGeometry::sl_t, T*>& elems
		){
    /**
     * Starting from all points, either start a new cluster or expand cluster
     * until unable.
     * @param elems some collection of gometric objects, stored as a
     *                  Map from a sublattice index to a pointer to the element.
     *  @param Lmin The minimum physical distance that a wrapping cluster must cover. 
     */


    std::stack<T*> stack;
    // Classic union-find algorithm. The Visitbale concept ensures that type T 
    // has a pointer "root" that can be used to keep track of cluster ownership.
    
    // Init: mark all unvisited
    for (auto& [_, v] : elems){
        v->root = nullptr;
    }

    std::vector<conn_components<T>> component_list;
    
    // Iterate through all elements, eznsuring that everyone gets a visit. 
    // We "colour" each element with a non-null pointer to a root element
    // Each nullptr node we visit is therefore the beginning of a new cluster.
    for (auto& [_, root_cell] : elems){
        if (root_cell->root != nullptr) continue; // already visited, skip

        component_list.push_back({});
        auto& cell_union = component_list.back();
        
        // start a DFS
        stack.push(root_cell);
        while(!stack.empty()){
            auto curr = stack.top();
            cell_union.elems.insert(curr);
            stack.pop();

            curr->root = root_cell;
            for(auto next : CellGeometry::get_neighbours<T>(curr)){
                if (next->root == nullptr){
                    stack.push(next);
                }
            }
        }
    }

    auto Lmin2 = calc_Lmin2(lat.cell_vectors);

    // All components identified. Second pass: determine if these wrap
    for (auto& component : component_list){
        if (component.elems.empty()) { continue; }
        auto root_cell = (*component.elems.begin())->root;
        auto x0 = root_cell -> position;
        
        for (auto el = component.elems.rbegin(); 
                el != component.elems.rend(); el++){
//        for (auto el : component.elems) {
            if (lat.d2((*el)->position, x0) > Lmin2/4) {
                component.wraps = true;
                break;
            }
        }
    }

    return component_list;
}


