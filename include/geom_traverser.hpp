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
};

template<Visitable T>
struct conn_components {
    std::set<T*> elems;
    bool wraps = false;
};

inline ipos_t min(ipos_t a, ipos_t b){
    ipos_t c(a);
    for (int i=0; i<3; i++){
        if (c[i] > b[i]) c[i] = b[i];
    }
    return c;
}


inline ipos_t max(ipos_t a, ipos_t b){
    ipos_t c(a);
    for (int i=0; i<3; i++){
        if (c[i] < b[i]) c[i] = b[i];
    }
    return c;
}

inline bool any_wraps(ipos_t X, double Lmin){
    bool rv = 0;
    for (int i=0; i<3; i++){
        rv |= (X[i]*X[i] > Lmin/2 );
    }
    return rv;
}


template<Visitable T>
inline std::vector<conn_components<T>> find_connected(
		CellGeometry::SparseMap<CellGeometry::sl_t, T*>& elems, 
		double Lmin){
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
        auto& cell_union = retval.back();

        ipos_t bb_min(root_cell->position);
        ipos_t bb_max(root_cell->position);

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
                } else if (next->root == root_cell) {
                    // we have linked with the hive! Update the bounding box
                    auto dx1 = curr->position - next->position;
                    if (any_wraps(dx1, Lmin)){
                        bb_min = min(bb_min, curr->position);
                        bb_min = min(bb_min, next->position);

                        bb_max = max(bb_max, curr->position);
                        bb_max = max(bb_max, next->position);
                    }
                }
            }
        }

        auto delta_bb = bb_max-bb_min;
        for (int i=0; i<3; i++){ 
            if (delta_bb[i] > Lmin - 2) {
                cell_union.wraps=true;
            }
        }        
    }

    return retval;
}


inline double calc_Lmin(const imat33_t& cell_vectors){
    int64_t l2[3] = {0,0,0};
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            l2[i] += cell_vectors(j,i)*cell_vectors(j,i);
        }
    }
    auto Lmin2 = std::min(l2[0], std::min(l2[1],l2[2]));
    return sqrt(1.0*Lmin2);
}
