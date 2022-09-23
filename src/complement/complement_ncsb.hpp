#pragma once

#include "cola.hpp"
#include "complement_class.hpp"

namespace cola
{
    // (C, S, B) for complementing DACs
    const int NCSB_C = 2;
    const int NCSB_S = 4;
    const int NCSB_B = 3;

    class ncsb_compl : public complement_class
    {
    public:
        ncsb_compl(const spot::const_twa_graph_ptr aut, std::vector<unsigned> scc_index, spot::scc_info &scc_info, compl_decomp_options &decomp_options, unsigned true_index, std::vector<std::pair<unsigned, unsigned>> &dir_sim, std::vector<std::set<int>> &reachable_vector, std::vector<bool> &is_accepting) : complement_class(aut, scc_index, scc_info, decomp_options, true_index, dir_sim, reachable_vector, is_accepting) {}

        complement_mstate get_init_track();
        complement_mstate get_init_active();

        std::vector<std::pair<complement_mstate, bool>> get_succ_active(complement_mstate mstate, bdd symbol);
        std::vector<std::pair<complement_mstate, bool>> get_succ_track(complement_mstate mstate, bdd symbol);
        std::vector<std::pair<complement_mstate, bool>> get_succ_track_to_active(complement_mstate mstate, bdd symbol);
    };
}