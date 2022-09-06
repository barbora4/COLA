#pragma once

#include "cola.hpp"
#include "complement_class.hpp"

namespace cola 
{
    class rank_compl : public complement_class
    {
        public:
        rank_compl(const spot::const_twa_graph_ptr aut, std::vector<unsigned> scc_index, spot::scc_info &scc_info, complement_mstate &mstate, compl_decomp_options decomp_options, bdd symbol, unsigned true_index, std::vector<std::pair<unsigned, unsigned>> dir_sim, std::vector<std::set<int>> reachable_vector) : complement_class(aut, scc_index, scc_info, mstate, decomp_options, symbol, true_index, dir_sim, reachable_vector) {}
    }; 
}