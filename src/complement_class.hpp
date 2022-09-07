#pragma once

#include "cola.hpp"
#include "complement_mstate.hpp"

namespace cola
{
    class complement_class
    {   
        protected:
        const spot::const_twa_graph_ptr aut_;
        std::vector<unsigned> scc_index_;
        spot::scc_info &scc_info_;
        complement_mstate &mstate_;
        compl_decomp_options decomp_options_;
        bdd symbol_;
        unsigned true_index_;
        std::vector<std::pair<unsigned, unsigned>> dir_sim_;
        std::vector<std::set<int>> reachable_vector_;

        public:
        complement_class(const spot::const_twa_graph_ptr aut, std::vector<unsigned> scc_index, spot::scc_info &scc_info, complement_mstate &mstate, compl_decomp_options decomp_options, bdd symbol, unsigned true_index, std::vector<std::pair<unsigned, unsigned>> dir_sim, std::vector<std::set<int>> reachable_vector) : aut_(aut), scc_index_(scc_index), scc_info_(scc_info), mstate_(mstate), decomp_options_(decomp_options), symbol_(symbol), true_index_(true_index), dir_sim_(dir_sim), reachable_vector_(reachable_vector) {}

        virtual complement_mstate getInit() = 0;
        std::vector<std::pair<complement_mstate, bool>> get_succ_track(){return std::vector<std::pair<complement_mstate, bool>>();}
        std::vector<std::pair<complement_mstate, bool>> get_succ_track_to_active(){return std::vector<std::pair<complement_mstate, bool>>();}
        std::vector<std::pair<complement_mstate, bool>> get_succ_active(){return std::vector<std::pair<complement_mstate, bool>>();}

        std::set<unsigned> get_all_successors(std::vector<unsigned> current_states, bdd symbol); 
        std::set<unsigned> get_all_successors_in_scc(std::vector<unsigned> current_states, bdd symbol); 
        std::set<unsigned> get_succ_acc_trans_scc(std::vector<unsigned> current_states, bdd symbol);
    };
}