#pragma once

#include "cola.hpp"
#include "rankings.hpp"

namespace cola
{
    class rank_complement
    {
        private:
            const spot::const_twa_graph_ptr aut_;
            spot::scc_info scc_info_;
            std::string scc_types_;
            compl_decomp_options decomp_options_;
            std::vector<std::pair<unsigned, unsigned>> dir_sim_;
            std::vector<std::set<int>> reachable_vector_;
            std::vector<std::set<unsigned>> predecessors_;
            std::vector<bool> is_accepting_;

        public:
            rank_complement(const spot::const_twa_graph_ptr &aut, spot::scc_info &scc_info, std::string scc_types, compl_decomp_options decomp_options, std::vector<std::pair<unsigned, unsigned>> dir_sim, std::vector<std::set<int>> reachable_vector, std::vector<bool> is_accepting) : aut_(aut), scc_info_(scc_info), scc_types_(scc_types), decomp_options_(decomp_options), dir_sim_(dir_sim), reachable_vector_(reachable_vector), is_accepting_(is_accepting) {
                
                std::vector<std::set<unsigned>> predecessors(scc_types_.size());
                for (unsigned i=0; i<reachable_vector.size(); i++)
                {
                    unsigned scc_index = scc_info.scc_of(i);
                    for (auto state : reachable_vector[i])
                    {
                        unsigned succ_index = scc_info.scc_of(state);
                        if (scc_index != succ_index)
                            predecessors[succ_index].insert(scc_index);
                    }
                }

                predecessors_ = predecessors;
            }

            std::set<unsigned> get_all_successors(std::set<unsigned> current_states, bdd symbol);
            std::set<int> get_successors_with_box(std::set<unsigned> reachable, rank_state state, bdd letter);

            void get_max(std::vector<ranking>& rankings);

            std::vector<std::pair<rank_state, bool>> get_succ_track(std::set<unsigned> reachable, rank_state state, bdd letter);
            std::vector<std::pair<rank_state, bool>> get_succ_track_to_active(std::set<unsigned> reachable, rank_state state, bdd letter);
            std::vector<std::pair<rank_state, bool>> get_succ_active(std::set<unsigned> reachable, rank_state state, bdd letter, bool one_scc);
    };
}