#pragma once

#include "cola.hpp"
#include "complement_class.hpp"
#include "rankings.hpp"

namespace cola
{   
    class rank_compl : public complement_class
    {
    private:
        std::vector<std::set<unsigned>> predecessors_;
        waiting &waiting_;
        std::map<std::set<unsigned>, unsigned> rank_restr_;

    public:
        rank_compl(const spot::const_twa_graph_ptr aut, std::vector<unsigned> scc_index, spot::scc_info &scc_info, compl_decomp_options &decomp_options, unsigned true_index, std::vector<std::pair<unsigned, unsigned>> &dir_sim, std::vector<std::set<int>> &reachable_vector, std::vector<bool> &is_accepting, waiting &wait) : complement_class(aut, scc_index, scc_info, decomp_options, true_index, dir_sim, reachable_vector, is_accepting), waiting_(wait)
        {
            // predecessors
            std::vector<std::set<unsigned>> predecessors(reachable_vector_.size());
            for (unsigned i = 0; i < reachable_vector.size(); i++)
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

            // rank restriction
            auto states_in_scc = scc_info.states_of(scc_index[0]).size();
            for (auto mstate : waiting_.get_states())
            {
                // initialization
                unsigned nonacc = 0;
                for (auto state : mstate)
                {
                    if (not is_accepting[state])
                        nonacc++;
                }
                rank_restr_.insert({mstate, 2*(nonacc + 1)});
            }

            // outer macrostate analysis
            auto wait_pred = waiting_.get_predecessors();
            bool change = true;

            while (change)
            {
                change = false;
                for (auto state : waiting_.get_states())
                {
                    unsigned bound = rank_restr_[state];
                    unsigned max = 0;
                    for (auto pred : wait_pred[state])
                    {
                        if (rank_restr_[pred] > max)
                        {
                            max = rank_restr_[pred];
                        }
                    }
                    if (max < bound)
                    {
                        rank_restr_[state] = max;
                        change = true;
                    }
                }
            }
        }

        complement_mstate get_init_track();
        complement_mstate get_init_active();
        
        std::vector<std::pair<complement_mstate, bool>> get_succ_active(complement_mstate mstate, bdd symbol);
        std::vector<std::pair<complement_mstate, bool>> get_succ_track(complement_mstate mstate, bdd symbol);
        std::vector<std::pair<complement_mstate, bool>> get_succ_track_to_active(complement_mstate mstate, bdd symbol);

        std::set<int> get_successors_with_box(std::set<unsigned> reachable, rank_state state, bdd letter, unsigned scc_index);

        void get_max(std::vector<ranking> &rankings);
        std::vector<ranking> get_maxrank(std::set<unsigned> reachable, rank_state state, bdd letter, unsigned scc_index);
        std::vector<ranking> get_succ_rankings(ranking r, std::vector<std::tuple<int, int, bool>> restr, std::set<unsigned> reachable, bdd letter, unsigned scc_index);
    };
}