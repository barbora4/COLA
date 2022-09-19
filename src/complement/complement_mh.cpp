#include "complement_mh.hpp"

namespace cola
{
    void mh_compl::prune_track(std::vector<unsigned> &succ_in_scc)
    {
        // remove smaller states from S
        std::set<unsigned> new_S(succ_in_scc.begin(), succ_in_scc.end());

        for (auto pr : dir_sim_)
        {
            if (pr.first != pr.second and new_S.find(pr.first) != new_S.end() and new_S.find(pr.second) != new_S.end())
            {
                // reachability check
                if (reachable_vector_[pr.second].find(pr.first) == reachable_vector_[pr.second].end())
                {
                    // both states in S -> we can remove the smaller one from S
                    new_S.erase(pr.first);
                }
            }
        }
        succ_in_scc = std::vector<unsigned>(new_S.begin(), new_S.end());
    }

    void mh_compl::prune_active(std::set<unsigned> &succ_in_scc, std::set<unsigned> &new_break_set)
    {
        // remove smaller states from S
        std::set<unsigned> new_S(succ_in_scc.begin(), succ_in_scc.end());
        for (auto pr : dir_sim_)
        {
            if (pr.first != pr.second and new_S.find(pr.first) != new_S.end() and new_S.find(pr.second) != new_S.end())
            {
                if (reachable_vector_[pr.second].find(pr.first) == reachable_vector_[pr.second].end())
                {
                    // both states in S -> we can remove the smaller one from S
                    new_S.erase(pr.first);

                    // remove the states also from new_break_set if present
                    if (new_break_set.find(pr.first) != new_break_set.end())
                        new_break_set.erase(pr.first);
                }
            }
        }
        succ_in_scc = new_S;
    }

    std::vector<std::pair<complement_mstate, bool>> mh_compl::get_succ_track(complement_mstate mstate, bdd symbol)
    {
        std::vector<std::pair<complement_mstate, bool>> succ;

        std::set<unsigned> all_succ = this->get_all_successors(mstate.curr_reachable_, symbol);
        std::vector<unsigned> succ_in_scc;
        std::set<unsigned> scc_states;
        for (auto i : scc_index_)
        {
            scc_states.insert(scc_info_.states_of(i).begin(), scc_info_.states_of(i).end());
        }
        std::set_intersection(all_succ.begin(), all_succ.end(), scc_states.begin(), scc_states.end(), std::back_inserter(succ_in_scc));

        // simulation
        if (decomp_options_.iw_sim)
            prune_track(succ_in_scc);

        complement_mstate succ_state(scc_info_);
        succ_state.iw_sccs_.push_back(succ_in_scc);
        succ.push_back({succ_state, false});

        return succ;
    }

    std::vector<std::pair<complement_mstate, bool>> mh_compl::get_succ_track_to_active(complement_mstate mstate, bdd symbol)
    {
        std::vector<std::pair<complement_mstate, bool>> succ;

        std::set<unsigned> all_succ = this->get_all_successors(mstate.curr_reachable_, symbol);
        std::set<unsigned> succ_in_scc;
        std::set<unsigned> scc_states;
        for (auto i : scc_index_)
        {
            scc_states.insert(scc_info_.states_of(i).begin(), scc_info_.states_of(i).end());
        }
        std::set_intersection(all_succ.begin(), all_succ.end(), scc_states.begin(), scc_states.end(), std::inserter(succ_in_scc, succ_in_scc.begin()));

        std::set<unsigned> new_break_set;
        for (auto i : scc_index_)
        {
            std::set<unsigned> states_in_scc(scc_info_.states_of(i).begin(), scc_info_.states_of(i).end());
            std::set<unsigned> inter2;
            std::set_intersection(succ_in_scc.begin(), succ_in_scc.end(), states_in_scc.begin(), states_in_scc.end(), std::inserter(inter2, inter2.begin()));
            new_break_set.insert(inter2.begin(), inter2.end());
        }

        // simulation
        if (decomp_options_.iw_sim)
            prune_active(succ_in_scc, new_break_set);

        complement_mstate succ_state(scc_info_);
        succ_state.iw_sccs_.push_back(std::vector<unsigned>(succ_in_scc.begin(), succ_in_scc.end()));
        succ_state.iw_break_set_ = std::vector<unsigned>(new_break_set.begin(), new_break_set.end());
        succ.push_back({succ_state, false});

        return succ;
    }

    std::vector<std::pair<complement_mstate, bool>> mh_compl::get_succ_active(complement_mstate mstate, bdd symbol)
    {
        std::vector<std::pair<complement_mstate, bool>> succ;

        if (mstate.iw_break_set_.empty())
        {
            // empty break set -> return TT and switch to other scc
            if ((decomp_options_.merge_iwa or mstate.iw_sccs_.size() == 1) and mstate.acc_detsccs_.size() == 0 and mstate.na_sccs_.size() == 0)
            {
                auto succ = get_succ_track_to_active(mstate, symbol);
                for (auto &state : succ)
                {
                    state.second = true;
                }
                return succ;
            }
            else
            {
                auto succ = get_succ_track(mstate, symbol);
                for (auto &state : succ)
                {
                    state.second = true;
                }
                return succ;
            }
        }
        else
        {
            // return AT and stay in the same scc
            std::set<unsigned> all_succ = this->get_all_successors(mstate.curr_reachable_, symbol);
            std::set<unsigned> succ_in_scc;
            std::set<unsigned> scc_states;
            for (auto i : scc_index_)
            {
                scc_states.insert(scc_info_.states_of(i).begin(), scc_info_.states_of(i).end());
            }
            std::set_intersection(all_succ.begin(), all_succ.end(), scc_states.begin(), scc_states.end(), std::inserter(succ_in_scc, succ_in_scc.begin()));

            std::set<unsigned> new_break_set;
            std::set<unsigned> states_in_sccs;
            std::set<unsigned> break_set_succ = this->get_all_successors(std::vector<unsigned>(mstate.iw_break_set_.begin(), mstate.iw_break_set_.end()), symbol);
            for (auto i : scc_index_)
            {
                states_in_sccs.insert(scc_info_.states_of(i).begin(), scc_info_.states_of(i).end());
            }

            std::set<unsigned> inter;
            std::set_intersection(break_set_succ.begin(), break_set_succ.end(), states_in_sccs.begin(), states_in_sccs.end(), std::inserter(new_break_set, new_break_set.begin()));

            // simulation
            if (decomp_options_.iw_sim)
                prune_active(succ_in_scc, new_break_set);

            // no state with empty break set (tba)
            if (new_break_set.size() == 0)
            {
                // empty break set -> return TT and switch to other scc
                if ((decomp_options_.merge_iwa or mstate.iw_sccs_.size() == 1) and mstate.acc_detsccs_.size() == 0 and mstate.na_sccs_.size() == 0)
                {
                    auto succ = get_succ_track_to_active(mstate, symbol);
                    for (auto &state : succ)
                    {
                        state.second = true;
                    }
                    return succ;
                }
                else
                {
                    auto succ = get_succ_track(mstate, symbol);
                    for (auto &state : succ)
                    {
                        state.second = true;
                    }
                    return succ;
                }
            }
            else
            {
                // stay in active scc
                complement_mstate succ_state(scc_info_);
                succ_state.iw_sccs_.push_back(std::vector<unsigned>(succ_in_scc.begin(), succ_in_scc.end()));
                succ_state.iw_break_set_ = std::vector<unsigned>(new_break_set.begin(), new_break_set.end());
                succ_state.active_index_ = scc_index_[0];
                succ.push_back({succ_state, false});
            }
        }

        return succ;
    }
}