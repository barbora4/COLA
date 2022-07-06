#include "rank_compl.hpp"

#include <set>
#include <algorithm>

namespace cola
{
    std::set<unsigned> rank_complement::get_all_successors(std::set<unsigned> current_states, bdd symbol)
    {
        std::set<unsigned> successors;

        for (unsigned s : current_states)
        {
        for (const auto &t : aut_->out(s))
        {
            if (!bdd_implies(symbol, t.cond))
            continue;

            successors.insert(t.dst);
        }
        }

        return successors;
    }
    
    std::set<int> rank_complement::get_successors_with_box(std::set<unsigned> reachable, rank_state state, bdd letter)
    {
        std::set<int> succ;
        
        if (state.reachable.size() == 0)
            return succ;
        unsigned scc_index = scc_info_.scc_of(*(state.reachable.begin()));

        // box
        if (state.reachable.find(-1) != state.reachable.end())
        {
            auto all_succ = get_all_successors(reachable, letter);
            for (auto state : all_succ)
            {
                if (std::find(scc_info_.states_of(scc_index).begin(), scc_info_.states_of(scc_index).end(), state) != scc_info_.states_of(scc_index).end())
                    succ.insert((int)state);
            }

            for (auto state : all_succ)
            {
                if (predecessors_[scc_index].find(scc_info_.scc_of(state)) != predecessors_[scc_index].end())
                {
                    succ.insert(-1);
                    break;
                }
            }
        }

        // no box
        else 
        {
            std::set<unsigned> tmp;
            for (auto s : state.reachable)
            {
                tmp.insert((int)s);
            }
            auto all_succ = get_all_successors(tmp, letter);
            for (auto state : all_succ)
            {
                if (std::find(scc_info_.states_of(scc_index).begin(), scc_info_.states_of(scc_index).end(), state) != scc_info_.states_of(scc_index).end())
                    succ.insert((int)state);
            }
        }

        return succ;
    }

    void rank_complement::get_max(std::vector<ranking>& rankings)
    {
        std::vector<ranking> tmp(rankings.begin(), rankings.end());

        for (auto r : rankings)
        {
            unsigned max_rank = r.get_max_rank();
            if (std::any_of(tmp.begin(), tmp.end(), [max_rank, r](ranking r2){
                return r != r2 and max_rank == r2.get_max_rank() and r2.is_bigger(r);
            }))
            {
                // there is bigger ranking
                tmp.erase(std::remove(tmp.begin(), tmp.end(), r), tmp.end());
            }
        }

        rankings = tmp;
    }
    
    std::vector<std::pair<rank_state, bool>> rank_complement::get_succ_track(std::set<unsigned> reachable, rank_state state, bdd letter)
    {
        if (state.track)
        {
            // tracking type
            // TODO
        }
        
        else 
        {
            // active type
            // TODO
        }
    }

    std::vector<std::pair<rank_state, bool>> rank_complement::get_succ_track_to_active(std::set<unsigned> reachable, rank_state state, bdd letter)
    {
        std::vector<std::pair<rank_state, bool>> succ;
        
        if (state.track)
        {
            // tracking type
            rank_state new_state;
            auto reach_succ = get_successors_with_box(reachable, state, letter);
            new_state.reachable = reach_succ;
            succ.push_back({new_state, false});

            // TODO
            std::vector<std::tuple<int, int, bool>> r;
            for (auto s : reach_succ)
            {
                bool accepting = (s != -1 ? is_accepting_[s] : false);
                r.push_back(std::make_tuple(s, 2*reach_succ.size(), accepting));
            }
            std::vector<ranking> rankings = get_tight_rankings(r);
            get_max(rankings);
            for (auto r : rankings)
            {
                rank_state new_state;
                new_state.f = r;
                new_state.i = 0;
                for (auto pr : r)
                {
                    if (pr.second == 0)
                        new_state.O.insert(pr.first);
                }
                succ.push_back({new_state, false});
            }
        }

        else
        {
            // active type
            // TODO
        }

        for (auto state : succ)
        {
            std::cerr << get_set_string_box(state.first.reachable) << std::endl;
        }

        return succ;
    }

    std::vector<std::pair<rank_state, bool>> rank_complement::get_succ_active(std::set<unsigned> reachable, rank_state state, bdd letter, bool one_scc)
    {
        std::vector<std::pair<rank_state, bool>> succ;

        if (state.track)
        {
            // tracking type
            if (state.reachable.size() == 0 or (state.reachable.size() == 1 and state.reachable.find(-1) != state.reachable.end()))
            {
                // switch to other scc + acc trans
                if (one_scc)
                {
                    // only one scc
                    // TODO getSuccTrackToActive + acc trans
                }
                else
                {
                    // TODO getSuccTrack + acc trans
                }
            }

            else 
            {
                // TODO getSuccTrackToActive
                std::cerr << "getSuccTrackToActive" << std::endl;
                succ = get_succ_track_to_active(reachable, state, letter);
            }
        }

        else 
        {
            // active type
            // TODO
        }

        return succ;
    }
}