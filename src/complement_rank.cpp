#include "complement_rank.hpp"

namespace cola
{
    std::set<int>
    rank_comp::get_successors_with_box(std::set<unsigned> reachable, rank_state state, bdd letter, unsigned scc_index)
    {
        std::set<int> succ;

        std::set<int> state_set;
        if (state.track)
            state_set = state.reachable;
        else
        {
            for (auto pr : state.f)
                state_set.insert(pr.first);
        }

        if (state_set.size() == 0)
            return succ;

        // box
        if (state_set.find(-1) != state_set.end())
        {
            auto all_succ = get_all_successors(std::vector<unsigned>(reachable.begin(), reachable.end()), letter);
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
            for (auto s : state_set)
            {
                tmp.insert((unsigned)s);
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

    void
    rank_comp::get_max(std::vector<ranking> &rankings)
    {
        std::vector<ranking> tmp(rankings.begin(), rankings.end());

        for (auto r : rankings)
        {
            unsigned max_rank = r.get_max_rank();
            if (std::any_of(tmp.begin(), tmp.end(), [max_rank, r](ranking r2)
                            { return r != r2 and max_rank == r2.get_max_rank() and r2.is_bigger(r); }))
            {
                // there is bigger ranking
                tmp.erase(std::remove(tmp.begin(), tmp.end(), r), tmp.end());
            }
        }

        rankings = tmp;
    }

    std::vector<ranking>
    rank_comp::get_succ_rankings(ranking r, std::vector<std::tuple<int, int, bool>> restr, std::set<unsigned> reachable, bdd letter, unsigned scc_index)
    {
        std::vector<ranking> rankings = get_tight_rankings(restr);
        std::vector<ranking> tmp(rankings.begin(), rankings.end());
        std::set<ranking> rankings_set(rankings.begin(), rankings.end());

        for (ranking r2 : tmp)
        {
            if (r2.get_max_rank() != r.get_max_rank())
            {
                rankings_set.erase(r2);
                continue;
            }

            bool skip = false;
            for (auto pr : r)
            {
                auto state = pr.first;
                std::set<int> sing;
                sing.insert(state);
                rank_state tmp;
                tmp.reachable = sing;
                tmp.track = true;
                std::set<int> succ = get_successors_with_box(reachable, tmp, letter, scc_index);

                for (auto s : succ)
                {
                    if (r2[s] > r[state])
                    {
                        skip = true;
                        rankings_set.erase(r2);
                        break;
                    }
                }
                if (skip)
                    break;

                if (state != -1)
                {
                    std::set<unsigned> sing;
                    sing.insert((unsigned)state);
                    std::set<int> succ = get_all_successors_acc(sing, letter, scc_index);

                    unsigned rank = (r[state] % 2 == 0 ? r[state] : r[state] - 1);
                    for (auto s : succ)
                    {
                        if (r2[s] > rank)
                        {
                            skip = true;
                            rankings_set.erase(r2);
                            break;
                        }
                    }
                    if (skip)
                        break;
                }
            }
        }

        return std::vector<ranking>(rankings_set.begin(), rankings_set.end());
    }

    std::vector<ranking>
    rank_comp::get_maxrank(std::set<unsigned> reachable, rank_state state, bdd letter, unsigned scc_index)
    {
        std::vector<std::tuple<int, int, bool>> restr;

        std::set<int> domain;
        for (auto pr : state.f)
            domain.insert(pr.first);

        auto succ_domain = get_successors_with_box(reachable, state, letter, scc_index);
        int size = scc_info_.states_of(scc_index).size();
        for (auto s : succ_domain)
        {
            restr.push_back(std::make_tuple(s, 2 * (size + 1), (s == -1) ? false : is_accepting_[s]));
        }

        std::vector<ranking> succ_rankings = get_succ_rankings(state.f, restr, reachable, letter, scc_index);
        get_max(succ_rankings);
        return succ_rankings;
    }

    std::vector<std::pair<complement_mstate, bool>>
    rank_comp::get_succ_track()
    {
        std::vector<std::pair<complement_mstate, bool>> succ;

        rank_state state = mstate_.na_sccs_[true_index_];
        std::set<unsigned> reachable = std::set<unsigned>(mstate_.curr_reachable_.begin(), mstate_.curr_reachable_.end());

        if (state.track)
        {
            // tracking type
            rank_state new_state;
            new_state.track = true;
            new_state.reachable = get_successors_with_box(reachable, state, symbol_, scc_index_[0]);
            complement_mstate tmp(scc_info_);
            tmp.na_sccs_.push_back(new_state);
            succ.push_back({tmp, false});
        }

        else
        {
            // active type
            std::vector<ranking> maxrank = get_maxrank(reachable, state, symbol_, scc_index_[0]);
            if (maxrank.size() > 0)
            {
                rank_state s;
                s.f = maxrank[0];
                s.track = false;
                complement_mstate tmp(scc_info_);
                tmp.na_sccs_.push_back(s);
                succ.push_back({tmp, false});
            }
        }

        return succ;
    }

    std::vector<std::pair<complement_mstate, bool>>
    rank_comp::get_succ_track_to_active()
    {
        std::vector<std::pair<complement_mstate, bool>> succ;

        rank_state state = mstate_.na_sccs_[true_index_];
        std::set<unsigned> reachable = std::set<unsigned>(mstate_.curr_reachable_.begin(), mstate_.curr_reachable_.end());

        if (state.track)
        {
            // tracking type
            rank_state new_state;
            auto reach_succ = get_successors_with_box(reachable, state, symbol_, scc_index_[0]);
            new_state.reachable = reach_succ;
            complement_mstate tmp(scc_info_);
            tmp.na_sccs_.push_back(new_state);
            succ.push_back({tmp, false});

            std::vector<std::tuple<int, int, bool>> r;
            int size = scc_info_.states_of(scc_index_[0]).size();
            for (auto s : reach_succ)
            {
                bool accepting = (s != -1 ? is_accepting_[s] : false);
                r.push_back(std::make_tuple(s, 2 * (size + 1), accepting));
            }
            std::vector<ranking> rankings = get_tight_rankings(r);
            get_max(rankings);
            for (auto r : rankings)
            {
                rank_state new_state;
                new_state.f = r;
                new_state.i = 0;
                new_state.track = false;
                for (auto pr : r)
                {
                    if (pr.second == 0)
                        new_state.O.insert(pr.first);
                }
                complement_mstate tmp(scc_info_);
                tmp.na_sccs_.push_back(new_state);
                succ.push_back({tmp, false});
            }
        }

        else
        {
            // active type
            std::vector<std::pair<complement_mstate, bool>> ret = get_succ_track();

            if (ret.size() > 0)
            {
                rank_state new_state;
                new_state.track = false;
                new_state.f = ret[0].first.na_sccs_[true_index_].f;
                new_state.i = 0;
                for (auto pr : new_state.f)
                {
                    if (pr.second == 0)
                        new_state.O.insert(pr.first);
                }
                complement_mstate tmp(scc_info_);
                tmp.na_sccs_.push_back(new_state);
                succ.push_back({tmp, false});
            }

            // otherwise no successor
        }

        return succ;
    }

    std::vector<std::pair<complement_mstate, bool>>
    rank_comp::get_succ_active()
    {
        std::vector<std::pair<complement_mstate, bool>> succ;

        rank_state state = mstate_.na_sccs_[true_index_];
        std::set<unsigned> reachable = std::set<unsigned>(mstate_.curr_reachable_.begin(), mstate_.curr_reachable_.end());
        bool one_scc = ((mstate_.iw_sccs_.size() == 0) and (mstate_.acc_detsccs_.size() == 0) and (mstate_.na_sccs_.size() == 1));

        if (state.track)
        {
            // tracking type
            if (state.reachable.size() == 0 or (state.reachable.size() == 1 and state.reachable.find(-1) != state.reachable.end()))
            {
                // switch to other scc + acc trans
                if (one_scc)
                {
                    // only one scc
                    std::vector<std::pair<complement_mstate, bool>> ret = get_succ_track_to_active();

                    if (ret.size() > 0)
                        succ.push_back({ret[0].first, true});
                }
                else
                {
                    rank_state new_state;
                    new_state.track = true;
                    new_state.reachable = get_successors_with_box(reachable, state, symbol_, scc_index_[0]);
                    complement_mstate tmp(scc_info_);
                    tmp.na_sccs_.push_back(new_state);
                    succ.push_back({tmp, true});
                }
            }

            else
            {
                succ = get_succ_track_to_active();
            }
        }

        else
        {
            // active type
            std::vector<std::pair<complement_mstate, bool>> ret = get_succ_track();
            if (ret.size() > 0)
            {
                ranking g = ret[0].first.na_sccs_[true_index_].f;

                std::vector<rank_state> eta_3;
                std::vector<rank_state> eta_4;

                // eta 3
                if (state.O.size() > 0)
                {
                    rank_state new_state;
                    new_state.track = false;
                    new_state.f = g;
                    new_state.i = state.i;

                    rank_state tmp;
                    tmp.reachable = state.O;
                    std::set<int> O_succ = get_successors_with_box(reachable, tmp, symbol_, scc_index_[0]);
                    std::set<int> g_rev;
                    for (auto pr : g)
                    {
                        if (pr.second == state.i)
                            g_rev.insert(pr.first);
                    }
                    std::set_intersection(O_succ.begin(), O_succ.end(), g_rev.begin(), g_rev.end(), std::inserter(new_state.O, new_state.O.begin()));

                    eta_3.push_back(new_state);
                }
                else
                {
                    rank_state new_state;
                    new_state.track = false;
                    new_state.f = g;
                    new_state.i = (state.i + 2) % (g.get_max_rank() + 1);

                    std::set<int> dom_succ = get_successors_with_box(reachable, state, symbol_, scc_index_[0]);
                    std::set<int> g_rev;
                    for (auto pr : g)
                    {
                        if (pr.second == new_state.i)
                            g_rev.insert(pr.first);
                    }
                    std::set_intersection(dom_succ.begin(), dom_succ.end(), g_rev.begin(), g_rev.end(), std::inserter(new_state.O, new_state.O.begin()));

                    eta_3.push_back(new_state);
                }

                // eta 4
                if (state.i != state.f.get_max_rank() - 1)
                {
                    rank_state new_state;

                    std::set<int> M;
                    rank_state tmp;
                    tmp.reachable = state.O;
                    std::set<int> O_succ = get_successors_with_box(reachable, tmp, symbol_, scc_index_[0]);
                    std::set<int> g_rev;
                    for (auto pr : g)
                    {
                        if (pr.second == state.i)
                            g_rev.insert(pr.first);
                    }
                    std::set_intersection(O_succ.begin(), O_succ.end(), g_rev.begin(), g_rev.end(), std::inserter(M, M.begin()));

                    new_state.track = false;
                    new_state.i = state.i;
                    for (auto s : M)
                    {
                        if (s != -1 and is_accepting_[s])
                            new_state.O.insert(s);
                    }

                    ranking g_prime = g;
                    for (auto &pr : g_prime)
                    {
                        if (pr.first != -1 and (not is_accepting_[pr.first] and M.find(pr.first) != M.end()))
                            pr.second = pr.second - 1;
                    }
                    new_state.f = g_prime;

                    eta_4.push_back(new_state);
                }

                std::vector<rank_state> U;
                if (eta_3.size() > 0 and eta_3[0].O.size() == 0 and eta_3[0].i == eta_3[0].f.get_max_rank() - 1)
                    U.push_back(eta_3[0]);
                if (eta_4.size() > 0 and eta_4[0].O.size() == 0 and eta_4[0].i == eta_4[0].f.get_max_rank() - 1)
                    U.push_back(eta_4[0]);

                if (eta_3.size() > 0 and std::find(U.begin(), U.end(), eta_3[0]) == U.end())
                {
                    complement_mstate tmp_mstate(scc_info_);
                    tmp_mstate.na_sccs_.push_back(eta_3[0]);
                    succ.push_back({tmp_mstate, false});
                }
                if (eta_4.size() > 0 and std::find(U.begin(), U.end(), eta_4[0]) == U.end())
                {
                    complement_mstate tmp_mstate(scc_info_);
                    tmp_mstate.na_sccs_.push_back(eta_4[0]);
                    succ.push_back({tmp_mstate, false}); 
                }

                for (rank_state s : U)
                {
                    if (not one_scc)
                    {
                        rank_state tmp;
                        tmp.track = false;
                        tmp.f = s.f;
                        complement_mstate tmp_mstate(scc_info_);
                        tmp_mstate.na_sccs_.push_back(tmp);
                        succ.push_back({tmp_mstate, true});
                    }
                    else
                    {
                        complement_mstate tmp(scc_info_);
                        tmp.na_sccs_.push_back(s);
                        succ.push_back({tmp, true});
                    }
                }
            }
        }

        return succ;
    }
}