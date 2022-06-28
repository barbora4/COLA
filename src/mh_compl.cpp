#include "mh_compl.hpp"

#include <stack>

std::vector<std::set<int>> 
  mh_complement::get_reachable_vector()
  {
    std::vector<std::set<int>> list_set(aut_->num_states());
    std::vector<std::vector<int>> list_vector(aut_->num_states());

    for (unsigned s = 0; s < aut_->num_states(); s++)
    {
      reachable_vector_.push_back(std::set<int>());
      list_set[s] = std::set<int>();

      // iterate over all transitions from s
      for (const auto &t : aut_->out(s))
      {
        list_set[s].insert(t.dst);
      }

      list_vector[s] = std::vector<int>(list_set[s].begin(), list_set[s].end());
    }

    for (int s = 0; s < aut_->num_states(); s++)
    {
      std::set<int> tmp({s});
      reachable_vector_[s] = reachable_vertices(list_vector, tmp);
    }

    return reachable_vector_;
  }

  std::set<int>
  mh_complement::reachable_vertices(std::vector<std::vector<int>> list, std::set<int> from)
  {
    std::set<int> all(from);
    std::stack<int> stack;
    int item;
    for (int it : all)
      stack.push(it);

    while (stack.size() > 0)
    {
      item = stack.top();
      stack.pop();
      for (int dst : list[item])
      {
        if (all.find(dst) == all.end())
        {
          stack.push(dst);
          all.insert(dst);
        }
      }
    }
    return all;
  }

  std::vector<std::pair<std::set<unsigned>, unsigned>>
  mh_complement::get_succ_track(std::set<unsigned> reachable, std::set<unsigned> reach_in_scc, bdd symbol, std::vector<unsigned> scc_index)
  {
    std::vector<std::pair<std::set<unsigned>, unsigned>> succ;

    std::set<unsigned> all_succ = this->get_all_successors(reachable, symbol);
    std::vector<unsigned> succ_in_scc;
    std::set<unsigned> scc_states;
    for (auto i : scc_index)
    {
      scc_states.insert(scc_info_.states_of(i).begin(), scc_info_.states_of(i).end());
    }
    std::set_intersection(all_succ.begin(), all_succ.end(), scc_states.begin(), scc_states.end(), std::back_inserter(succ_in_scc));

    // simulation
    if (decomp_options_.iw_sim)
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

    succ.push_back({std::set<unsigned>(succ_in_scc.begin(), succ_in_scc.end()), 0});

    return succ;
  }

  std::vector<std::pair<std::pair<std::set<unsigned>, std::set<unsigned>>, unsigned>>
  mh_complement::get_succ_track_to_active(std::set<unsigned> reachable, std::set<unsigned> reach_in_scc, bdd symbol, std::vector<unsigned> scc_index)
  { 
    std::vector<std::pair<std::pair<std::set<unsigned>, std::set<unsigned>>, unsigned>> succ;

    std::set<unsigned> all_succ = this->get_all_successors(reachable, symbol);
    std::set<unsigned> succ_in_scc;
    std::set<unsigned> scc_states;
    for (auto i : scc_index)
    {
      scc_states.insert(scc_info_.states_of(i).begin(), scc_info_.states_of(i).end());
    }
    std::set_intersection(all_succ.begin(), all_succ.end(), scc_states.begin(), scc_states.end(), std::inserter(succ_in_scc, succ_in_scc.begin()));

    std::set<unsigned> new_break_set;
    for (auto i : scc_index)
    {
      std::set<unsigned> states_in_scc(scc_info_.states_of(i).begin(), scc_info_.states_of(i).end());
      if (cola::is_accepting_weakscc(scc_types_, i))
      {
        std::set<unsigned> inter2;
        std::set_intersection(succ_in_scc.begin(), succ_in_scc.end(), states_in_scc.begin(), states_in_scc.end(), std::inserter(inter2, inter2.begin()));
        new_break_set.insert(inter2.begin(), inter2.end());
      }
    }

    // simulation
    if (decomp_options_.iw_sim)
    {
      // remove smaller states from S
      std::set<unsigned> new_S(succ_in_scc.begin(), succ_in_scc.end());
      for (auto pr : dir_sim_)
      {
        // std::cerr << pr.first << " " << pr.second << std::endl;
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
      // std::cerr << std::endl;
      succ_in_scc = new_S;
    }

    succ.push_back({{std::set<unsigned>(succ_in_scc.begin(), succ_in_scc.end()), new_break_set}, 0});

    return succ;
  }

  std::pair<std::vector<std::pair<std::set<unsigned>, unsigned>>, std::vector<std::pair<std::pair<std::set<unsigned>, std::set<unsigned>>, unsigned>>>
  mh_complement::get_succ_active(std::set<unsigned> reachable, std::set<unsigned> reach_in_scc, bdd symbol, std::vector<unsigned> scc_index, std::vector<unsigned> break_set, bool one_scc)
  {
    std::vector<std::pair<std::set<unsigned>, unsigned>> succ_tt;
    std::vector<std::pair<std::pair<std::set<unsigned>, std::set<unsigned>>, unsigned>> succ_at;

    if (break_set.empty())
    {
      // empty break set -> return TT and switch to other scc
      succ_tt = get_succ_track(reachable, reach_in_scc, symbol, scc_index);
      for (auto &succ : succ_tt)
      {
        succ.second = 1;
      }
    }
    else
    {
      // return AT and stay in the same scc
      std::set<unsigned> all_succ = this->get_all_successors(reachable, symbol);
      std::set<unsigned> succ_in_scc;
      std::set<unsigned> scc_states;
      for (auto i : scc_index)
      {
        scc_states.insert(scc_info_.states_of(i).begin(), scc_info_.states_of(i).end());
      }
      std::set_intersection(all_succ.begin(), all_succ.end(), scc_states.begin(), scc_states.end(), std::inserter(succ_in_scc, succ_in_scc.begin()));

      std::set<unsigned> new_break_set;
      std::set<unsigned> states_in_sccs;
      std::set<unsigned> break_set_succ = this->get_all_successors(std::set<unsigned>(break_set.begin(), break_set.end()), symbol);
      for (auto i : scc_index)
      {
        states_in_sccs.insert(scc_info_.states_of(i).begin(), scc_info_.states_of(i).end()); 
      }
        
      std::set<unsigned> inter;
      std::set_intersection(break_set_succ.begin(), break_set_succ.end(), states_in_sccs.begin(), states_in_sccs.end(), std::inserter(new_break_set, new_break_set.begin()));

      // simulation
      if (decomp_options_.iw_sim)
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


      // no state with empty break set (tba) 
      if (new_break_set.size() == 0)
      {
        // empty break set -> return TT and switch to other scc
        if (one_scc)
        {
          succ_at = get_succ_track_to_active(reachable, reach_in_scc, symbol, scc_index);
          for (auto &succ : succ_at)
          {
            succ.second = 1;
          }
        }
        else
        {
          succ_tt = get_succ_track(reachable, reach_in_scc, symbol, scc_index);
          for (auto &succ : succ_tt)
          {
            succ.second = 1;
          }
        }
      }
      else
      {
        succ_at.push_back({{std::set<unsigned>(succ_in_scc.begin(), succ_in_scc.end()), new_break_set}, 0});
      }
    }

    return {succ_tt, succ_at};
  }

  std::set<unsigned>
  mh_complement::get_all_successors(std::set<unsigned> current_states, bdd symbol)
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