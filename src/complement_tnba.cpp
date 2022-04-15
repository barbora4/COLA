// Copyright (C) 2017-2019 Laboratoire de Recherche et Développement
// de l'Epita.
// Copyright (C) 2022  The COLA Authors
//
// COLA is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// COLA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

//#include "optimizer.hpp"
#include "cola.hpp"
#include "simulation.hpp"
#include "types.hpp"
//#include "struct.hpp"

#include <deque>
#include <map>
#include <set>

#include <spot/misc/hashfunc.hh>
#include <spot/twaalgos/isdet.hh>
#include <spot/twaalgos/sccinfo.hh>
#include <spot/twaalgos/isunamb.hh>
#include <spot/twa/acc.hh>
#include <spot/twaalgos/degen.hh>
#include <spot/twaalgos/simulation.hh>
#include <spot/twaalgos/determinize.hh>
#include <spot/twaalgos/parity.hh>
#include <spot/twaalgos/cleanacc.hh>
#include <spot/twaalgos/postproc.hh>
#include <spot/misc/bddlt.hh>
#include <spot/parseaut/public.hh>
#include <spot/twaalgos/complement.hh>
#include <spot/twaalgos/hoa.hh>
#include <spot/misc/version.hh>
#include <spot/twa/acc.hh>

// Complementation of Buchi automara based on SCC decomposition
// We classify three types of SCCs in the input NBA:
// 1. inherently weak SCCs (IWCs): every cycle in the SCC will not visit accepting transitions or every cycle visits an accepting transition
// 2. deterministic accepting SCCs (DACs): states in the SCC have at most one successor remain in the same SCC for a letter
// 3. nondeterministic accepting SCCs (NACs): has an accepting transition and nondeterministic

namespace cola
{
  // (C, S, B) for complementing DACs
  const int NCSB_C = 2;
  const int NCSB_S = 4;
  const int NCSB_B = 3;

  enum slice_mark
  {
    None = -1,
    Inf = 0,
    New = 1,
    Fin = 2
  };

  class complement_mstate
  {
  public:
    // the number of states num, default values, and number of NACs
    complement_mstate(spot::scc_info &si, unsigned num_det_sccs)
        : si_(si)
    {
      for (unsigned i = 0; i < num_det_sccs; i++)
      {
        detscc_ranks_.emplace_back(std::vector<state_rank>());
      }
    }

    complement_mstate(const complement_mstate &other)
        : si_(other.si_)
    {
      this->break_set_.clear();
      this->break_set_.insert(other.break_set_.begin(), other.break_set_.end());
      this->weak_set_.clear();
      this->weak_set_.insert(other.weak_set_.begin(), other.weak_set_.end());

      this->detscc_ranks_.clear();
      for (unsigned i = 0; i < other.detscc_ranks_.size(); i++)
      {
        std::vector<state_rank> copy = other.detscc_ranks_[i];
        this->detscc_ranks_.emplace_back(copy);
      }

      this->nondetscc_ranks_.clear();
      this->nondetscc_marks_.clear();
      for (unsigned i = 0; i < other.nondetscc_ranks_.size(); i++)
      {
        std::set<unsigned> copy = other.nondetscc_ranks_[i];
        this->nondetscc_ranks_.emplace_back(copy);
        this->nondetscc_marks_.push_back(other.nondetscc_marks_[i]);
      }

      this->set_iw_sccs(other.iw_sccs_);
      this->set_iw_break_set(other.iw_break_set_);
      this->det_break_set_ = other.det_break_set_;
      this->set_active_index(other.active_index_);
      this->set_acc_detsccs(other.acc_detsccs_);
      this->curr_reachable_ = other.curr_reachable_;
    }

    std::set<unsigned>
    get_reach_set() const;

    std::set<unsigned>
    get_weak_set() const;

    bool is_empty() const;

    int get_max_rank() const;

    bool operator<(const complement_mstate &other) const;
    bool operator==(const complement_mstate &other) const;

    complement_mstate &
    operator=(const complement_mstate &other)
    {
      this->si_ = other.si_;
      this->break_set_.clear();
      this->break_set_.insert(other.break_set_.begin(), other.break_set_.end());
      this->weak_set_.clear();
      this->weak_set_.insert(other.weak_set_.begin(), other.weak_set_.end());

      this->detscc_ranks_.clear();
      for (unsigned i = 0; i < other.detscc_ranks_.size(); i++)
      {
        std::vector<state_rank> copy = other.detscc_ranks_[i];
        this->detscc_ranks_.emplace_back(copy);
      }

      this->nondetscc_ranks_.clear();
      this->nondetscc_marks_.clear();
      for (unsigned i = 0; i < other.nondetscc_ranks_.size(); i++)
      {
        std::set<unsigned> copy = other.nondetscc_ranks_[i];
        this->nondetscc_ranks_.emplace_back(copy);
        this->nondetscc_marks_.push_back(other.nondetscc_marks_[i]);
      }

      this->set_iw_sccs(other.iw_sccs_);
      this->set_iw_break_set(other.iw_break_set_);
      this->det_break_set_ = other.det_break_set_;
      this->set_active_index(other.active_index_);
      this->set_acc_detsccs(other.acc_detsccs_);
      this->curr_reachable_ = other.curr_reachable_;

      return *this;
    }

    size_t hash() const;

    // SCC information
    spot::scc_info &si_;
    // 1. NAC by slice-based complementation
    std::vector<std::set<unsigned>> nondetscc_ranks_;
    std::vector<slice_mark> nondetscc_marks_; // marks for in O or not?

    // 2. DAC by determinization or NCSB
    std::vector<std::vector<state_rank>> detscc_ranks_;

    // 3. IWC states point to RANK_WEAK
    // breakpoint construction for weak accepting SCCs
    std::set<unsigned> weak_set_;
    std::set<unsigned> break_set_;

    std::vector<std::vector<unsigned>> iw_sccs_;
    std::vector<std::pair<std::vector<unsigned>, std::vector<unsigned>>> acc_detsccs_;
    std::vector<unsigned> iw_break_set_;
    std::vector<unsigned> det_break_set_;
    int active_index_ = 0;
    std::vector<unsigned> curr_reachable_;

    void
    set_iw_sccs(std::vector<std::vector<unsigned>> iw_sccs)
    {
      this->iw_sccs_ = iw_sccs;
    }

    void
    set_acc_detsccs(std::vector<std::pair<std::vector<unsigned>, std::vector<unsigned>>> acc_detsccs)
    {
      this->acc_detsccs_ = acc_detsccs;
    }

    void
    set_iw_break_set(std::vector<unsigned> iw_break_set)
    {
      this->iw_break_set_ = iw_break_set;
    }

    void
    set_active_index(int index)
    {
      this->active_index_ = index;
    }
  };

  struct complement_mstate_hash
  {
    size_t
    operator()(const complement_mstate &s) const noexcept
    {
      return s.hash();
    }
  };

  bool
  complement_mstate::operator<(const complement_mstate &other) const
  {
    if (active_index_ == other.active_index_)
    {
      if (iw_sccs_ == other.iw_sccs_)
      {
        if (iw_break_set_ == other.iw_break_set_)
        {
          if (det_break_set_ == other.det_break_set_)
          {
            if (acc_detsccs_ == other.acc_detsccs_)
            {
              if (curr_reachable_ == other.curr_reachable_)
              {
                return false;
              }
              else
              {
                return curr_reachable_ < other.curr_reachable_;
              }
            }
            else
            {
              return acc_detsccs_ < other.acc_detsccs_;
            }
          }
          else
          {
            return det_break_set_ < other.det_break_set_;
          }
        }
        else
        {
          return iw_break_set_ < other.iw_break_set_;
        }
      }
      else
      {
        return iw_sccs_ < other.iw_sccs_;
      }
    }
    else
    {
      return active_index_ < other.active_index_;
    }
  }

  bool
  complement_mstate::operator==(const complement_mstate &other) const
  {
    if (this->active_index_ != other.active_index_)
      return false;
    if (this->iw_sccs_ != other.iw_sccs_)
      return false;
    if (this->iw_break_set_ != other.iw_break_set_)
      return false;
    if (this->det_break_set_ != other.det_break_set_)
      return false;
    if (this->acc_detsccs_ != other.acc_detsccs_)
      return false;
    if (this->curr_reachable_ != other.curr_reachable_)
      return false;
    return true;
  }

  int complement_mstate::get_max_rank() const
  {
    return -1;
  }

  std::set<unsigned>
  complement_mstate::get_reach_set() const
  {
    return std::set<unsigned>(curr_reachable_.begin(), curr_reachable_.end());
  }

  bool complement_mstate::is_empty() const
  {
    if (!weak_set_.empty())
    {
      return false;
    }
    for (unsigned i = 0; i < detscc_ranks_.size(); i++)
    {
      if (!detscc_ranks_[i].empty())
      {
        return false;
      }
    }

    if (!nondetscc_ranks_.empty())
    {
      return false;
    }

    return true;
  }

  std::set<unsigned>
  complement_mstate::get_weak_set() const
  {
    return weak_set_;
  }

  size_t
  complement_mstate::hash() const
  {
    size_t res = 0;
    res = (res << 3) ^ active_index_;
    for (auto v : iw_sccs_)
    {
      for (auto s : v)
      {
        res ^= (res << 3) ^ s;
      }
    }
    for (auto s : iw_break_set_)
    {
      res ^= (res << 3) ^ s;
    }
    for (auto s : det_break_set_)
    {
      res ^= (res << 3) ^ s;
    }
    for (auto pr : acc_detsccs_)
    {
      for (auto v : pr.first)
      {
        res ^= (res << 3) ^ v;
      }
      for (auto v : pr.second)
      {
        res ^= (res << 3) ^ v;
      }
    }
    for (auto s : curr_reachable_)
    {
      res ^= (res << 3) ^ s;
    }

    return res;
  }

  // computation of deterministic successor
  class compute_det_succ
  {
    public:
    spot::scc_info &si_;
    // current ranking values of the DAC states
    const std::vector<state_rank> &curr_ranks_;
    // the reachable states at this level inside this SCC
    std::set<unsigned> &next_level_;
    // transitions
    std::unordered_map<unsigned, std::vector<std::pair<bool, unsigned>>> &det_trans_;
    // DAC number
    unsigned scc_;

    compute_det_succ(spot::scc_info &si, unsigned scc, const std::vector<state_rank> &curr_ranks, std::set<unsigned> &next_level, std::unordered_map<unsigned, std::vector<std::pair<bool, unsigned>>> &det_trans)
        : si_(si), scc_(scc), curr_ranks_(curr_ranks), next_level_(next_level), det_trans_(det_trans)
    {
    }

    std::vector<state_rank> next_ranks_;

    void
    compute()
    {
      next_ranks_.clear();
      // list of deterministic states, already ordered by its labelling
      std::map<unsigned, int> succ_nodes;
      int max_rnk = -1;
      // print_label_vec(acc_det_states);
      for (unsigned j = 0; j < curr_ranks_.size(); j++)
      {
        unsigned s = curr_ranks_[j].first;
        int curr_rnk = curr_ranks_[j].second;
        max_rnk = std::max(max_rnk, curr_rnk);
        assert(curr_rnk == j);
        // states and ranking
        for (const auto &t : det_trans_[s])
        {
          unsigned succ_scc = si_.scc_of(t.second);
          // ignore the states that go to other SCCs
          if (scc_ != succ_scc)
            continue;
          next_level_.erase(t.second);
          // Stay in the same accepting deterministic SCC or just enter this SCC
          // All DAC-states already have assigned with MAX_RANK
          auto it = succ_nodes.emplace(t.second, curr_rnk);
          if (!it.second) // already there
          {
            int prev_rnk = it.first->second;
            it.first->second = std::min(curr_rnk, prev_rnk);
          }
        }
      }
      ++max_rnk;
      // put them into succ
      for (unsigned p : next_level_)
      {
        // insertion failed is possible
        succ_nodes.emplace(p, max_rnk);
        ++max_rnk;
      }
      //succ.detscc_labels_[i].clear();
      for (auto &node : succ_nodes)
      {
        next_ranks_.emplace_back(node.first, node.second);
      }
    }

    int
    get_color()
    {
      int min_acc = -1;
      int min_dcc = -1;
      std::map<unsigned, int> succ_nodes;
      for (auto &p : next_ranks_)
      {
        succ_nodes[p.first] = p.second;
      }

      for (unsigned j = 0; j < curr_ranks_.size(); j++)
      {
        bool has_succ = false;
        bool has_acc = false;
        unsigned s = curr_ranks_[j].first;
        int curr_rnk = curr_ranks_[j].second;
        assert(curr_rnk == j);
        for (const auto &t : det_trans_[s])
        {
          // ignore the states that are not existing
          if (succ_nodes.find(t.second) == succ_nodes.end())
          {
            continue;
          }
          // 1. first they should be in the same SCC
          // 2. second the label should be equal
          if (si_.scc_of(s) == si_.scc_of(t.second) && succ_nodes[t.second] == curr_rnk)
          {
            has_succ = true;
            has_acc = has_acc || t.first;
          }
        }
        if (!has_succ)
        {
          // i. no successor, record the smaller label
          if (min_dcc == -1)
          {
            min_dcc = 2 * j + 1;
          }
        }
        else if (has_acc && min_acc == -1)
        {
          // ii. see an accepting transition
          min_acc = 2 * (j + 1);
        }
        // number
      }
      // reorganize the indices
      std::sort(next_ranks_.begin(), next_ranks_.end(), rank_compare);
      for (int k = 0; k < next_ranks_.size(); k++)
      {
        next_ranks_[k].second = k;
      }
      // compute the color
      return std::min(min_acc, min_dcc);
    }
  };


  // Miyano-Hayashi complementation for IW components
  class mh_complement
  {
  private:
    // source automaton
    const spot::const_twa_graph_ptr aut_;
    spot::scc_info scc_info_;
    std::string scc_types_;

  public:
    mh_complement(const spot::const_twa_graph_ptr& aut, spot::scc_info& scc_info, std::string scc_types) : aut_(aut), scc_info_(scc_info), scc_types_(scc_types){}

    std::vector<std::pair<std::set<unsigned>, unsigned>> get_succ_track(std::set<unsigned> reachable, std::set<unsigned> reach_in_scc, bdd symbol, std::vector<unsigned> scc_index);
    std::vector<std::pair<std::pair<std::set<unsigned>, std::set<unsigned>>, unsigned>> get_succ_track_to_active(std::set<unsigned> reachable, std::set<unsigned> reach_in_scc, bdd symbol, std::vector<unsigned> scc_index);
    std::pair<std::vector<std::pair<std::set<unsigned>, unsigned>>, std::vector<std::pair<std::pair<std::set<unsigned>, std::set<unsigned>>, unsigned>>> get_succ_active(std::set<unsigned> reachable, std::set<unsigned> reach_in_scc, bdd symbol, std::vector<unsigned> scc_index, std::vector<unsigned> break_set);

    std::set<unsigned> get_all_successors(std::set<unsigned> current_states, bdd symbol);
  };

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
      if (is_accepting_weakscc(scc_types_, i))
      {
        std::set<unsigned> inter2;
        std::set_intersection(succ_in_scc.begin(), succ_in_scc.end(), states_in_scc.begin(), states_in_scc.end(), std::inserter(inter2, inter2.begin()));
        new_break_set.insert(inter2.begin(), inter2.end());
      }
    }    

    succ.push_back({{std::set<unsigned>(succ_in_scc.begin(), succ_in_scc.end()), new_break_set}, 0});

    return succ;
  }

  std::pair<std::vector<std::pair<std::set<unsigned>, unsigned>>, std::vector<std::pair<std::pair<std::set<unsigned>, std::set<unsigned>>, unsigned>>>
  mh_complement::get_succ_active(std::set<unsigned> reachable, std::set<unsigned> reach_in_scc, bdd symbol, std::vector<unsigned> scc_index, std::vector<unsigned> break_set)
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

      // std::set<unsigned> break_set_succ = this->get_all_successors(std::set<unsigned>(break_set.begin(), break_set.end()), symbol);
      // std::set<unsigned> inter;
      // std::set_intersection(break_set_succ.begin(), break_set_succ.end(), succ_in_scc.begin(), succ_in_scc.end(), std::inserter(inter, inter.begin()));

      std::set<unsigned> new_break_set;
      for (auto i : scc_index)
      {
        std::set<unsigned> states_in_scc(scc_info_.states_of(i).begin(), scc_info_.states_of(i).end());
        std::set<unsigned> break_set_in_scc;
        std::set_intersection(states_in_scc.begin(), states_in_scc.end(), break_set.begin(), break_set.end(), std::inserter(break_set_in_scc, break_set_in_scc.begin()));
        std::set<unsigned> break_set_succ = this->get_all_successors(std::set<unsigned>(break_set_in_scc.begin(), break_set_in_scc.end()), symbol);
        std::set<unsigned> inter;
        std::set_intersection(break_set_succ.begin(), break_set_succ.end(), states_in_scc.begin(), states_in_scc.end(), std::inserter(inter, inter.begin()));
        if (is_accepting_weakscc(scc_types_, i))
        {
          std::set<unsigned> inter2;
          std::set_intersection(inter.begin(), inter.end(), states_in_scc.begin(), states_in_scc.end(), std::inserter(inter2, inter2.begin()));
          new_break_set.insert(inter2.begin(), inter2.end());
        }
      }
      
      succ_at.push_back({{std::set<unsigned>(succ_in_scc.begin(), succ_in_scc.end()), new_break_set}, 0});
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


  // complementation Buchi automata
  class tnba_complement
  {
  private:
    // The source automaton.
    const spot::const_twa_graph_ptr aut_;

    // SCCs information of the source automaton.
    spot::scc_info &si_;

    // Complement decomposition options
    compl_decomp_options decomp_options_;

    // Number of states in the input automaton.
    unsigned nb_states_;

    // state_simulator
    state_simulator simulator_;

    // delayed simulation
    delayed_simulation delayed_simulator_;

    // The parity automata being built.
    spot::twa_graph_ptr res_;

    // the number of indices
    unsigned sets_ = 0;

    unsigned num_colors_;

    spot::option_map &om_;

    // use ambiguous
    bool use_unambiguous_;

    bool use_scc_;

    // use stutter
    bool use_stutter_;

    bool use_simulation_;

    // Association between labelling states and state numbers of the
    // DPA.
    std::unordered_map<complement_mstate, unsigned, complement_mstate_hash> rank2n_;

    // States to process.
    std::deque<std::pair<complement_mstate, unsigned>> todo_;

    // Support for each state of the source automaton.
    std::vector<bdd> support_;

    // Propositions compatible with all transitions of a state.
    std::vector<bdd> compat_;

    // is accepting for states
    std::vector<bool> is_accepting_;

    // Whether a SCC is deterministic or not
    std::string scc_types_;

    // State names for graphviz display
    std::vector<std::string> *names_;

    // the index of each weak SCCs
    std::vector<unsigned> weaksccs_;
    // the index of each deterministic accepting SCCs
    std::vector<unsigned> acc_detsccs_;
    // the index of each deterministic accepting SCCs
    std::vector<unsigned> acc_nondetsccs_;

    // Show Rank states in state name to help debug
    bool show_names_;

    // maximal ranking in a labelling

    std::string
    get_det_string(const std::vector<state_rank> &states)
    {
      std::string res = "[";
      bool first_state = true;
      for (unsigned p = 0; p < states.size(); p++)
      {
        if (!first_state)
          res += " < ";
        first_state = false;
        res += std::to_string(states[p].first);
      }
      res += "]";
      return res;
    }

    std::string
    get_name(const complement_mstate &ms)
    {
      std::string name;
      name += "(";
      name += get_set_string(std::set<unsigned>(ms.curr_reachable_.begin(), ms.curr_reachable_.end()));
      name += ",";
      for (auto partial : ms.iw_sccs_)
      {
        const std::set<unsigned> tmp(partial.begin(), partial.end());
        name += get_set_string(tmp);
      }
      name += ",";
      for (auto partial : ms.acc_detsccs_)
      {
        const std::set<unsigned> tmp(partial.first.begin(), partial.first.end());
        name += get_set_string(tmp);
        name += "+";
        const std::set<unsigned> aux(partial.second.begin(), partial.second.end());
        name += get_set_string(aux);
      }
      name += ",";
      const std::set<unsigned> breakset(ms.iw_break_set_.begin(), ms.iw_break_set_.end());
      name += get_set_string(breakset);
      name += ",";
      const std::set<unsigned> det_breakset(ms.det_break_set_.begin(), ms.det_break_set_.end());
      name += get_set_string(det_breakset);
      name += ",";
      name += std::to_string(ms.active_index_);
      name += ")";
      return name;
    }

    // From a Rank state, looks for a duplicate in the map before
    // creating a new state if needed.
    unsigned
    new_state(complement_mstate &s)
    {
      complement_mstate dup(s);
      auto p = rank2n_.emplace(dup, 0);
      if (p.second) // This is a new state
      {
        p.first->second = res_->new_state();
        if (show_names_)
        {
          names_->push_back(get_name(p.first->first));
        }
        todo_.emplace_back(dup, p.first->second);
      }
      return p.first->second;
    }

    bool exists(complement_mstate &s)
    {
      return rank2n_.end() != rank2n_.find(s);
    }

    void remove_rank(std::vector<state_rank> &nodes, std::set<unsigned> &to_remove)
    {
      std::vector<state_rank> tmp;
      auto it1 = nodes.begin();
      while (it1 != nodes.end())
      {
        auto old_it1 = it1++;
        if (to_remove.find(old_it1->first) != to_remove.end())
        {
          it1 = nodes.erase(old_it1);
        }
      }
    }

    void merge_redundant_states(complement_mstate &ms, std::vector<state_rank> &nodes, bool nondet)
    {
      // auto it1 = nodes.begin();
      // while (it1 != nodes.end())
      // {
      //   auto old_it1 = it1++;
      //   for (auto it2 = nodes.begin(); it2 != nodes.end(); ++it2)
      //   {
      //     if (old_it1 == it2)
      //       continue;
      //     unsigned i = old_it1->first;
      //     unsigned j = it2->first;
      //     if (!(simulator_.simulate(j, i) || delayed_simulator_.simulate(j, i)))
      //     {
      //       continue;
      //     }
      //     int brace_i = old_it1->second;
      //     int brace_j = it2->second;
      //     bool remove = false;
      //     if (nondet)
      //     {
      //       // need to compare there nesting pattern
      //       unsigned scc_i = si_.scc_of(i);
      //       int scc_index = get_nondetscc_index(scc_i);
      //       std::vector<int> &braces = ms.nondetscc_breaces_[scc_index];
      //       // std::cout << "State " << i << " brace: " << brace_i << std::endl;
      //       // std::cout << "State " << j << " brace: " << brace_j << std::endl;
      //       //print_pattern_vec(braces, braces.size());
      //       if (compare_braces(braces, brace_j, brace_i))
      //       {
      //         remove = true;
      //       }
      //     }
      //     else if (brace_j < brace_i)
      //     {
      //       remove = true;
      //     }
      //     if (remove)
      //     {
      //       it1 = nodes.erase(old_it1);
      //       break;
      //     }
      //   }
      // }
    }

    // remove a state i if it is simulated by a state j
    void
    make_simulation_state(complement_mstate &ms, std::set<unsigned> &level_states, std::vector<std::vector<state_rank>> &det_succs, std::vector<std::vector<state_rank>> &nondet_succs)
    {
      std::set<unsigned> det_remove;
      std::set<unsigned> nondet_remove;
      for (unsigned i : level_states)
      {
        for (unsigned j : level_states)
        {
          if (i == j)
            continue;
          unsigned scc_i = si_.scc_of(i);
          // j simulates i and j cannot reach i
          if ((simulator_.simulate(j, i) || delayed_simulator_.simulate(j, i)) && simulator_.can_reach(j, i) == 0)
          {
            // std::cout << j << " simulated " << i << std::endl;
            // std::cout << "can_reach = " << simulator_.can_reach(j, i) << std::endl;
            if (is_weakscc(scc_types_, scc_i))
            {
              ms.weak_set_.erase(i);
              ms.break_set_.erase(i);
            }
            else if (is_accepting_detscc(scc_types_, scc_i))
            {
              det_remove.insert(i);
            }
            else if (is_accepting_nondetscc(scc_types_, scc_i))
            {
              nondet_remove.insert(i);
            }
          }
        }
      }
      for (std::vector<state_rank> &succ : det_succs)
      {
        remove_rank(succ, det_remove);
      }
      // for (unsigned i = 0; i < det_remove.size(); i++)
      // {
      //   remove_label(ms.detscc_labels_[i], det_remove[i]);
      //   // now remove more
      //   merge_redundant_states(ms, ms.detscc_labels_[i], false);
      // }
      // for (unsigned i = 0; i < nondet_remove.size(); i++)
      // {
      //   remove_label(ms.nondetscc_labels_[i], nondet_remove[i]);
      //   merge_redundant_states(ms, ms.nondetscc_labels_[i], true);
      // }
    }

    void
    csb_successors(const std::vector<state_rank> &curr_det_states, int scc_index, std::vector<int> &next_scc_indices, std::vector<std::map<unsigned, int>> &succ_maps, std::vector<bool> &acc_succs, std::set<unsigned> &next_detstates, std::unordered_map<unsigned, std::vector<std::pair<bool, unsigned>>> &det_cache, unsigned active_index)
    {
      // std::cout << "csb_successor scc " << scc_index << " rank: " << get_rank_string(curr_det_states) << std::endl;
      //1. Handle S states.
      // Treated first because we can escape early if the letter
      // leads to an accepting transition for a Safe state.
      // std::vector<std::map<unsigned, int>> succ_maps;
      std::map<unsigned, int> succ_nodes;
      for (unsigned i = 0; i < curr_det_states.size(); ++i)
      {
        // ignore other states
        if (curr_det_states[i].second != NCSB_S)
          continue;

        unsigned curr_s = curr_det_states[i].first;
        // std::cout << "S curr_s: " << curr_s << std::endl;
        for (auto &p : det_cache[curr_s])
        {
          // only care about the transitions in the same SCC
          if (si_.scc_of(curr_s) != si_.scc_of(p.second))
          {
            continue;
          }
          // accepting state or accepting transition, abort
          if (p.first || is_accepting_[p.second])
            // Exit early; accepting transition is forbidden for safe state.
            return;
          // std::cout << "S succ: " << p.second << std::endl;
          // states are already safe will stay safe forever
          next_detstates.erase(p.second);
          succ_nodes[p.second] = NCSB_S;
          // No need to look for other compatible transitions
          // for this state; it's in the deterministic in the same SCC
          break;
        }
      }

      std::set<unsigned> scc_indices;
      //2. Handle C states.
      for (unsigned i = 0; i < curr_det_states.size(); ++i)
      {
        // including B-states
        if (!(curr_det_states[i].second & NCSB_C))
          continue;

        unsigned curr_s = curr_det_states[i].first;
        // std::cout << "C curr_s: " << curr_s << std::endl;
        for (auto &p : det_cache[curr_s])
        {
          // only care about the transitions in the same SCC
          if (si_.scc_of(curr_s) != si_.scc_of(p.second))
          {
            continue;
          }

          next_detstates.erase(p.second);
          // Ignore states that goes to S
          if (succ_nodes.find(p.second) == succ_nodes.end())
          {
            // std::cout << "C succ: " << p.second << std::endl;
            succ_nodes[p.second] = NCSB_C;
            scc_indices.insert(si_.scc_of(p.second));
          }
          break;
        }
      }

      //3. Handle incoming states.
      for (unsigned p : next_detstates)
      {
        if (succ_nodes.find(p) == succ_nodes.end())
        {
          // all incoming states will be set in C
          succ_nodes[p] = NCSB_C;
          // std::cout << "C incoming curr_s: " << p << std::endl;
          scc_indices.insert(si_.scc_of(p));
        }
      }


      //4. Handle B-states
      bool is_pre_b_empty = true;
      bool is_b_empty = true;
      for (unsigned i = 0; i < curr_det_states.size(); ++i)
      {
        // including B-states
        if (curr_det_states[i].second != NCSB_B)
          continue;

        is_b_empty = false; 
        is_pre_b_empty = false;
        unsigned curr_s = curr_det_states[i].first;
        for (auto &p : det_cache[curr_s])
        {
          if (si_.scc_of(curr_s) != si_.scc_of(p.second))
          {
            continue;
          }
          if (succ_nodes[p.second] == NCSB_C)
          {
            // std::cout << "B succ: " << p.second << std::endl;
            succ_nodes[p.second] = NCSB_B;
            is_b_empty = false;
          }
          break;
        }
      }

      int curr_scc_index;
      if (is_pre_b_empty)
      {
        // the DACs has just been reached
        curr_scc_index = ((int)acc_detsccs_.size()) - 1;
      }
      else
      {
        curr_scc_index = scc_index;
      }
      // std::cout << "curr_scc_index = " << curr_scc_index << std::endl;

      int next_scc_index;
      // round rabin strategy for next DAC, we should select existing one if it is not 0
      if (curr_scc_index == 0)
      {
        next_scc_index = acc_detsccs_.size() - 1;
      }
      else
      {
        next_scc_index = curr_scc_index - 1;
      }
      // std::cout << "next_scc_index: " << next_scc_index << std::endl;
      unsigned curr_scc = acc_detsccs_[next_scc_index];
      // std::cout << "Current scc: " << curr_scc << std::endl;

      if (next_scc_index != 0 && scc_indices.find(curr_scc) == scc_indices.end())
      {
        // need to find an index inside scc_indices
        int max_lower = -1;
        int max_upper = -1;
        // already sorted
        for (unsigned index : scc_indices)
        {
          if (index < curr_scc)
          {
            max_lower = std::max(max_lower, (int)index);
          }
          if (index > curr_scc)
          {
            max_upper = std::max(max_upper, (int)index);
          }
        }
        // std::cout << "C scc: " << get_set_string(scc_indices) << std::endl;
        // std::cout << "max_lower = " << max_lower << " max_upper = " << max_upper << std::endl;
        if (max_lower != -1)
        {
          next_scc_index = get_detscc_index(max_lower);
        }
        else if (max_upper != -1)
        {
          next_scc_index = get_detscc_index(max_upper);
        }
        else
        {
          // C' maybe empty, so set it to 0
          next_scc_index = 0;
        }
      }
      std::map<unsigned, int> tmp(succ_nodes);
      //5. Now add the first successor
      // next_scc_index = true_index;
      if (is_b_empty)
      {
        acc_succs.emplace_back(true);
        // round rabin for checking next DAC
        for (auto &p : succ_nodes)
        {
          if (p.second == NCSB_C && si_.scc_of(p.first) == active_index)
          {
            tmp[p.first] = NCSB_B;
          }
        }
        next_scc_indices.emplace_back(next_scc_index);
      }
      else
      {
        acc_succs.emplace_back(false);
        next_scc_indices.emplace_back(scc_index);
      }
      // std::cout << "Add first map: " << std::endl;
      // for (auto & p : tmp)
      // {
      //   std::cout << "s = " << p.first << " r = " << p.second << std::endl;
      // }
      succ_maps.emplace_back(tmp);

      //6. MaxRank - another successor
      if (!is_b_empty)
      {
        for (auto &p : succ_nodes)
        {
          // B' has accepting states
          if (is_accepting_[p.first] && p.second == NCSB_B)
          {
            is_b_empty = true;
            break;
          }
        }
      }
     
      if (!is_b_empty)
      {
        for (auto mp : det_cache)
        {
          std::pair<unsigned, int> pr{mp.first, NCSB_B};
          if (std::find(curr_det_states.begin(), curr_det_states.end(), pr) != curr_det_states.end()){
            for (auto dst : mp.second)
            {
              // we care only about accepting transitions in the same scc
              if (si_.scc_of(mp.first) == si_.scc_of(dst.second))
              {
                if (dst.first and succ_nodes[dst.second] == NCSB_B)
                {
                  is_b_empty = true;
                  break;
                }
              }
            }
          }
        }
      }
      
      if (!is_b_empty)
      {
        std::map<unsigned, int> tmp2(succ_nodes);
        for (auto &p : succ_nodes)
        {
          // Move all B'-states to S
          if (p.second == NCSB_B)
          {
            tmp2[p.first] = NCSB_S;
          }
        }
        for (auto &p : succ_nodes)
        {
          // Move all B'-states to S
          if (p.second == NCSB_C && si_.scc_of(p.first) == active_index)
          {
            tmp2[p.first] = NCSB_B;
          }
        }
        // std::cout << "Add second map: " << std::endl;
        // for (auto & p : tmp2)
        // {
        //   std::cout << "s = " << p.first << " r = " << p.second << std::endl;
        // }
        succ_maps.emplace_back(tmp2);
        acc_succs.emplace_back(true);
        next_scc_indices.emplace_back(next_scc_index);
      }
    }

    // adapted CSB complementation, every part may have at most two successors
    void det_succ(const complement_mstate &ms, std::vector<std::vector<state_rank>> &succs, std::vector<bool> &acc_succs, std::vector<int> &next_scc_index, std::set<unsigned> &next_detstates, std::unordered_map<unsigned, std::vector<std::pair<bool, unsigned>>> &det_cache)
    {
      std::vector<std::map<unsigned, int>> succ_maps;
      //csb_successors(ms.detscc_ranks_, ms.detscc_index_, next_scc_index, succ_maps, acc_succs, next_detstates, det_cache);
      // std::cout << "map size: " << succ_maps.size() << std::endl;
      for (unsigned i = 0; i < succ_maps.size(); i++)
      {
        std::vector<state_rank> succ;
        for (auto &p : succ_maps[i])
        {
          succ.emplace_back(p.first, p.second);
        }
        // std::cout << "next " << next_scc_index[i] << " rank: " << get_rank_string(succ) << std::endl;
        succs.emplace_back(succ);
      }
    }

    // compute the successor P={nondeterministic states and nonaccepting SCCs} O = {breakpoint for weak SCCs}
    // and labelling states for each SCC
    void
    compute_successors(const complement_mstate &ms, unsigned origin, bdd letter)
    {
      // std::cout << "current state: " << get_name(ms) << " src: " << origin << " letter: " << letter << std::endl;
      complement_mstate succ(si_, acc_detsccs_.size());
      // used for unambiguous automaton
      std::vector<bool> incoming(nb_states_, false);
      std::vector<bool> ignores(nb_states_, false);

      // this function is used for unambiguous NBAs
      auto can_ignore = [&incoming, &ignores](bool use_ambiguous, unsigned dst) -> bool
      {
        if (use_ambiguous)
        {
          if (incoming[dst])
          {
            // this is the second incoming transitions
            ignores[dst] = true;
          }
          else
          {
            incoming[dst] = true;
          }
          return ignores[dst];
        }
        else
        {
          return false;
        }
      };

      std::set<unsigned> next_level_states;
      std::set<unsigned> acc_weak_coming_states;
      // states at current level
      std::set<unsigned> current_states = ms.get_reach_set();
      // states at next level
      std::set<unsigned> next_nondetstates;
      std::vector<std::set<unsigned>> next_detstates;
      for (unsigned i = 0; i < acc_detsccs_.size(); i++)
      {
        next_detstates.emplace_back(std::set<unsigned>());
      }

      // deterministic transitions
      std::unordered_map<unsigned, std::vector<std::pair<bool, unsigned>>> det_cache;
      // nondeterministic transitions
      std::unordered_map<unsigned, std::vector<std::pair<bool, unsigned>>> nondet_cache;

      //1. first handle inherently weak states
      for (unsigned s : current_states)
      {
        // nondeterministic states or states in nonaccepting SCCs
        bool in_break_set = (ms.break_set_.find(s) != ms.break_set_.end());
        bool in_acc_det = is_accepting_detscc(scc_types_, si_.scc_of(s));
        if (in_acc_det)
        {
          det_cache.emplace(s, std::vector<std::pair<bool, unsigned>>());
        }
        bool in_acc_nondet = is_accepting_nondetscc(scc_types_, si_.scc_of(s));
        if (in_acc_nondet)
        {
          nondet_cache.emplace(s, std::vector<std::pair<bool, unsigned>>());
        }
        for (const auto &t : aut_->out(s))
        {
          if (!bdd_implies(letter, t.cond))
            continue;
          // it is legal to ignore the states have two incoming transitions
          // in unambiguous Buchi automaton
          if (can_ignore(use_unambiguous_, t.dst))
            continue;
          next_level_states.insert(t.dst);
          unsigned scc_id = si_.scc_of(t.dst);
          // we move the states in accepting det SCC to ordered states
          if (is_accepting_detscc(scc_types_, scc_id))
          {
            int scc_index = get_detscc_index(scc_id);
            next_detstates[scc_index].insert(t.dst);
            if (in_acc_det)
            {
              det_cache[s].emplace_back(t.acc.has(0), t.dst);
            }
          }
          else if (is_weakscc(scc_types_, scc_id))
          {
            // weak states or nondeterministic or nonaccepting det scc
            succ.weak_set_.insert(t.dst);
            // be accepting and weak
            bool in_acc_set = (scc_types_[scc_id] & SCC_ACC) > 0;
            // in breakpoint and it is accepting
            if (in_break_set && in_acc_set)
            {
              succ.break_set_.insert(t.dst);
            }
            // in accepting weak SCCs
            if (in_acc_set)
            {
              acc_weak_coming_states.insert(t.dst);
            }
          }
          else if (is_accepting_nondetscc(scc_types_, scc_id))
          {
            next_nondetstates.insert(t.dst);
            if (in_acc_nondet)
            {
              nondet_cache[s].emplace_back(t.acc.has(0), t.dst);
            }
          }
          else
          {
            assert(false);
          }
        }
      }
      // std::cout << "det: " << get_set_string(next_detstates) << std::endl;
      // std::cout << "nondet: " << get_set_string(next_nondetstates) << std::endl;
      // std::cout << "After weak: " << get_name(succ) << std::endl;
      //2. Compute the successors in deterministic SCCs
      std::vector<int> det_colors;
      for (unsigned i = 0; i < acc_detsccs_.size(); i++)
      {
        compute_det_succ compute_det(si_, acc_detsccs_[i], ms.detscc_ranks_[i], next_detstates[i], det_cache);
        compute_det.compute();
        succ.detscc_ranks_[i] = compute_det.next_ranks_;
        det_colors.emplace_back(compute_det.get_color());
      }
      // std::cout << "After deterministic part = " << get_name(succ) << std::endl;

      //3. Compute the successors for nondeterministic SCCs
      // compute_nondet_succ compute_nondet(si_, acc_nondetsccs_, is_accepting_, ms.nondetscc_ranks_, ms.nondetscc_marks_, next_nondetstates, nondet_cache);
      // compute_nondet.compute();
      // at most two successors
      // std::vector<std::vector<std::set<unsigned>>> &nondet_successors = compute_nondet.next_slices_;
      // std::vector<std::vector<slice_mark>> &nondet_marks = compute_nondet.next_marks_;
      // std::vector<int> nondet_colors = compute_nondet.acc_colors_;
      // compute_nondeterministic_successors(ms, succ, next_nondetstates, nondet_cache);
      // std::cout << "After nondeterministic part = " << get_name(succ) << std::endl;

      // remove redudant states
      if (use_simulation_)
      {
        //make_simulation_state(succ, next_level_states, det_successors, nondet_successors);
        //merge_redundant_states(succ, det_successors, nondet_successors);
      }

      bool break_empty = succ.break_set_.empty();
      // now determine the break set
      if (break_empty)
      {
        // if the breakpoint is empty, then fill it with newly-incoming accepting weak SCC states
        std::set<unsigned> result;
        std::set<unsigned> reach_sucss = succ.weak_set_; // copy
        std::set_intersection(reach_sucss.begin(), reach_sucss.end(), acc_weak_coming_states.begin(), acc_weak_coming_states.end(), std::inserter(result, result.begin()));
        succ.break_set_ = result;
      }
    }

    int get_nondetscc_index(unsigned scc)
    {
      for (int idx = 0; idx < acc_nondetsccs_.size(); idx++)
      {
        if (acc_nondetsccs_[idx] == scc)
        {
          return idx;
        }
      }
      return -1;
    }
    int get_detscc_index(unsigned scc)
    {
      for (int idx = 0; idx < acc_detsccs_.size(); idx++)
      {
        if (acc_detsccs_[idx] == scc)
        {
          return idx;
        }
      }
      return -1;
    }

    std::vector<complement_mstate>
    get_succ_active_CSB(std::set<unsigned> reachable, bdd letter, std::vector<state_rank> det_ranks, int active_index, complement_mstate &new_succ, complement_mstate &new_succ2, std::vector<std::pair<std::vector<unsigned>, std::vector<unsigned>>>& acc_det_succ, unsigned true_index)
    {
      std::vector<complement_mstate> succ;

      std::vector<std::set<unsigned>> next_detstates;
      for (unsigned i = 0; i < acc_detsccs_.size(); i++)
      {
        next_detstates.emplace_back(std::set<unsigned>());
      }
      // deterministic transitions
      std::unordered_map<unsigned, std::vector<std::pair<bool, unsigned>>> det_cache;

      for (unsigned s : reachable)
      {
        det_cache.emplace(s, std::vector<std::pair<bool, unsigned>>());
        for (const auto &t : aut_->out(s))
        {
          if (!bdd_implies(letter, t.cond))
            continue;
          unsigned scc_id = si_.scc_of(t.dst);
          // we move the states in accepting det SCC to ordered states
          if (is_accepting_detscc(scc_types_, scc_id))
          {
            int scc_index = get_detscc_index(scc_id);
            next_detstates[scc_index].insert(t.dst);
            det_cache[s].emplace_back(t.acc.has(0), t.dst);
          }
        }
      }
      
      std::vector<int> next_scc_indices;
      std::vector<std::map<unsigned, int>> succ_maps;
      std::vector<bool> acc_succs;
      
      csb_successors(det_ranks, active_index, next_scc_indices, succ_maps, acc_succs, next_detstates[get_detscc_index(active_index)], det_cache, active_index);

      bool first = true;
      bool repeat = false;
      for (auto s : succ_maps)
      {
        std::set<unsigned> C;
        std::set<unsigned> S;
        std::set<unsigned> B;
        for (auto pr : s)
        {
          if (first)
            new_succ.detscc_ranks_[get_detscc_index(active_index)].push_back({pr.first, pr.second});
          else
            new_succ2.detscc_ranks_[get_detscc_index(active_index)].push_back({pr.first, pr.second});

          if (pr.second == NCSB_B)
          {
            B.insert(pr.first);
            C.insert(pr.first);
          }
          if (pr.second == NCSB_C)
          {
            C.insert(pr.first);
          }
          else if (pr.second == NCSB_S)
          {
            S.insert(pr.first);
          }
        }

        if (C.empty())
          repeat = true;
        
        acc_det_succ[true_index] = {std::vector<unsigned>(C.begin(), C.end()), std::vector<unsigned>(S.begin(), S.end())};

        if (first)
        {
          auto B_orig = B;
          auto C_orig = C;
          auto S_orig = S;

          B.insert(new_succ.det_break_set_.begin(), new_succ.det_break_set_.end());
          new_succ.det_break_set_ = std::vector<unsigned>(B.begin(), B.end());

          if (new_succ.acc_detsccs_.empty())
            new_succ.set_acc_detsccs(acc_det_succ);

          C.insert(new_succ.acc_detsccs_[true_index].first.begin(), new_succ.acc_detsccs_[true_index].first.end());
          S.insert(new_succ.acc_detsccs_[true_index].second.begin(), new_succ.acc_detsccs_[true_index].second.end());
          new_succ.acc_detsccs_[true_index] = {std::vector<unsigned>(C.begin(), C.end()), std::vector<unsigned>(S.begin(), S.end())};
          //new_succ.acc_detsccs_[true_index] = {std::vector<unsigned>(C.begin(), C.end()), std::vector<unsigned>(S.begin(), S.end())};
          
          succ.push_back(new_succ);
          first = false;

          if (repeat and succ_maps.size() == 1)
          {
            B_orig.insert(new_succ2.det_break_set_.begin(), new_succ2.det_break_set_.end());
            new_succ2.det_break_set_ = std::vector<unsigned>(B_orig.begin(), B_orig.end());
            
            if (new_succ2.acc_detsccs_.empty())
              new_succ2.set_acc_detsccs(acc_det_succ);
            
            C_orig.insert(new_succ2.acc_detsccs_[true_index].first.begin(), new_succ2.acc_detsccs_[true_index].first.end());
            S_orig.insert(new_succ2.acc_detsccs_[true_index].second.begin(), new_succ2.acc_detsccs_[true_index].second.end());
            new_succ2.acc_detsccs_[true_index] = {std::vector<unsigned>(C_orig.begin(), C_orig.end()), std::vector<unsigned>(S_orig.begin(), S_orig.end())};
            
            succ.push_back(new_succ2);
          }
        }
        else 
        {
          B.insert(new_succ2.det_break_set_.begin(), new_succ2.det_break_set_.end());
          new_succ2.det_break_set_ = std::vector<unsigned>(B.begin(), B.end());

          if (new_succ2.acc_detsccs_.empty())
            new_succ2.set_acc_detsccs(acc_det_succ);

          C.insert(new_succ2.acc_detsccs_[true_index].first.begin(), new_succ2.acc_detsccs_[true_index].first.end());
          S.insert(new_succ2.acc_detsccs_[true_index].second.begin(), new_succ2.acc_detsccs_[true_index].second.end());
          new_succ2.acc_detsccs_[true_index] = {std::vector<unsigned>(C.begin(), C.end()), std::vector<unsigned>(S.begin(), S.end())};

          succ.push_back(new_succ2);
        }
      }

      if (succ.size() == 1)
        new_succ2.set_active_index(-2);

      return succ;
    }

    std::vector<complement_mstate>
    get_succ_track_CSB(std::set<unsigned> reachable, bdd letter, std::vector<state_rank> det_ranks, int active_index, complement_mstate &new_succ, complement_mstate &new_succ2, std::vector<std::pair<std::vector<unsigned>, std::vector<unsigned>>>& acc_det_succ, unsigned true_index)
    {
      std::vector<complement_mstate> succ;

      std::vector<std::set<unsigned>> next_detstates;
      for (unsigned i = 0; i < acc_detsccs_.size(); i++)
      {
        next_detstates.emplace_back(std::set<unsigned>());
      }
      // deterministic transitions
      std::unordered_map<unsigned, std::vector<std::pair<bool, unsigned>>> det_cache;

      for (unsigned s : reachable)
      {
        det_cache.emplace(s, std::vector<std::pair<bool, unsigned>>());
        for (const auto &t : aut_->out(s))
        {
          if (!bdd_implies(letter, t.cond))
            continue;
          unsigned scc_id = si_.scc_of(t.dst);
          // we move the states in accepting det SCC to ordered states
          if (is_accepting_detscc(scc_types_, scc_id))
          {
            int scc_index = get_detscc_index(scc_id);
            next_detstates[scc_index].insert(t.dst);
            det_cache[s].emplace_back(t.acc.has(0), t.dst);
          }
        }
      }
      
      std::vector<int> next_scc_indices;
      std::vector<std::map<unsigned, int>> succ_maps;
      std::vector<bool> acc_succs;
      
      csb_successors(det_ranks, active_index, next_scc_indices, succ_maps, acc_succs, next_detstates[get_detscc_index(active_index)], det_cache, active_index);
      // complement_mstate new_succ2(new_succ);
      bool first = true;
      for (auto s : succ_maps)
      {
        std::set<unsigned> C;
        std::set<unsigned> S;
        std::set<unsigned> B;
        for (auto pr : s)
        {
          if (first)
            new_succ.detscc_ranks_[get_detscc_index(active_index)].push_back({pr.first, pr.second});
          else
            new_succ2.detscc_ranks_[get_detscc_index(active_index)].push_back({pr.first, pr.second});
          if (pr.second == NCSB_B)
          {
            B.insert(pr.first);
            C.insert(pr.first);
          }
          if (pr.second == NCSB_C)
          {
            C.insert(pr.first);
          }
          else if (pr.second == NCSB_S)
          {
            S.insert(pr.first);
          }
        }
        
        if (first)
        {
          if (new_succ.acc_detsccs_.empty())
            new_succ.set_acc_detsccs(acc_det_succ); 

          if (new_succ2.acc_detsccs_.empty())
            new_succ2.set_acc_detsccs(acc_det_succ);
          
          C.insert(new_succ.acc_detsccs_[true_index].first.begin(), new_succ.acc_detsccs_[true_index].first.end());
          S.insert(new_succ.acc_detsccs_[true_index].second.begin(), new_succ.acc_detsccs_[true_index].second.end());
          new_succ.acc_detsccs_[true_index] = {std::vector<unsigned>(C.begin(), C.end()), std::vector<unsigned>(S.begin(), S.end())};
          new_succ2.acc_detsccs_[true_index] = {std::vector<unsigned>(C.begin(), C.end()), std::vector<unsigned>(S.begin(), S.end())};
          
          succ.push_back(new_succ);
          first = false;
        }
        // else 
        // {
        //   if (new_succ2.acc_detsccs_.empty()) 
        //     new_succ2.set_acc_detsccs(acc_det_succ);
        //   succ.push_back(new_succ2);
        // }
      }

      return succ;
    }

  public:
    tnba_complement(const spot::const_twa_graph_ptr &aut, spot::scc_info &si, spot::option_map &om, std::vector<bdd> &implications, compl_decomp_options& decomp_options)
        : aut_(aut),
          om_(om),
          decomp_options_(decomp_options),
          use_simulation_(om.get(USE_SIMULATION) > 0),
          use_scc_(om.get(USE_SCC_INFO) > 0),
          use_stutter_(om.get(USE_STUTTER) > 0),
          use_unambiguous_(om.get(USE_UNAMBIGUITY) > 0),
          si_(si),
          nb_states_(aut->num_states()),
          support_(nb_states_),
          compat_(nb_states_),
          is_accepting_(aut->num_states(), false),
          simulator_(aut, si, implications, om.get(USE_SIMULATION) > 0),
          delayed_simulator_(aut, om),
          show_names_(om.get(VERBOSE_LEVEL) >= 1)
    {
      
      if (om.get(VERBOSE_LEVEL) >= 2)
      {
        simulator_.output_simulation();
      }
      res_ = spot::make_twa_graph(aut->get_dict());
      res_->copy_ap_of(aut);
      res_->prop_copy(aut,
                        {
                            false,        // state based
                            false,        // inherently_weak
                            false, false, // deterministic
                            true,         // complete
                            false         // stutter inv
                        });
      // Generate bdd supports and compatible options for each state.
      // Also check if all its transitions are accepting.
      for (unsigned i = 0; i < nb_states_; ++i)
      {
        bdd res_support = bddtrue;
        bdd res_compat = bddfalse;
        bool accepting = true;
        bool has_transitions = false;
        for (const auto &out : aut->out(i))
        {
          has_transitions = true;
          res_support &= bdd_support(out.cond);
          res_compat |= out.cond;
          if (!out.acc)
            accepting = false;
        }
        support_[i] = res_support;
        compat_[i] = res_compat;
        is_accepting_[i] = accepting && has_transitions;
      }
      // obtain the types of each SCC
      scc_types_ = get_scc_types(si_);
      // find out the DACs and NACs
      unsigned nonacc_weak = 0;
      for (unsigned i = 0; i < scc_types_.size(); i++)
      {
        if (is_accepting_detscc(scc_types_, i))
        {
          acc_detsccs_.push_back(i);
        }
        else if (is_accepting_weakscc(scc_types_, i))
        {
          weaksccs_.push_back(i);
          // if (not is_accepting_weakscc(scc_types_, i))
          //   nonacc_weak++;
        }
        else if (is_accepting_nondetscc(scc_types_, i))
        {
          // accepting nondeterministic scc
          acc_nondetsccs_.emplace_back(i);
        }
      }

      std::cerr << "SCCs: " << acc_detsccs_.size() + weaksccs_.size() << ", DET: " << acc_detsccs_.size() << ", IWA: " << weaksccs_.size() << std::endl << std::endl; 
    }

    unsigned
    get_num_states()
    {
      return this->nb_states_;
    }

    spot::scc_info&
    get_scc_info()
    {
      return this->si_;
    }

    spot::twa_graph_ptr
    run()
    {
      if (show_names_)
      {
        names_ = new std::vector<std::string>();
        res_->set_named_prop("state-names", names_);
      }

      if (this->weaksccs_.size() == 0)
        decomp_options_.merge_iwa = false;
      if (this->acc_detsccs_.size() == 0)
      {
        decomp_options_.merge_det = false;
        decomp_options_.tgba = false; // no TGBA for IW SCCs only
      }

      // complementation algorithm
      // auto acc = res_->set_buchi();
      if (decomp_options_.tgba)
        res_->set_generalized_buchi(2);
      else 
        res_->set_generalized_buchi(1);

      spot::print_hoa(std::cerr, aut_);
      std::cerr << std::endl << std::endl;

      // initial macrostate
      auto scc_info = get_scc_info();
      complement_mstate init_state(scc_info, acc_detsccs_.size());
      unsigned orig_init = aut_->get_init_state_number();
      int active_index = scc_info.scc_of(orig_init);

      if (decomp_options_.tgba)
      {
         init_state.set_active_index(this->acc_detsccs_[0]);
         active_index = this->acc_detsccs_[0];
      }
      else if (decomp_options_.merge_iwa and is_weakscc(scc_types_, active_index))
      {
        init_state.set_active_index(this->weaksccs_[0]);
      }
      else if (decomp_options_.merge_det and is_accepting_detscc(scc_types_, active_index))
      {
        init_state.set_active_index(this->acc_detsccs_[0]);
      }
      else
      {
        if (is_weakscc(scc_types_, active_index) and (not is_accepting_weakscc(scc_types_, active_index)))
        {
          // initial state in nonaccepting scc
          unsigned tmp_index = active_index;
          do 
          {
            active_index = (active_index + 1) % si_.scc_count();

            if (active_index == tmp_index)
              break;
          } while (is_weakscc(scc_types_, active_index) and (not is_accepting_weakscc(scc_types_, active_index)));
          init_state.set_active_index(active_index);
        }
        else
        {
          init_state.set_active_index(active_index);
        }
      }

      std::vector<complement_mstate> all_states;

      std::vector<std::vector<unsigned>> iw_sccs;
      std::vector<std::pair<std::vector<unsigned>, std::vector<unsigned>>> acc_detsccs;
      
      // weak SCCs
      for (unsigned index : weaksccs_)
      {
        if (index != active_index or active_index != si_.scc_of(orig_init))
          iw_sccs.push_back(std::vector<unsigned>());
        else
          iw_sccs.push_back(std::vector<unsigned>(1, orig_init));
      }
      if (decomp_options_.merge_iwa)
      {
        iw_sccs.clear();
        if (is_weakscc(scc_types_, active_index))
          iw_sccs.push_back(std::vector<unsigned>(1, orig_init));
        else  
          iw_sccs.push_back(std::vector<unsigned>());
      }

      // det SCCs
      if (not decomp_options_.merge_det)
      {
        for (unsigned index : acc_detsccs_)
        {
          if (index != active_index or scc_info.scc_of(orig_init) != active_index)
            acc_detsccs.push_back({std::vector<unsigned>(), std::vector<unsigned>()});
          else
          {
            acc_detsccs.push_back({std::vector<unsigned>(1, orig_init), std::vector<unsigned>()});
            init_state.detscc_ranks_[get_detscc_index(active_index)].emplace_back(orig_init, NCSB_C); 
            init_state.detscc_ranks_[get_detscc_index(active_index)].emplace_back(orig_init, NCSB_B); 
          }
        }
      }
      else 
      {
        if (is_accepting_detscc(scc_types_, active_index))
        {
          acc_detsccs.push_back({std::vector<unsigned>(1, orig_init), std::vector<unsigned>()});
          for (unsigned index : acc_detsccs_)
          {
            init_state.detscc_ranks_[get_detscc_index(active_index)].emplace_back(orig_init, NCSB_C); 
            init_state.detscc_ranks_[get_detscc_index(active_index)].emplace_back(orig_init, NCSB_B); 
          }
        }
        else
        {
          acc_detsccs.push_back({std::vector<unsigned>(), std::vector<unsigned>()});
        }
      }
      
      init_state.set_iw_sccs(iw_sccs);
      init_state.set_acc_detsccs(acc_detsccs);
      auto acc_detsccs_orig = acc_detsccs;
      init_state.curr_reachable_.push_back(orig_init);

      // get break set for active scc
      if (is_weakscc(scc_types_, active_index) and active_index == scc_info.scc_of(orig_init))
      { 
        init_state.set_iw_break_set(std::vector<unsigned>(1, orig_init)); 
        init_state.det_break_set_ = std::vector<unsigned>();
      }
      else if (is_accepting_detscc(scc_types_, active_index) and active_index == scc_info.scc_of(orig_init))
      {
        init_state.det_break_set_ = std::vector<unsigned>(1, orig_init); 
        init_state.iw_break_set_ = std::vector<unsigned>();
      }
      else
      {
        init_state.set_iw_break_set(std::vector<unsigned>());
        init_state.det_break_set_ = std::vector<unsigned>();
      }

      if (scc_info.scc_of(orig_init) != active_index)
        init_state.set_iw_break_set(std::vector<unsigned>());

      std::cerr << "Initial: " << get_name(init_state) << std::endl;
      auto init = new_state(init_state);
      res_->set_init_state(init); 

      all_states.push_back(init_state);

      mh_complement mh(aut_, scc_info, scc_types_);
      bool sink_state = false;

      while (!todo_.empty())
      { 
        auto top = todo_.front();
        todo_.pop_front();
        complement_mstate ms = top.first;
        std::cerr << std::endl << "State: " << get_name(ms) << std::endl;
        active_index = ms.active_index_;

        if (active_index >= 0 and is_weakscc(scc_types_, active_index) and (not is_accepting_weakscc(scc_types_, active_index)))
        {
          ms.active_index_ = (ms.active_index_ + 1) % si_.scc_count();
          todo_.emplace_back(ms, top.second);
          continue;
        }

        // reachable states
        std::set<unsigned> reachable = std::set<unsigned>(ms.curr_reachable_.begin(), ms.curr_reachable_.end());

        // Compute support of all available states.
        bdd msupport = bddtrue;
        bdd n_s_compat = bddfalse;
        const std::set<unsigned> &reach_set = ms.get_reach_set();
        // compute the occurred variables in the outgoing transitions of ms, stored in msupport
        for (unsigned s : reach_set)
        {
          msupport &= support_[s];
          n_s_compat |= compat_[s];
        }

        bdd all = n_s_compat;
        if (all != bddtrue)
        {
          // direct the rest to sink state
          complement_mstate succ(si_, acc_detsccs_.size());
          succ.active_index_ = -1;
          auto sink = new_state(succ);
          // empty state use 0 as well as the weak ones
          res_->new_edge(top.second, sink, !all);
          if (not sink_state)
          {
            res_->new_edge(sink, sink, !all, {0});
            res_->new_edge(sink, sink, all, {0});
            if (decomp_options_.tgba)
            {
              res_->new_edge(sink, sink, !all, {1});
              res_->new_edge(sink, sink, all, {1});
            }
            sink_state = true;
          }
        }
        
        while (all != bddfalse)
        {
          bdd letter = bdd_satoneset(all, msupport, bddfalse);
          all -= letter;
          std::cerr << "Current symbol: " << letter << std::endl;

          std::set<unsigned> all_succ = mh.get_all_successors(reachable, letter);

          bool active_type = true;
          bool active_type2 = true;
          bool no_succ = false;
          
          std::vector<complement_mstate> new_succ;
          complement_mstate new_succ1(scc_info, acc_detsccs_.size());
          complement_mstate new_succ2(scc_info, acc_detsccs_.size()); 
          new_succ.push_back(new_succ1);
          new_succ.push_back(new_succ2);

          std::vector<std::vector<unsigned>> iw_succ(ms.iw_sccs_.size());
          if (decomp_options_.merge_iwa)
          {
            iw_succ.clear();
            iw_succ.push_back(std::vector<unsigned>());
          }
          std::vector<std::pair<std::vector<unsigned>, std::vector<unsigned>>> acc_det_succ;
          if (not decomp_options_.merge_det)
          {
            for (unsigned i=0; i<ms.acc_detsccs_.size(); i++)
            {
              acc_det_succ.push_back({std::vector<unsigned>(), std::vector<unsigned>()});
            }
          }
          else
          {
            acc_det_succ.push_back({std::vector<unsigned>(), std::vector<unsigned>()});
          }

          // scc indices
          std::vector<unsigned> indices;
          if (not decomp_options_.merge_iwa)
            indices.insert(indices.end(), this->weaksccs_.begin(), this->weaksccs_.end());
          else if (this->weaksccs_.size() > 0)
            indices.push_back(this->weaksccs_[0]);
          if (not decomp_options_.merge_det)
            indices.insert(indices.end(), this->acc_detsccs_.begin(), this->acc_detsccs_.end());
          else if (this->acc_detsccs_.size() > 0)
            indices.push_back(this->acc_detsccs_[0]);
          
          // index of value active_index
          auto it = std::find(indices.begin(), indices.end(), active_index);
          unsigned true_index = std::distance(indices.begin(), it);
          unsigned orig_index = true_index;

          std::vector<complement_mstate> succ_det;

          if (ms.iw_break_set_.size() == 0 and ms.det_break_set_.size() == 0)
            active_type = false;

          
          bool iwa_done = false;
          bool det_done = false;

          for (unsigned i=0; i<indices.size(); i++)
          {
            true_index = (orig_index+i)%indices.size();
            
            std::vector<unsigned> index;
            if (decomp_options_.merge_iwa and is_weakscc(scc_types_, indices[true_index]))
              index.insert(index.end(), weaksccs_.begin(), weaksccs_.end());
            else if (decomp_options_.merge_det and is_accepting_detscc(scc_types_, indices[true_index]))
              index.insert(index.end(), acc_detsccs_.begin(), acc_detsccs_.end());
            else
              index.push_back(indices[true_index]);

            // merge iwa
            if (decomp_options_.merge_iwa and is_weakscc(scc_types_, index[0]) and iwa_done)
              continue; 
            // merge det
            if (decomp_options_.merge_det and is_accepting_detscc(scc_types_, index[0]) and det_done)
              continue; 

            // reachable states in this scc
            std::set<unsigned> reach_track;
            std::set<unsigned> scc_states;
            if (decomp_options_.merge_iwa and is_weakscc(scc_types_, index[0]))
            {
              for (auto i : weaksccs_)
              {
                scc_states.insert(scc_info.states_of(i).begin(), scc_info.states_of(i).end());
              }
            }
            else if (decomp_options_.merge_det and is_accepting_detscc(scc_types_, index[0]))
            {
              for (auto i : acc_detsccs_)
              {
                scc_states.insert(scc_info.states_of(i).begin(), scc_info.states_of(i).end());
              }
            }
            else 
            {
              scc_states.insert(scc_info.states_of(index[0]).begin(), scc_info.states_of(index[0]).end());
            }
            std::set_intersection(scc_states.begin(), scc_states.end(), reachable.begin(), reachable.end(), std::inserter(reach_track, reach_track.begin()));

            // successors in this scc
            std::set<unsigned> succ_in_scc;
            std::set_intersection(scc_states.begin(), scc_states.end(), all_succ.begin(), all_succ.end(), std::inserter(succ_in_scc, succ_in_scc.begin()));

            // non-active component
            if (std::find(index.begin(), index.end(), active_index) == index.end() and (not decomp_options_.tgba or not is_weakscc(scc_types_, index[0]))) // TODO check!
            {
              std::cerr << "Non-active" << std::endl;
              // non-active scc
              if (is_weakscc(scc_types_, index[0]))
              {
                iwa_done = true;
                // getSuccTrack
                if (active_type or ((index[0] != (indices[(orig_index + 1)%indices.size()]))) or (not is_accepting_weakscc(scc_types_, index[0])))
                {
                  auto succ_track = mh.get_succ_track(reachable, reach_track, letter, index);
                  for (auto succ : succ_track)
                  {
                    if (decomp_options_.merge_iwa)
                      iw_succ[0] = std::vector<unsigned>(succ.first.begin(), succ.first.end());
                    else
                      iw_succ[true_index] = std::vector<unsigned>(succ.first.begin(), succ.first.end());
                  }
                }

                // getSuccTrackToActive
                else
                {
                  auto succ_track_to_active = mh.get_succ_track_to_active(reachable, reach_track, letter, index);
                  for (auto succ : succ_track_to_active) 
                  {
                    if (decomp_options_.merge_iwa)
                      iw_succ[0] = std::vector<unsigned>(succ.first.first.begin(), succ.first.first.end());
                    else
                      iw_succ[true_index] = std::vector<unsigned>(succ.first.first.begin(), succ.first.first.end());
                    for (unsigned i=0; i<new_succ.size(); i++)
                    {
                      new_succ[i].set_iw_break_set(std::vector<unsigned>(succ.first.second.begin(), succ.first.second.end()));
                    }
                  }
                }
              }

              else if (is_accepting_detscc(scc_types_, index[0]))
              {
                if (active_type or (index[0] != (indices[(orig_index + 1)%indices.size()])))
                {
                  // getSuccTrack
                  std::vector<state_rank> ranks;
                  if (not decomp_options_.merge_det)
                  {
                    for (auto s : ms.acc_detsccs_[true_index - iw_succ.size()].first)
                    {
                      ranks.push_back({s, NCSB_C}); 
                      ranks.push_back({s, NCSB_B});
                    }
                    for (auto s : ms.acc_detsccs_[true_index - iw_succ.size()].second)
                    {
                      ranks.push_back({s, NCSB_S});
                    }
                  }
                  else 
                  {
                    for (unsigned i=0; i<acc_det_succ.size(); i++)
                    {
                      for (auto s : ms.acc_detsccs_[i].first)
                      {
                        ranks.push_back({s, NCSB_C}); 
                        ranks.push_back({s, NCSB_B});
                      }
                      for (auto s : ms.acc_detsccs_[i].second)
                      {
                        ranks.push_back({s, NCSB_S});
                      }
                    }
                  }

                  for (auto i : index)
                  {
                    auto succ = get_succ_track_CSB(reachable, letter, ranks, i, new_succ[0], new_succ[1], acc_det_succ, true_index - iw_succ.size());
                    //get_succ_track_CSB(reachable, letter, ranks, i, new_succ[1], acc_det_succ, true_index - iw_succ.size());
                    if (succ.empty())
                    {
                      new_succ[0].set_active_index(-2);
                      new_succ[1].set_active_index(-2);
                    }
                  }
                  
                }
                else
                {
                  // getSuccTrackToActive
                  std::vector<state_rank> ranks;
                  for (auto s : ms.acc_detsccs_[true_index - iw_succ.size()].first)
                  {
                    ranks.push_back({s, NCSB_C}); 
                    ranks.push_back({s, NCSB_B});
                  }
                  for (auto s : ms.acc_detsccs_[true_index - iw_succ.size()].second)
                  {
                    ranks.push_back({s, NCSB_S});
                  }

                  for (auto i : index)
                  {
                    // get successors for every deterministic component
                    auto succ = get_succ_active_CSB(std::set<unsigned>(ms.curr_reachable_.begin(), ms.curr_reachable_.end()), letter, ranks, i, new_succ[0], new_succ[1], acc_det_succ, true_index - iw_succ.size());

                    if (succ.empty())
                    {
                      new_succ[0].set_active_index(-2);
                      new_succ[1].set_active_index(-2);
                    }
                  }
                  
                }
              }

              else
              {
                  std::cerr << "ERROR" << std::endl;
                  throw std::runtime_error("nondeterministic accepting sccs not supported yet");
              }
            }

            else
            {
              // active component
              std::cerr << "Active" << std::endl;
              if (is_weakscc(scc_types_, index[0]))
              {
                iwa_done = true;
                // getSuccActive
                if ((not decomp_options_.tgba and ((not decomp_options_.merge_iwa and this->weaksccs_.size() + this->acc_detsccs_.size() > 1) or this->acc_detsccs_.size() > 0 or not ms.iw_break_set_.empty())) or not ms.iw_break_set_.empty())
                {
                  std::cerr << "MH succ active" << std::endl;
                  auto succ_active = mh.get_succ_active(reachable, reach_track, letter, index, ms.iw_break_set_);
                  for (auto succ_tt : succ_active.first)
                  {
                    if (decomp_options_.merge_iwa)
                      iw_succ[0] = std::vector<unsigned>(succ_tt.first.begin(), succ_tt.first.end());
                    else
                      iw_succ[true_index] = std::vector<unsigned>(succ_tt.first.begin(), succ_tt.first.end());
                  }
                  for (auto succ_at : succ_active.second)
                  {
                    if (decomp_options_.merge_iwa)
                      iw_succ[0] = std::vector<unsigned>(succ_at.first.first.begin(), succ_at.first.first.end());
                    else
                      iw_succ[true_index] = std::vector<unsigned>(succ_at.first.first.begin(), succ_at.first.first.end());
                    for (unsigned i=0; i<new_succ.size(); i++)
                    {
                      new_succ[i].set_iw_break_set(std::vector<unsigned>(succ_at.first.second.begin(), succ_at.first.second.end()));
                    }
                  }

                  if (succ_active.first.size() > 0)
                    active_type = false;
                }

                // getSuccTrackToActive 
                else 
                {
                  auto succ_track_to_active = mh.get_succ_track_to_active(reachable, reach_track, letter, index);
                  for (auto succ : succ_track_to_active) 
                  {
                    if (decomp_options_.merge_iwa)
                      iw_succ[0] = std::vector<unsigned>(succ.first.first.begin(), succ.first.first.end());
                    else
                      iw_succ[true_index] = std::vector<unsigned>(succ.first.first.begin(), succ.first.first.end());
                    new_succ[0].set_iw_break_set(std::vector<unsigned>(succ.first.second.begin(), succ.first.second.end()));
                    new_succ[1].set_iw_break_set(std::vector<unsigned>(succ.first.second.begin(), succ.first.second.end()));
                    //new_succ[1].set_active_index(-2);
                  }
                  active_type = true; 
                }
                
              }

              else if (is_accepting_detscc(scc_types_, index[0]))
              { 
                // getSuccActive 
                // currently sampled components
                std::cerr << "getSuccActive" << std::endl;
                std::set<unsigned> indices;
                for (auto s : ms.acc_detsccs_[true_index - iw_succ.size()].first)
                {
                  indices.insert(scc_info.scc_of(s));
                }
                for (auto s : ms.acc_detsccs_[true_index - iw_succ.size()].second)
                {
                  indices.insert(scc_info.scc_of(s));
                }

                if (decomp_options_.merge_det)
                {
                  std::set<unsigned> tmp;
                  for (auto s : reach_track)
                  {
                    if (indices.find(scc_info.scc_of(s)) != indices.end())
                      tmp.insert(s);
                  }
                  reach_track = tmp;
                }

                for (auto i : index)
                {
                  std::vector<state_rank> ranks;
                  for (auto s : ms.acc_detsccs_[true_index - iw_succ.size()].first)
                  {
                    ranks.push_back({s, NCSB_C}); 
                  }
                  for (auto s : ms.acc_detsccs_[true_index - iw_succ.size()].second)
                  {
                    ranks.push_back({s, NCSB_S});
                  }
                  for (auto s : ms.det_break_set_)
                  {
                    ranks.push_back({s, NCSB_B});
                  }

                  //if (indices.find(i) != indices.end() or ms.iw_break_set_.empty()) // TODO
                  //{
                    succ_det = get_succ_active_CSB((ms.det_break_set_.empty()) ? std::set<unsigned>(ms.curr_reachable_.begin(), ms.curr_reachable_.end()) : reach_track, letter, ranks, i, new_succ[0], new_succ[1], acc_det_succ, true_index - iw_succ.size());

                    // std::cerr << "Name: " << get_name(new_succ[0]) << std::endl;

                    if (succ_det.empty())
                    {
                      new_succ[0].set_active_index(-2); 
                      new_succ[1].set_active_index(-2);
                    } 
                  //}
                }
              }

              else
              {
                  std::cerr << "ERROR" << std::endl;
                  throw std::runtime_error("nondeterministic accepting sccs not supported yet");
              }
              
            }
            //std::cerr << "End" << std::endl;
          }
          

          for (unsigned i=0; i<new_succ.size(); i++)
          {
            if (new_succ[i].active_index_ != -2)
            {
              active_type = true;
              if (ms.iw_break_set_.size() == 0 and ms.det_break_set_.size() == 0)
                active_type = false;
              if (decomp_options_.tgba and ms.det_break_set_.size() == 0)
                active_type = false;
              if (not active_type){ 
                new_succ[i].set_active_index((indices[(orig_index + 1)%indices.size()]));
                
                // TODO !!!
                if (decomp_options_.tgba)
                {
                  // no active index for all inherently weak sccs
                  while (not is_accepting_detscc(scc_types_, new_succ[i].active_index_))
                  {
                    orig_index = (orig_index + 1) % indices.size();
                    new_succ[i].active_index_ = indices[orig_index];
                  }
                }
                else
                {
                  // no active index for nonaccepting scc
                  while (is_weakscc(scc_types_, new_succ[i].active_index_) and (not is_accepting_weakscc(scc_types_, new_succ[i].active_index_)))
                    new_succ[i].active_index_ = (new_succ[i].active_index_ + 1) % indices.size();
                }
              }
              else
                new_succ[i].set_active_index(active_index);
              new_succ[i].set_iw_sccs(iw_succ);

              new_succ[i].curr_reachable_ = std::vector<unsigned>(all_succ.begin(), all_succ.end());

              std::cerr << "New succ: " << get_name(new_succ[i]) << std::endl;
              if (std::find(all_states.begin(), all_states.end(), new_succ[i]) == all_states.end())
              {
                all_states.push_back(new_succ[i]);
                auto s = new_state(new_succ[i]);
              }

              auto p = rank2n_.emplace(new_succ[i], 0);
              if (active_type)
                res_->new_edge(top.second, p.first->second, letter);
              else 
                res_->new_edge(top.second, p.first->second, letter, {0});
              
              if (decomp_options_.tgba and ms.iw_break_set_.size() == 0)
                res_->new_edge(top.second, p.first->second, letter, {1});
            }
          }
        }
      }

      spot::print_hoa(std::cerr, res_);
      std::cerr << std::endl;
      
      return res_;

    }
  };

  spot::twa_graph_ptr
  complement_tnba(const spot::const_twa_graph_ptr &aut, spot::option_map &om, compl_decomp_options decomp_options)
  {
    const int trans_pruning = om.get(NUM_TRANS_PRUNING);
    // now we compute the simulator
    spot::const_twa_graph_ptr aut_reduced;
    std::vector<bdd> implications;
    spot::twa_graph_ptr aut_tmp = nullptr;
    if (om.get(USE_SIMULATION) > 0)
    {
      aut_tmp = spot::scc_filter(aut);
      auto aut2 = spot::simulation(aut_tmp, &implications, trans_pruning);
      aut_tmp = aut2;
    }
    if (aut_tmp)
      aut_reduced = aut_tmp;
    else
      aut_reduced = aut;
    spot::scc_info scc(aut_reduced, spot::scc_info_options::ALL);

    auto comp = cola::tnba_complement(aut_reduced, scc, om, implications, decomp_options);
    return comp.run();
  }
}
