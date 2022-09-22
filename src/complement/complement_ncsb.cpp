#include "complement_ncsb.hpp"

#include <algorithm>

namespace cola
{
  complement_mstate
  ncsb_compl::get_init_track()
  {
    complement_mstate mstate(scc_info_);

    unsigned orig_init = aut_->get_init_state_number();
    mstate.curr_reachable_.push_back(orig_init);

    mstate.acc_detsccs_.push_back({std::vector<unsigned>(), std::vector<unsigned>()});

    return mstate;
  }

  complement_mstate
  ncsb_compl::get_init_active()
  {
    complement_mstate mstate(scc_info_);

    unsigned orig_init = aut_->get_init_state_number();
    mstate.curr_reachable_.push_back(orig_init);

    mstate.acc_detsccs_.push_back({std::vector<unsigned>(1, orig_init), std::vector<unsigned>()});
    mstate.det_break_set_ = std::vector<unsigned>(1, orig_init); 

    return mstate;
  }

  std::vector<std::pair<complement_mstate, bool>>
  ncsb_compl::get_succ_track(complement_mstate mstate, bdd symbol)
  {
    std::vector<std::pair<complement_mstate, bool>> succ;

    if (get_succ_acc_trans_scc(mstate.acc_detsccs_[true_index_].second, symbol).empty())
    {
      complement_mstate succ_state(scc_info_);

      auto S_succ = get_all_successors_in_scc(mstate.acc_detsccs_[true_index_].second, symbol);
      std::vector<unsigned> S_prime = std::vector<unsigned>(S_succ.begin(), S_succ.end());

      auto C_succ = get_all_successors_in_scc(mstate.acc_detsccs_[true_index_].first, symbol);
      std::vector<unsigned> C_prime;
      std::set_difference(C_succ.begin(), C_succ.end(), S_succ.begin(), S_succ.end(), std::inserter(C_prime, C_prime.begin()));

      succ_state.acc_detsccs_.push_back({C_prime, S_prime});
      succ.push_back({succ_state, false});
    }

    return succ;
  }

  std::vector<std::pair<complement_mstate, bool>> ncsb_compl::get_succ_track_to_active(complement_mstate mstate, bdd symbol)
  {
    std::vector<std::pair<complement_mstate, bool>> succ;

    if (get_succ_acc_trans_scc(mstate.acc_detsccs_[true_index_].second, symbol).empty())
    {
      complement_mstate succ_state(scc_info_);

      auto S_succ = get_all_successors_in_scc_same(mstate.acc_detsccs_[true_index_].second, symbol);
      std::vector<unsigned> S_prime = std::vector<unsigned>(S_succ.begin(), S_succ.end());

      auto C_succ = get_all_successors_in_scc(mstate.curr_reachable_, symbol);
      std::vector<unsigned> C_prime;
      std::set_difference(C_succ.begin(), C_succ.end(), S_succ.begin(), S_succ.end(), std::inserter(C_prime, C_prime.begin()));

      succ_state.acc_detsccs_.push_back({C_prime, S_prime});
      succ_state.det_break_set_ = C_prime;
      succ.push_back({succ_state, false});
    }

    return succ;
  }

  std::vector<std::pair<complement_mstate, bool>>
  ncsb_compl::get_succ_active(complement_mstate mstate, bdd symbol)
  {
    std::vector<std::pair<complement_mstate, bool>> succ;

    if (mstate.det_break_set_.empty())
    {
      // empty break set -> return TT and switch to other scc
      if ((decomp_options_.merge_det or mstate.acc_detsccs_.size() == 1) and mstate.iw_sccs_.size() == 0 and mstate.na_sccs_.size() == 0)
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

    if (get_succ_acc_trans_scc(mstate.acc_detsccs_[true_index_].second, symbol).empty())
    {
      complement_mstate succ1(scc_info_);
      complement_mstate succ2(scc_info_);

      std::set<unsigned> S_prime = get_all_successors_in_scc_same(mstate.acc_detsccs_[true_index_].second, symbol);

      std::set<unsigned> C_succ = get_all_successors_in_scc(mstate.acc_detsccs_[true_index_].first, symbol);
      std::set<unsigned> C_prime;
      std::set_difference(C_succ.begin(), C_succ.end(), S_prime.begin(), S_prime.end(), std::inserter(C_prime, C_prime.begin()));

      if (mstate.det_break_set_.empty())
      {
        succ1.acc_detsccs_.push_back({std::vector<unsigned>(C_prime.begin(), C_prime.end()), std::vector<unsigned>(S_prime.begin(), S_prime.end())});
        succ.push_back({succ1, true});
        return succ;
      }

      std::set<unsigned> B_succ = get_all_successors_in_scc_same(mstate.det_break_set_, symbol);
      std::set<unsigned> B_prime;
      std::set_intersection(B_succ.begin(), B_succ.end(), C_prime.begin(), C_prime.end(), std::inserter(B_prime, B_prime.begin()));

      succ1.acc_detsccs_.push_back({std::vector<unsigned>(C_prime.begin(), C_prime.end()), std::vector<unsigned>(S_prime.begin(), S_prime.end())});
      succ1.det_break_set_ = std::vector<unsigned>(B_prime.begin(), B_prime.end());
      succ.push_back({succ1, false});

      if (std::any_of(B_prime.begin(), B_prime.end(), [this](unsigned state)
                      { return is_accepting_[state]; }) or
          not get_succ_acc_trans_scc(mstate.det_break_set_, symbol).empty())
      {
        return succ;
      }

      std::set<unsigned> S_prime2;
      std::set_union(S_prime.begin(), S_prime.end(), B_prime.begin(), B_prime.end(), std::inserter(S_prime2, S_prime2.begin()));

      std::set<unsigned> C_prime2;
      std::set_difference(C_prime.begin(), C_prime.end(), S_prime2.begin(), S_prime2.end(), std::inserter(C_prime2, C_prime2.begin()));

      succ2.acc_detsccs_.push_back({std::vector<unsigned>(C_prime2.begin(), C_prime2.end()), std::vector<unsigned>(S_prime2.begin(), S_prime2.end())});
      succ2.det_break_set_ = succ2.acc_detsccs_.back().first;
      succ.push_back({succ2, false});
    }

    return succ;
  }
}