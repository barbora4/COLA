#pragma once

#include <set>
#include <vector>

#include "cola.hpp"

#include <spot/misc/hashfunc.hh>
#include <spot/twaalgos/isdet.hh>
#include <spot/twaalgos/sccinfo.hh>
#include <spot/twaalgos/isunamb.hh>
#include <spot/twaalgos/product.hh>
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

// Miyano-Hayashi complementation for IW components
namespace cola 
{
  class mh_complement
  {
  private:
    // source automaton
    const spot::const_twa_graph_ptr aut_;
    spot::scc_info scc_info_;
    std::string scc_types_;
    compl_decomp_options decomp_options_;
    std::vector<std::pair<unsigned, unsigned>> dir_sim_;
    std::vector<std::set<int>> reachable_vector_;
  
  public:
    mh_complement(const spot::const_twa_graph_ptr &aut, spot::scc_info &scc_info, std::string scc_types, compl_decomp_options decomp_options, std::vector<std::pair<unsigned, unsigned>> dir_sim) : aut_(aut), scc_info_(scc_info), scc_types_(scc_types), decomp_options_(decomp_options), dir_sim_(dir_sim) {}
  
    std::vector<std::pair<std::set<unsigned>, unsigned>> get_succ_track(std::set<unsigned> reachable, std::set<unsigned> reach_in_scc, bdd symbol, std::vector<unsigned> scc_index);
    std::vector<std::pair<std::pair<std::set<unsigned>, std::set<unsigned>>, unsigned>> get_succ_track_to_active(std::set<unsigned> reachable, std::set<unsigned> reach_in_scc, bdd symbol, std::vector<unsigned> scc_index);
    std::pair<std::vector<std::pair<std::set<unsigned>, unsigned>>, std::vector<std::pair<std::pair<std::set<unsigned>, std::set<unsigned>>, unsigned>>> get_succ_active(std::set<unsigned> reachable, std::set<unsigned> reach_in_scc, bdd symbol, std::vector<unsigned> scc_index, std::vector<unsigned> break_set, bool one_scc=false);
    std::set<unsigned> get_all_successors(std::set<unsigned> current_states, bdd symbol);
    std::vector<std::set<int>>  get_reachable_vector();
    std::set<int> reachable_vertices(std::vector<std::vector<int>> list, std::set<int> from);
  };
}