// Copyright (C) 2019-2020  The Seminator Authors
// Copyright (C) 2022  The COLA Authors
//
// This file is a part of cola, a tool for complementation and determinization
// of omega automata.
//
// cola is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// cola is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "config.h"

#include "cola.hpp"
#include "composer.hpp"
#include "optimizer.hpp"
#include "decomposer.hpp"
#include "simulation.hpp"
// #include "postproc.hpp"

#include <unistd.h>
#include <fstream>
#include <ctime>
#include <string>
#include <sstream>

#include <spot/twaalgos/simulation.hh>
#include <spot/parseaut/public.hh>
#include <spot/twaalgos/isunamb.hh>
#include <spot/twaalgos/isdet.hh>
#include <spot/twaalgos/degen.hh>
#include <spot/twaalgos/hoa.hh>
#include <spot/twaalgos/sccfilter.hh>
#include <spot/twaalgos/complement.hh>
#include <spot/twaalgos/minimize.hh>
#include <spot/twaalgos/totgba.hh>
#include <spot/twaalgos/determinize.hh>
#include <spot/twaalgos/zlktree.hh>
#include <spot/twaalgos/dualize.hh>
#include <spot/twaalgos/word.hh>
#include <spot/misc/version.hh>

void print_usage(std::ostream &os)
{
  os << "Usage: cola [OPTION...] [FILENAME...]\n";
}

void print_help()
{
  print_usage(std::cout);
  std::cout <<
      R"(The tool obtains the equivalent deterministic automaton or complement automaton of input Buchi automaton.

By default, it reads a (generalized) Büchi automaton from standard input
and converts it into deterministic automata.

Input options:
    -f FILENAME reads the input from FILENAME instead of stdin
    --algo=[cola|dc|iar|comp|ncsb|congr]
            Use determinization or complementation algorithms (comp or ncsb) to obtain the output
              cola     Default setting for determinization (--algo=dc --simulation --use-scc --stutter --parity)
              dc       Divide-and-Conquer determinization based on SCC decomposition
              iar      Specialized algorithm for limit-deterministic Buchi automata in TACAS'17 paper by Esparza et al.
              comp     Complementation algorithm based on SCC decomposition
              ncsb     NCSB complementation variants for limit deterministic Buchi automata 
              congr    Congruence-based algorithm for containment checking (with --contain=...)
    --type 
            Output the type of the input Buchi automaton: deterministic, limit-deterministic, elevator, unambiguous or none of them
    --print-scc
            Output the information about the SCCs in the input NBA
    --contain=[FILENAME]
            Test whether the language of the input contains the language of [FILENAME] 

Output options:
    --verbose=[INT] Output verbose level (0 = minimal level, 1 = meduim level, 2 = debug level)
    -o FILENAME     Write the output to FILENAME instead of stdout
    --generic       Output the automaton with Emenson-Lei acceptance condition
    --rabin         Output the automaton with generalized Rabin condition
    --parity        Output the automaton with Pairty acceptance condition (Default)
    --complement    Output the complement Buchi automaton of the input (default with --algo=comp)


Optimizations:
    --simulation          Use direct simulation for determinization/complementation
    --stutter             Use stutter invariance for determinization
    --use-scc             Use SCC information for macrostates merging
    --more-acc-egdes      Enumerate elementary cycles for obtaining more accepting egdes 
    --trans-pruning=[INT] Number to limit the transition pruning in simulation (default=512) 
    --unambiguous         Check whether the input is unambiguous and use this fact in determinization
    --rerank              Rearrange the labelling for NAC-states

Decomposition-based complementation:
    --merge-iwa           Merge all IWA SCCs
    --merge-det           Merge all deterministic SCCs
    --tgba                Outputs a TGBA with two colours
    --iw-sim              Simulation on IWA SCCs
    --det-sim             Simulation on DET SCCs
    --scc-compl           Complementation for each SCC separately
    --scc-high            SCC compl with high postprocessing before intersection
    --no-sat              No saturation of accepting states/transitions
    --dataflow            Data flow analysis in rank-based complementation
    --debug               Print state names
    --print-hoa           Print preprocessed automaton in HOA

Pre- and Post-processing:
    --preprocess=[0|1|2|3]       Level for simplifying the input automaton (default=1)
    --postprocess-det[=0|1|2|3]  Level for simplifying the output of the determinization (default=1)
    --num-states=[INT]           Simplify the output with number of states less than INT (default=30000)

Miscellaneous options:
  -h, --help    Print this help
  --version     Print program version
)";
}

void check_cout()
{
  std::cout.flush();
  if (!std::cout)
  {
    std::cerr << "cola: error writing to standard output\n";
    exit(2);
  }
}

unsigned
parse_int(const std::string &arg)
{
  unsigned result;
  // obtain the substring after '='
  std::size_t idx = arg.find('=');
  //std::cout << "Index of = : " << idx << std::endl;
  std::string number = arg.substr(idx + 1, arg.length());
  std::istringstream iss(number);
  iss >> result;
  return result;
}

// determinization
enum determinize_algo
{
  NoDeterminize = 0,
  DC,
  IAR // induction appearance record
};

// determinization
enum complement_algo
{
  NoComplement = 0,
  COMP,
  NCSB,
  CONGR
};

enum postprocess_level
{
  None = 0,
  Low,
  Medium,
  High
};

//transition-based automata
enum output_aut_type
{
  Generic = 0,
  Parity,
  Rabin,
  Buchi,
  GeneralizedBuchi
};

// We may provide multiple algorithms for comparison
spot::twa_graph_ptr
to_deterministic(spot::twa_graph_ptr aut, spot::option_map &om, unsigned aut_type, determinize_algo algo)
{
  // determinization
  spot::twa_graph_ptr res;
  if (algo == DC)
  {
    if (aut_type & INHERENTLY_WEAK)
      res = cola::determinize_twba(aut, om);
    else
      res = cola::determinize_tnba(aut, om);
  }
  else if (algo == IAR)
  {
    res = cola::determinize_tldba(aut, om);
  }
  else
  {
    res = cola::determinize_tnba(aut, om);
  }
  return res;
}

void output_input_type(spot::twa_graph_ptr aut)
{
  bool type = false;
  if (spot::is_deterministic(aut))
  {
    type = true;
    std::cout << "deterministic" << std::endl;
  }
  if (spot::is_semi_deterministic(aut))
  {
    type = true;
    std::cout << "limit-deterministic" << std::endl;
  }
  if (cola::is_elevator_automaton(aut))
  {
    std::cout << "elevator" << std::endl;
  }
  if (cola::is_weak_automaton(aut))
  {
    std::cout << "inherently weak" << std::endl;
  }
  if (spot::is_unambiguous(aut))
  {
    std::cout << "unambiguous" << std::endl;
  }
  if (!type)
  {
    std::cout << "nondeterministic" << std::endl;
  }
}

void output_scc_info(spot::twa_graph_ptr aut)
{
  // strengther
  spot::scc_info si(aut, spot::scc_info_options::ALL);
  unsigned num_iwcs = 0;
  unsigned num_acc_iwcs = 0;
  unsigned num_iwcs_states = 0;
  unsigned num_max_iwcs_states = 0;
  unsigned num_acciwcs_states = 0;
  unsigned num_max_acciwcs_states = 0;
  unsigned num_dacs = 0;
  unsigned num_dacs_states = 0;
  unsigned num_max_dacs_states = 0;
  unsigned num_nacs = 0;
  unsigned num_nacs_states = 0;
  unsigned num_max_nacs_states = 0;

  std::string types = cola::get_scc_types(si);
  for (unsigned sc = 0; sc < si.scc_count(); sc++)
  {
    unsigned num = si.states_of(sc).size();
    if (cola::is_weakscc(types, sc))
    {
      num_iwcs_states += num;
      num_iwcs++;
      num_max_iwcs_states = std::max(num_max_iwcs_states, num);
    }
    if (cola::is_accepting_weakscc(types, sc))
    {
      num_acciwcs_states += num;
      num_acc_iwcs++;
      num_max_acciwcs_states = std::max(num_max_acciwcs_states, num);
    }

    if (cola::is_accepting_detscc(types, sc))
    {
      num_dacs_states += num;
      num_dacs++;
      num_max_dacs_states = std::max(num_max_dacs_states, num);
    }
    if (cola::is_accepting_nondetscc(types, sc))
    {
      num_nacs_states += num;
      num_nacs++;
      num_max_nacs_states = std::max(num_max_nacs_states, num);
    }
  }
  std::cout << "Number of IWCs: " << num_iwcs << " with " << num_iwcs_states << " states, in which max IWC with " << num_max_iwcs_states << " states\n";
  std::cout << "Number of ACC_IWCs: " << num_acc_iwcs << " with " << num_acciwcs_states << " states, in which max IWC with " << num_max_acciwcs_states << " states\n";
  std::cout << "Number of DACs: " << num_dacs << " with " << num_dacs_states << " states, in which max DAC with " << num_max_dacs_states << " states\n";
  std::cout << "Number of NACs: " << num_nacs << " with " << num_nacs_states << " states, in which max NAC with " << num_max_nacs_states << " states\n";
}

int main(int argc, char *argv[])
{
  // Declaration for input options. The rest is in cola.hpp
  // as they need to be included in other files.
  bool cd_check = false;
  bool high = false;
  std::vector<std::string> path_to_files;

  spot::option_map om;
  // default setting
  om.set(USE_SIMULATION, 0);
  om.set(USE_STUTTER, 0);
  om.set(USE_UNAMBIGUITY, 0);
  om.set(USE_SCC_INFO, 0);
  om.set(VERBOSE_LEVEL, 0);
  om.set(USE_DELAYED_SIMULATION, 0);
  om.set(MORE_ACC_EDGES, 0);
  om.set(NUM_TRANS_PRUNING, 512);

  om.set(SCC_REACH_MEMORY_LIMIT, 0);
  om.set(NUM_SCC_LIMIT_MERGER, 0);
  om.set(MSTATE_REARRANGE, 0);

  determinize_algo determinize = NoDeterminize;

  complement_algo complement = NoComplement;

  // options
  bool use_simulation = false;
  //bool merge_transitions = false;
  bool debug = false;
  bool aut_type = false;
  bool use_unambiguous = false;
  bool use_stutter = false;
  bool decompose = false;
  bool use_acd = false;
  bool print_scc = false;
  bool comp = false;
  bool contain = false;
  bool congr = false;
  std::string file_to_contain;

  compl_decomp_options decomp_options;

  postprocess_level preprocess = Low;
  postprocess_level post_process = Low; 
  bool use_scc = false;
  unsigned num_post = 30000;

  output_aut_type output_type = Generic;

  std::string output_filename = "";

  for (int i = 1; i < argc; i++)
  {
    std::string arg = argv[i];
    if (arg.find("--preprocess=") != std::string::npos)
    {
      unsigned level = parse_int(arg);
      if (level == 0)
      {
        preprocess = None;
      }
      else if (level == 1)
      {
        preprocess = Low;
      }
      else if (level == 2)
      {
        preprocess = Medium;
      }
      else if (level == 3)
      {
        preprocess = High;
      }
    }
    else if (arg == "--print-scc")
    {
      print_scc = true;
    }
    else if (arg == "--postprocess-det=0")
      post_process = None;
    else if (arg == "--postprocess-det=1")
      post_process = Low;
    else if (arg == "--postprocess-det=2")
      post_process = Medium;
    else if (arg == "--postprocess-det=3")
      post_process = High;
    else if (arg == "--generic")
    {
      output_type = Generic;
    }
    else if (arg == "--parity")
    {
      output_type = Parity;
    }
    else if (arg == "--rabin")
    {
      output_type = Rabin;
    }
    else if (arg == "--complement")
    {
      comp = true;
      output_type = Buchi;
    }
    else if (arg == "--simulation")
    {
      use_simulation = true;
      om.set(USE_SIMULATION, 1);
    }else if (arg == "--rerank")
    {
      om.set(MSTATE_REARRANGE, 1);
    }
    else if (arg.find("--trans-pruning=") != std::string::npos)
    {
      int trans_pruning = parse_int(arg);
      om.set(NUM_TRANS_PRUNING, trans_pruning);
    }
    else if (arg == "--delayed-sim")
    {
      om.set(USE_DELAYED_SIMULATION, 1);
    }
    else if (arg == "--use-scc")
    {
      use_scc = true;
      om.set(USE_SCC_INFO, 1);
    }
    else if (arg == "--more-acc-edges")
    {
      om.set(MORE_ACC_EDGES, 1);
    }
    else if (arg == "--decompose")
    {
      decompose = true;
      om.set(NUM_NBA_DECOMPOSED, -1);
    }
    else if (arg.find("--decompose=") != std::string::npos)
    {
      decompose = true;
      unsigned num_scc = parse_int(arg);
      om.set(NUM_NBA_DECOMPOSED, num_scc);
    }
    else if (arg == "--acd")
    {
      use_acd = true;
    }
    // Prefered output
    else if (arg == "--d")
      debug = true;
    else if (arg == "--type")
      aut_type = true;
    else if (arg == "--unambiguous")
    {
      use_unambiguous = true;
      om.set(USE_UNAMBIGUITY, 1);
    }
    else if (arg == "--stutter")
    {
      use_stutter = true;
      om.set(USE_STUTTER, 1);
    }
    else if (arg == "--algo=dc")
      determinize = DC;
    else if (arg == "--algo=iar")
      determinize = IAR;
    else if (arg == "--algo=cola")
    {
      determinize = DC;
      // default settings
      om.set(USE_SIMULATION, 1);
      om.set(USE_SCC_INFO, 1);
      om.set(USE_STUTTER, 1);
      use_acd = true;
      output_type = Parity;
    }
    else if (arg == "--algo=comp")
    {
      comp = true;
      complement = COMP;
    }
    else if (arg == "--algo=ncsb")
    {
      comp = true;
      complement = NCSB;
    }
    else if (arg == "--algo=congr")
    {
      complement = CONGR;
    }
    else if (arg == "--merge-iwa")
    {
      decomp_options.merge_iwa = true;
    }
    else if (arg == "--merge-det")
    {
      decomp_options.merge_det = true;
    }
    else if (arg == "--tgba")
    {
      decomp_options.tgba = true;
    }
    else if (arg == "--iw-sim")
    {
      decomp_options.iw_sim = true;
    }
    else if (arg == "--det-sim")
    {
      decomp_options.det_sim = true;
    }
    else if (arg == "--scc-compl")
    {
      decomp_options.scc_compl = true;
    }
    else if (arg == "--scc-high")
    {
      decomp_options.scc_compl_high = true;
    }
    else if (arg == "--no-sat")
    {
      decomp_options.sat = false;
    }
    else if (arg == "--dataflow")
    {
      decomp_options.dataflow = true;
    }
    else if (arg == "--debug")
    {
      decomp_options.debug = true;
    }
    else if (arg == "--print-hoa")
    {
      decomp_options.print_hoa = true;
    }
    else if (arg == "-f")
    {
      if (argc < i + 1)
      {
        std::cerr << "cola: Option -f requires an argument.\n";
        return 1;
      }
      else
      {
        path_to_files.emplace_back(argv[i + 1]);
        i++;
      }
    }
    else if (arg == "-o")
    {
      if (argc < i + 1)
      {
        std::cerr << "cola: Option -o requires an argument.\n";
        return 1;
      }
      else
      {
        std::string str(argv[i + 1]);
        output_filename = str;
        i++;
      }
    }
    else if (arg.find("--contain=") != std::string::npos)
    {
      contain = true;
      std::size_t idx = arg.find('=');
      file_to_contain = arg.substr(idx + 1, arg.length());
    }
    else if (arg.find("--num-states=") != std::string::npos)
    {
      // obtain the substring after '='
      num_post = parse_int(arg);
      //std::cout << "Input number : " << num_post << std::endl;
    }
    else if (arg.find("--verbose=") != std::string::npos)
    {
      om.set(VERBOSE_LEVEL, parse_int(arg));
      decomp_options.dir_sim = false;
    }
    else if (arg.find("--scc-mem-limit=") != std::string::npos)
    {
      om.set(SCC_REACH_MEMORY_LIMIT, parse_int(arg));
    }
    else if (arg.find("--scc-num-limit=") != std::string::npos)
    {
      om.set(NUM_SCC_LIMIT_MERGER, parse_int(arg));
    }
    else if ((arg == "--help") || (arg == "-h"))
    {
      print_help();
      check_cout();
      return 0;
    }
    else if (arg == "--version")
    {
      std::cout << "cola " PACKAGE_VERSION
                   " (using Spot "
                << spot::version() << ")\n\n"
                                      "Copyright (C) 2020  The cola Authors.\n"
                                      "License GPLv3+: GNU GPL version 3 or later"
                                      " <http://gnu.org/licenses/gpl.html>.\n"
                                      "This is free software: you are free to change "
                                      "and redistribute it.\n"
                                      "There is NO WARRANTY, to the extent permitted by law.\n"
                << std::flush;
      return 0;
    }
    // Detection of unsupported options
    else if (arg[0] == '-')
    {
      std::cerr << "cola: Unsupported option " << arg << '\n';
      return 2;
    }
    else
    {
      path_to_files.emplace_back(argv[i]);
    }
  }

  if (path_to_files.empty())
  {
    if (isatty(STDIN_FILENO))
    {
      std::cerr << "cola: No automaton to process? "
                   "Run 'cola --help' for more help.\n";
      print_usage(std::cerr);
      return 1;
    }
    else
    {
      // Process stdin by default.
      path_to_files.emplace_back("-");
    }
  }

  auto dict = spot::make_bdd_dict();

  spot::twa_graph_ptr aut_to_contain = nullptr;
  // contain
  if (contain)
  {
    spot::automaton_stream_parser parser(file_to_contain);
    spot::parsed_aut_ptr parsed_aut = parser.parse(dict);

    if (parsed_aut->format_errors(std::cerr))
    {
      std::runtime_error("File " + file_to_contain + " is not in valid HOA format");
      return 1;
    }
    aut_to_contain = parsed_aut->aut;
  }

  for (std::string &path_to_file : path_to_files)
  {
    if (om.get(VERBOSE_LEVEL))
      std::cout << "File: " << path_to_file << " Algo: " << determinize << std::endl;
    spot::automaton_stream_parser parser(path_to_file);

    for (;;)
    {
      spot::parsed_aut_ptr parsed_aut = parser.parse(dict);

      if (parsed_aut->format_errors(std::cerr))
        return 1;

      // input automata
      spot::twa_graph_ptr aut = parsed_aut->aut;

      if (!aut)
        break;

      // Check if input is TGBA
      if (aut->acc().is_generalized_buchi())
      {
        aut = spot::degeneralize_tba(aut);
      }

      if (!aut->acc().is_buchi())
      {
        std::cerr << "cola requires Buchi condition on input.\n";
        return 1;
      }

      if (aut_type)
      {
        output_input_type(aut);
        break;
      }

      if (print_scc)
      {
        output_scc_info(aut);
        break;
      }

      if (om.get(MORE_ACC_EDGES) > 0)
      {
        const unsigned num = 200;
        // strengther
        spot::scc_info si(aut, spot::scc_info_options::ALL);
        cola::edge_strengther e_strengther(aut, si, 200);
        for (unsigned sc = 0; sc < si.scc_count(); sc++)
        {
          if (si.is_accepting_scc(sc))
          {
            e_strengther.fix_scc(sc);
          }
        }
      }

      //1. preprocess
      clock_t c_start = clock();
      unsigned aut_type = NONDETERMINISTIC;
      if (cola::is_weak_automaton(aut))
      {
        aut_type |= INHERENTLY_WEAK;
      }
      if (spot::is_semi_deterministic(aut))
      {
        aut_type |= LIMIT_DETERMINISTIC;
      }
      if (cola::is_elevator_automaton(aut))
      {
        aut_type |= ELEVATOR;
      }

      // preprocessing for the input.
      if (preprocess)
      {
        spot::postprocessor preprocessor;
        // only a very low level of preprocessing is allowed
        if (preprocess == Low)
          preprocessor.set_level(spot::postprocessor::Low);
        else if (preprocess == Medium)
          preprocessor.set_level(spot::postprocessor::Medium);
        else if (preprocess == High)
          preprocessor.set_level(spot::postprocessor::High);

        if (decomp_options.scc_compl)
          preprocessor.set_type(spot::postprocessor::Buchi); 
        aut = preprocessor.run(aut);
      }

      if (om.get(VERBOSE_LEVEL) >= 2)
      {
        cola::output_file(aut, "sim_aut.hoa");
        std::cout << "Output processed automaton (" << aut->num_states() << ", " << aut->num_edges() << ") to sim_aut.hoa\n";
      }
      clock_t c_end = clock();
      if (om.get(VERBOSE_LEVEL) > 0)
      {
        std::cout << "Done for preprocessing the input automaton in " << 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC << " ms..." << std::endl;
      }
      if (aut->acc().is_all())
      {
        // trivial acceptance condition
        aut = spot::minimize_monitor(aut);
      }
      if (!spot::is_deterministic(aut) && determinize)
      {
        if (decompose && aut->acc().is_buchi() && !spot::is_deterministic(aut))
        {
          cola::decomposer nba_decomposer(aut, om);
          std::vector<spot::twa_graph_ptr> subnbas = nba_decomposer.run();
          std::vector<spot::twa_graph_ptr> dpas;
          for (unsigned i = 0; i < subnbas.size(); i++)
          {
            spot::twa_graph_ptr dpa = to_deterministic(subnbas[i], om, aut_type, determinize);
            dpas.push_back(dpa);
          }
          cola::composer dpa_composer(dpas, om);
          aut = dpa_composer.run();
        }
        else if (aut->acc().is_buchi())
        {
          spot::twa_graph_ptr res = nullptr;
          c_start = clock();
          res = to_deterministic(aut, om, aut_type, determinize);
          c_end = clock();
          if (om.get(VERBOSE_LEVEL) > 0)
          {
            std::cout << "Done for determinizing the input automaton in " << 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC << " ms..." << std::endl;
          }
          aut = res;
          // spot::print_hoa(std::cout, aut); 
          // std::cout << std::endl;
        }
      }
      if (complement && !determinize)
      {
        if (complement == CONGR && contain)
        {
          if (!aut_to_contain)
            std::cout << "Contained" << std::endl;
          cola::congr_contain(aut, aut_to_contain, om);
        }
        else if (complement == COMP)
        {
          aut = cola::complement_tnba(aut, om, decomp_options);
          output_type = GeneralizedBuchi; 
        }
        else
        {
          // set NCSB algorithm later
          aut = cola::complement_tnba(aut, om, decomp_options);
        }
      }else if (comp && determinize)
      {
        aut = spot::dualize(aut);
      }
      const char *opts = nullptr;
      aut->merge_edges();
      if (om.get(VERBOSE_LEVEL) > 0) 
      {
        std::cout << "Number of (states, transitions, colors) in the result automaton: ("
                  << aut->num_states() << "," << aut->num_edges() << "," << aut->num_sets() << ")" << std::endl;
        spot::print_hoa(std::cout, aut);
        std::cout << std::endl;
      }
      // postprocessing, remove dead states
      //aut->purge_unreachable_states();
      if (post_process != None && !decompose)
      {
        clock_t c_start = clock();
        if (aut->acc().is_all())
        {
          aut = spot::minimize_monitor(aut);
        }
        else if (output_type != Buchi and output_type != GeneralizedBuchi)//if (aut->num_states() < num_post)
        {
          spot::postprocessor p;
          if (output_type == Parity)
          {
            if (use_acd)
            {
              p.set_type(spot::postprocessor::Generic);
            }
            else
            {
              p.set_type(spot::postprocessor::Parity);
            }
          }
          else if (output_type == Generic || output_type == Rabin)
          {
            p.set_type(spot::postprocessor::Generic);
          }
          p.set_pref(spot::postprocessor::Deterministic);
          // set postprocess level
          if (post_process == Low)
          {
            p.set_level(spot::postprocessor::Low);
          }
          else if (post_process == Medium)
          {
            p.set_level(spot::postprocessor::Medium);
          }
          else if (post_process == High)
          {
            p.set_level(spot::postprocessor::High);
          }
          aut = p.run(aut);
        }
        if (output_type == Rabin)
        {
          aut = spot::to_generalized_rabin(aut, true);
        }
        else if (output_type == Parity && use_acd)
        {
          // call the alternating cycle decomposition to translate our rabin automaton
          // to parity automaton
          aut = spot::acd_transform(aut);
        }
        // now post processing again since we may not do postprocessing above
        {
          spot::postprocessor p;
          if (post_process == Low)
          {
            p.set_level(spot::postprocessor::Low);
          }
          else if (post_process == Medium)
          {
            p.set_level(spot::postprocessor::Medium);
          }
          else if (post_process == High)
          {
            p.set_level(spot::postprocessor::High);
          }
          if (output_type != Buchi and output_type != GeneralizedBuchi) p.set_pref(spot::postprocessor::Deterministic);
          if (output_type == Generic)
          {
            p.set_type(spot::postprocessor::Generic);
          }
          else if (output_type == Parity)
          {
            p.set_type(spot::postprocessor::Parity);
          }
          else if (output_type == Buchi)
          {
            p.set_type(spot::postprocessor::Buchi);
          }
          else if (output_type == GeneralizedBuchi)
          {
            p.set_type(spot::postprocessor::GeneralizedBuchi);
          }
          aut = p.run(aut);
        }
        clock_t c_end = clock();
        if (om.get(VERBOSE_LEVEL) > 0)
          std::cout << "Done for postprocessing the result automaton in " << 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC << " ms..." << std::endl;
      }

      if (contain)
      {
        //now check whether the output is complementation
        if (!aut_to_contain)
        {
          std::cout << "Contained" << std::endl;
          break;
        }
        std::stringstream ss;
        bool has_counterexample = false;
        if (comp)
        {
          spot::twa_word_ptr word = aut->intersecting_word(aut_to_contain);
          if (word != nullptr)
          {
            ss << (*word);
            has_counterexample = true;
          }
        }
        else
        {
          // not complement, now the automaton should be determinized
          aut = spot::complement(aut);
          spot::twa_word_ptr word = aut->intersecting_word(aut_to_contain);
          if (word != nullptr)
          {
            ss << (*word);
            has_counterexample = true;
          }
        }
        if (!has_counterexample)
        {
          std::cout << "Contained" << std::endl;
        }
        else
        {
          std::cout << "Not contained: " << ss.str() << std::endl;
        }
      }
      else if (output_filename != "")
      {
        cola::output_file(aut, output_filename.c_str());
      }
      else
      {
        spot::print_hoa(std::cout, aut, opts);
        std::cout << "\n";
      }
    }
  }

  check_cout();

  return 0;
}
