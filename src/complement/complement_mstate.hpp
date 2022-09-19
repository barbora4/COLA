#pragma once

#include "cola.hpp"
#include "types.hpp"
#include "rankings.hpp"

namespace cola
{
    class complement_mstate
    {
    public:
        complement_mstate(spot::scc_info &si) : si_(si){}

        complement_mstate(const complement_mstate &other)
            : si_(other.si_)
        {
            this->set_iw_sccs(other.iw_sccs_);
            this->set_iw_break_set(other.iw_break_set_);
            this->det_break_set_ = other.det_break_set_;
            this->set_active_index(other.active_index_);
            this->set_acc_detsccs(other.acc_detsccs_);
            this->curr_reachable_ = other.curr_reachable_;

            this->na_sccs_ = other.na_sccs_;
        }

        std::set<unsigned>
        get_reach_set() const;

        bool operator<(const complement_mstate &other) const;
        bool operator==(const complement_mstate &other) const;

        complement_mstate &
        operator=(const complement_mstate &other)
        {
            this->si_ = other.si_;

            this->set_iw_sccs(other.iw_sccs_);
            this->set_iw_break_set(other.iw_break_set_);
            this->det_break_set_ = other.det_break_set_;
            this->set_active_index(other.active_index_);
            this->set_acc_detsccs(other.acc_detsccs_);
            this->curr_reachable_ = other.curr_reachable_;

            this->na_sccs_ = other.na_sccs_;

            return *this;
        }

        size_t hash() const;

        // SCC information
        spot::scc_info &si_;

        std::vector<std::vector<unsigned>> iw_sccs_;
        std::vector<std::pair<std::vector<unsigned>, std::vector<unsigned>>> acc_detsccs_;
        std::vector<unsigned> iw_break_set_;
        std::vector<unsigned> det_break_set_;
        int active_index_ = 0;
        std::vector<unsigned> curr_reachable_;

        std::vector<rank_state> na_sccs_;

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

   
}