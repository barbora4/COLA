#include "complement_mstate.hpp"

namespace cola
{
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
                                return na_sccs_ < other.na_sccs_;
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
        if (this->na_sccs_ != other.na_sccs_)
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
        for (auto v : na_sccs_)
        {
            for (auto s : v.reachable)
            {
                res ^= (res << 3) ^ s;
            }
            res ^= (res << 3) ^ v.i;
            for (auto s : v.O)
            {
                res ^= (res << 3) ^ s;
            }
            for (auto pr : v.f)
            {
                res ^= (res << 3) ^ (pr.first ^ pr.second);
            }
        }

        return res;
    }
}