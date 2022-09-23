#include "rankings.hpp"

#include <algorithm>

namespace cola
{
    std::string ranking::get_name() 
    {
        std::string name;

        for (auto it=this->begin(); it!=this->end(); it++)
        {
            name += std::to_string(std::get<0>(*it)) + ":" + std::to_string(std::get<1>(*it)) + ",";
        }

        return name;
    }
    
    static bool compare_ranks(std::tuple<int, int, bool> first, std::tuple<int, int, bool> second)
    {
        return (std::get<1>(first) < std::get<1>(second));
    }

    std::vector<ranking> cart_product(std::vector<ranking> rankings, std::tuple<int, int, bool> state)
    {
        int state_name = std::get<0>(state);
        int max_rank = std::get<1>(state);
        bool accepting = std::get<2>(state);

        std::vector<ranking> result;

        for (int i=0; i<max_rank; (accepting ? i+=2 : i++))
        {
            std::vector<ranking> new_rankings(rankings.begin(), rankings.end());
            for (auto r : new_rankings)
            {
                r[state_name] = i;
                if (r.get_max_rank() < i)
                    r.set_max_rank(i);
                result.push_back(r);
            }
        }

        return result;
    }

    void check_tight(std::vector<ranking>& rankings)
    {
        std::vector<ranking> result;
        for (auto r : rankings)
        {
            bool skip = false;
            if (r.get_max_rank() % 2 == 1)
            {
                // check box
                int max_rank = -1;
                if (r.find(BOX) != r.end())
                {
                    if (r.get_max_rank() != r[BOX])
                        skip = true;
                    max_rank = r[BOX];
                }

                if (not skip)
                {
                    std::set<int> ranks;
                    for (auto pr : r)
                    {
                        if (pr.first != BOX and pr.second == max_rank)
                        {
                            skip = true;
                            break;
                        }
                        ranks.insert(pr.second);
                    }
                    if (not skip)
                    {
                        for (int i=1; i<r.get_max_rank(); i+=2)
                        {
                            if (ranks.find(i) == ranks.end())
                            {
                                skip = true;
                                break;
                            }    
                        }
                        if (not skip)
                            result.push_back(r);
                    }
                }
            }
        }

        rankings = result;
    }
    
    std::vector<ranking> get_tight_rankings(std::vector<std::tuple<int, int, bool>> mp)
    {
        std::vector<ranking> rankings;

        // get max rank
        if (mp.size() > 0)
        {
            auto max = std::max_element(mp.begin(), mp.end(), compare_ranks);
            int max_rank = std::get<1>(*max);
            if (max_rank > 2*mp.size())
                max_rank = 2*mp.size();

            // odd rank
            if (max_rank == 0)
                return rankings;
            if (max_rank % 2 == 0)
                max_rank--;

            for (auto state : mp)
            {
                if (rankings.size() == 0)
                {
                    for (int i=0; i<std::get<1>(state); (std::get<2>(state) ? i+=2 : i++))
                    {
                        ranking r;
                        r[std::get<0>(state)] = i;
                        r.set_max_rank(i);
                        rankings.push_back(r); 
                    }
                } 
                else
                {
                    rankings = cart_product(rankings, state);
                }
            }

            check_tight(rankings);
        }

        return rankings;
    }

    bool ranking::is_bigger(ranking other)
    {
        for (auto it=this->begin(); it!=this->end(); it++)
        {
            if (it->second < other[it->first])
                return false;
        }
        return true;
    }
}