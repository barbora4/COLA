#pragma once

#include <set>
#include <map>
#include <vector>
#include <string>

#include "cola.hpp"

namespace cola
{
    class ranking : public std::map<int, int>
    {   
        private:
            unsigned max_rank = 0;
        
        public:
            ranking() : std::map<int, int>(){};
            std::string get_name();
            unsigned get_max_rank(){return max_rank;};
            void set_max_rank(unsigned max_rank){this->max_rank = max_rank;};
            void check_tight(std::vector<ranking> rankings);
    };

    static bool compare_ranks(std::tuple<int, int, bool> first, std::tuple<int, int, bool> second);
    std::vector<ranking> get_tight_rankings(std::vector<std::tuple<int, int, bool>> mp);
    std::vector<ranking> cart_product(std::vector<ranking> rankings, std::tuple<int, int, bool> state);
}