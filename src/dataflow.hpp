#pragma once

#include "cola.hpp"

namespace cola
{
    typedef std::vector<std::set<unsigned>> wait_state;
    
    class waiting
    {
        private:
            std::set<wait_state> states_;
            std::map<wait_state, std::set<wait_state>> trans_;

        public:
            waiting(std::set<wait_state> states, std::map<wait_state, std::set<wait_state>> trans) : states_(states), trans_(trans) {}
    };
}