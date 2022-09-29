#pragma once

#include "cola.hpp"

namespace cola
{   
    class waiting
    {
        private:
            std::set<std::set<unsigned>> states_;
            std::map<std::set<unsigned>, std::set<std::set<unsigned>>> trans_;
            std::map<std::set<unsigned>, std::set<std::set<unsigned>>> predecessors_;

        public:
            waiting(std::set<std::set<unsigned>> states, std::map<std::set<unsigned>, std::set<std::set<unsigned>>> trans) : states_(states), trans_(trans)
            {
                std::map<std::set<unsigned>, std::set<std::set<unsigned>>> predecessors;

                for (auto state : states)
                    predecessors.insert({state, std::set<std::set<unsigned>>()});

                for (auto src : states)
                {
                    for (auto dst : trans[src])
                    {
                        predecessors[dst].insert(src);
                    }
                }

                predecessors_ = predecessors;
            }

            std::set<std::set<unsigned>> get_states()
            {
                return states_;
            }

            std::map<std::set<unsigned>, std::set<std::set<unsigned>>> get_predecessors()
            {
                return predecessors_;
            }
    };
}