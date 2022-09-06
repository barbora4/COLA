#include "complement_class.hpp"

namespace cola
{
    std::set<unsigned> 
    complement_class::get_all_successors(std::vector<unsigned> current_states, bdd symbol)
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
}
