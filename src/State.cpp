#include "State.hpp"
#include <iostream>

std::ostream& operator<<(std::ostream& os, const State& state)
{
    os << "setX: " << state.setX << std::endl;
    os << "currentIndex: " << state.currentIndex << std::endl;
    os << "xCount: " << state.xCount << std::endl;
    os << "cut: " << state.cut << std::endl;
    return os;
}

bool State::operator < (const State& other) const {
    return cut < other.cut;
}