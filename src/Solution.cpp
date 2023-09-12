#include <iostream>
#include "Solution.hpp"

Solution::Solution(int min, std::bitset<MAX_VERTEX> setX) :
m_minCut(min), m_result(setX)
{

}

bool Solution::operator<(const Solution &other) const {
    return m_minCut < other.m_minCut;
}

std::ostream& operator<<(std::ostream& os, const Solution& solution)
{
    os << "minCut: " << solution.m_minCut << std::endl;
    os << "setX: " << solution.m_result << std::endl;
    return os;
}