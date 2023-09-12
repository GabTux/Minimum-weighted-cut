#ifndef SEMESTRALKA_SOLUTION_HPP
#define SEMESTRALKA_SOLUTION_HPP


#include <bitset>
#include <limits>


class Solution {
public:
    static const int MAX_VERTEX = 140;

    Solution() = default;
    Solution(int min, std::bitset<MAX_VERTEX> setX);
    bool operator<(const Solution& other) const;
    friend std::ostream& operator<<(std::ostream& os, const Solution& solution);

    int m_minCut = std::numeric_limits<int>::max();
    std::bitset<MAX_VERTEX> m_result;
};


#endif //SEMESTRALKA_SOLUTION_HPP
