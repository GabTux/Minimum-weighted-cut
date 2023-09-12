#ifndef SEMESTRALKA_STATE_HPP
#define SEMESTRALKA_STATE_HPP

#include <bitset>
#include "Solution.hpp"

class State {
public:
    std::bitset<Solution::MAX_VERTEX> setX;
    int currentIndex{};
    int xCount{};
    int cut{};
    Solution currentBest{};

    bool operator < (const State& other) const;

    friend std::ostream& operator<<(std::ostream& os, const State& state);
};


#endif //SEMESTRALKA_STATE_HPP
