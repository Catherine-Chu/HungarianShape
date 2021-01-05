//
// Created by 褚文杰 on 2020/2/24.
//

#ifndef HUNGARIANSHAPE2_AGENT_H
#define HUNGARIANSHAPE2_AGENT_H


#include<unordered_map>
#include<tuple>
#include<cstdlib>
#include<ctime>
#include <chrono>
#include <mutex>
#include <vector>
#include <random>
#include <iostream>
#include <string>
#include <iterator>
#include <limits>
#include <list>
#include <type_traits>
#include <fstream>
#include <algorithm>

const int width = 80;
const int height = 80;

class Agent {
public:
    int index;
    int pos_x;
    int pos_y;
    int sense_r;
    int tar[2];
    bool is_running;
    bool not_end;
    bool not_move;
    std::unique_lock<std::mutex> guard;

    Agent();
    Agent(int id, int px, int py, int r = 2);
    void set_config(int id, int px, int py, int r = 2);
    void set_tar(int tar_x, int tar_y);
    void parallel_running(bool block_flag);
};


#endif //HUNGARIANSHAPE2_AGENT_H
