////
//// Created by 褚文杰 on 2020/2/24.
////
//
#include "Agent.h"
using namespace std;
extern mutex locks[width][height];
extern vector<vector<int>> agent_poses;
extern vector<vector<int>> agent_tars;

Agent::Agent() {
    index = -1;
    tar[0]= -1;
    tar[1]= -1;
    is_running = true;
    not_end = false;
    not_move = true;
}

Agent::Agent(int id, int px, int py, int r) {
    index = id;
    pos_x = px;
    pos_y = py;
    if (guard.owns_lock()) {
        guard.unlock();
        guard.release();
    }
    guard = unique_lock<mutex>(locks[pos_x][pos_y], defer_lock);
    guard.lock();
    sense_r = r;
    tar[0]= -1;
    tar[1]= -1;
    is_running = true;
    not_end = false;
    not_move = true;
}

void Agent::set_config(int id, int px, int py, int r) {
    index = id;
    pos_x = px;
    pos_y = py;
    if (guard.owns_lock()) {
        guard.unlock();
    }
    guard = unique_lock<mutex>(locks[pos_x][pos_y], defer_lock);
    guard.lock();
    sense_r = r;
}

void Agent::set_tar(int tar_x, int tar_y){
    tar[0]=tar_x;
    tar[1]=tar_y;
}

void Agent::parallel_running(bool block_flag) {

    int x_direct = 0;
    int y_direct = 0;
    if(tar[0]-pos_x>0){
        x_direct = 1;
    }else if (tar[0]-pos_x<0){
        x_direct = -1;
    }
    if(tar[1]-pos_y>0){
        y_direct = 1;
    }else if (tar[1]-pos_y<0){
        y_direct = -1;
    }

    int lean_step = min(abs(tar[1]-pos_y), abs(tar[0]-pos_x));
    int y_step = 0;
    int x_step = 0;
    if(abs(tar[1]-pos_y)>abs(tar[0]-pos_x)){
        y_step = abs(tar[1]-pos_y)-abs(tar[0]-pos_x);
    }else if(abs(tar[1]-pos_y)<abs(tar[0]-pos_x)){
        x_step = abs(tar[0]-pos_x)-abs(tar[1]-pos_y);
    }

    std::vector<std::vector<int>> optional_actions;//all actions to approach to the goal.
    std::vector<std::vector<int>> rand_less;//all actions to approach to the goal.
    if(lean_step>0){
        if(pos_x+x_direct>=0 && pos_x+x_direct<width && pos_y+y_direct>=0 && pos_y+y_direct<height) {
            optional_actions.push_back({pos_x + x_direct, pos_y + y_direct});
        }
        if(pos_x-x_direct>=0 && pos_x-x_direct<width)
            rand_less.push_back({pos_x-x_direct,pos_y+y_direct});
        if(pos_y-y_direct>=0 && pos_y-y_direct<height)
            rand_less.push_back({pos_x+x_direct,pos_y-y_direct});
        if(pos_x-x_direct>=0 && pos_x-x_direct<width && pos_y-y_direct>=0 && pos_y-y_direct<height)
            rand_less.push_back({pos_x-x_direct,pos_y-y_direct});
    }else{
        if(pos_x+1<width && pos_y+1 < height)
            rand_less.push_back({pos_x+1,pos_y+1});
        if(pos_x-1>=0 && pos_y+1 < height)
            rand_less.push_back({pos_x-1,pos_y+1});
        if(pos_x+1<width && pos_y-1 >= 0)
            rand_less.push_back({pos_x+1,pos_y-1});
        if(pos_x-1>=0 && pos_y-1 >=0)
            rand_less.push_back({pos_x-1,pos_y-1});
    }

    if(x_step>0){
        if(pos_x+x_direct>=0 && pos_x+x_direct<width) {
            optional_actions.push_back({pos_x + x_direct, pos_y});
        }
        if(pos_x-x_direct>=0 && pos_x-x_direct<width)
            rand_less.push_back({pos_x-x_direct,pos_y});
    }else{
        if(pos_x+1<width)
            rand_less.push_back({pos_x+1,pos_y});
        if(pos_x-1>=0)
            rand_less.push_back({pos_x-1,pos_y});
    }
    if(y_step>0){
        if(pos_y+y_direct>=0 && pos_y+y_direct<height) {
            optional_actions.push_back({pos_x, pos_y + y_direct});
        }
        if(pos_y-y_direct>=0 && pos_y-y_direct<height)
            rand_less.push_back({pos_x,pos_y-y_direct});
    }else{
        if(pos_y+1<height)
            rand_less.push_back({pos_x,pos_y+1});
        if(pos_y-1>=0)
            rand_less.push_back({pos_x,pos_y-1});
    }
    optional_actions.push_back({pos_x,pos_y});

    // conflict avoidance
    int next_x = pos_x;
    int next_y = pos_y;

    unique_lock<mutex> temp_guard;

    for (int i = 0; i <= optional_actions.size() - 1; i++) {
        temp_guard = unique_lock<mutex>(locks[optional_actions[i][0]][optional_actions[i][1]], defer_lock);
        temp_guard.try_lock();
        if (temp_guard.owns_lock()) {
            next_x = optional_actions[i][0];
            next_y = optional_actions[i][1];
            agent_poses[pos_x][pos_y] = -1;
            agent_poses[next_x][next_y] = index;
            guard.unlock();
            guard.swap(temp_guard);
            break;
        } else {
            continue;
        }
    }


    if (next_x == pos_x && next_y == pos_y) {
        if(is_running) {
//        cout<<"agent "<<index<<", position: "<<pos_x<<" "<<pos_y<<endl;
//        cout<<"opt acs: "<<optional_actions.size()<<endl;
//        cout<<"goal is: "<<tar[0]<<" "<<tar[1]<<endl;
            not_end = false;
            for (int i = 0; i < optional_actions.size() - 1; i++) {
//            cout<<optional_actions[i][0]<<" "<<optional_actions[i][1]<<endl;
                int oc_id = agent_poses[optional_actions[i][0]][optional_actions[i][1]];
                if (oc_id > -1 && agent_tars[oc_id][0] == optional_actions[i][0] &&
                    agent_tars[oc_id][1] == optional_actions[i][1]) {
                    not_end = not_end || false;
                } else {
                    not_end = not_end || true;
                }
            }
//        cout<<"not end: "<<not_end<<endl;
        }else{
            not_end = false;
        }

        if((!is_running) && block_flag) {
            srand((unsigned) time(NULL));
            double p = (float) rand() / RAND_MAX;
            double prob = 1;
            if (p < prob) {
                if (rand_less.size() > 0) {
                    default_random_engine generator{random_device{}()};
                    shuffle(rand_less.begin(), rand_less.end(), generator);
                    for (int t = 0; t <= (rand_less.size() - 1); t++) {
                        temp_guard = unique_lock<mutex>(locks[rand_less[t][0]][rand_less[t][1]], defer_lock);
                        temp_guard.try_lock();
                        if (temp_guard.owns_lock()) {
                            next_x = rand_less[t][0];
                            next_y = rand_less[t][1];
                            not_end = true;
                            agent_poses[pos_x][pos_y] = -1;
                            agent_poses[next_x][next_y] = index;
                            guard.unlock();
                            guard.swap(temp_guard);
                            break;
                        } else {
                            continue;
                        }
                    }
                }
            }
        }
    }else{
        not_end = true;
    }
    //transmit to new positions and record new state
    if(next_x == pos_x && next_y == pos_y)
        not_move = true;
    else
        not_move = false;
    pos_x = next_x;
    pos_y = next_y;

    if (is_running && tar[0] == pos_x && tar[1] == pos_y) {
        is_running = false;
    }else if((!is_running) && (tar[0]!=pos_x || tar[1]!= pos_y)){
        is_running = true;
    }
}
