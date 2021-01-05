/* This is an implementation of the Hungarian algorithm in C++
 * The Hungarian algorithm, also know as Munkres or Kuhn-Munkres
 * algorithm is usefull for solving the assignment problem.
 *
 * Assignment problem: Let C be an n x n matrix
 * representing the costs of each of n workers to perform any of n jobs.
 * The assignment problem is to assign jobs to workers so as to
 * minimize the total cost. Since each worker can perform only one job and
 * each job can be assigned to only one worker the assignments constitute
 * an independent set of the matrix C.
 *
 * It is a port heavily based on http://csclab.murraystate.edu/~bob.pilgrim/445/munkres.html
 *
 * This version is written by Fernando B. Giannasi */


#include <cmath>
#include <iostream>
#include <iterator>
#include <limits>
#include <list>
#include <string>
#include <type_traits>
#include <vector>
#include <random>
#include <thread>
#include <dirent.h>
#include <fstream>
#include <queue>
#include "Agent.h"
//#include <dlib/optimization/max_cost_assignment.h>

using namespace std;
//using namespace dlib;

//const string RUN_ENV = "MAC";
//const int TEMP_THREAD_NUM = 7;
#define TEMP_THREAD_NUM 64
const string RUN_ENV = "LINUX";
#define inf 0x7f7f7f7f
vector<vector<vector<int>>> agent_init_pos_in_shapes; //7*agent_num*2
//vector<Agent> swarm;
//const int width = 80;
//const int height = 80;
mutex locks[width][height];
vector<vector<int>> agent_poses(width, vector<int>(height, -1));
vector<vector<int>> agent_tars;
vector<Agent> swarm;
enum Agent_num_mode {
    less_a = -1, equal_a = 0, more_a = 1
};
Agent_num_mode a_mode = more_a;

struct node {
    int v, w;

    node() {}

    node(int V, int W) {
        v = V;
        w = W;
    }

    bool operator<(const node &rhs) const {
        return rhs.w > w;
    }
};


const int N_NUM = width * height; //节点总数
vector<node> G[N_NUM]; //邻接矩阵,G[i]为一个node列表,表示出发点为i,到每个node.v的距离为node.w


int dijkstra(int s, int t) {
    int dis[N_NUM];
    for (int i = 1; i <= N_NUM; i++)
        dis[i] = inf;
    dis[s] = 0;
    priority_queue<node> q;
    q.push(node(s, 0));
    while (!q.empty()) {
        node x = q.top();
        q.pop();
        for (int i = 0; i < G[x.v].size(); i++) {
            node y = G[x.v][i];
            if (dis[y.v] > x.w + y.w) {
                dis[y.v] = x.w + y.w;
                q.push(node(y.v, dis[y.v]));
            }
        }
    }
    return dis[t];
}

void parallel_swarms(int t_id, int s_id, int e_id, bool block_flag) {
    for (int i = s_id; i <= e_id; i++) {
        swarm[i].parallel_running(block_flag);
    }
}

double getMold(const vector<vector<int>> &vec) {   //求向量的模长
    int n = vec.size();
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        int m = vec[i].size();
        for (int j = 0; j < m; ++j) {
            sum += vec[i][j] * vec[i][j];
        }
    }
    return sqrt(sum);
}

double getSimilarity(const vector<vector<int>> &lhs, const vector<vector<int>> &rhs) {
    int n = lhs.size();
    if (n == rhs.size()) {
        double tmp = 0.0;  //内积
        for (int i = 0; i < n; ++i) {
            int m = lhs[i].size();
            if (m == rhs[i].size()) {
                for (int j = 0; j < m; ++j) {
                    tmp += lhs[i][j] * rhs[i][j];
                }
            } else {
                return -1;
            }
        }
        return tmp / (getMold(lhs) * getMold(rhs));
    } else {
        return -1;
    }
}


bool initialize_no_seed_agent_positions(vector<vector<int>> &grids, int f, int min_i, int min_j, int shape_agent_num) {
    int cnt = 4;
    bool is_enough = false;
    if (min_j - 3 < 4 or min_i < 4 or min_i >= width - 4) {
        cout << "Too narrow space for edge following." << endl;
        return is_enough;
    }

    for (int j = min_j - 3; j >= 4; j--) {
//        swarm[cnt].set_config(min_i,j,2*sqrt(2));
        agent_init_pos_in_shapes[f].push_back({min_i, j});
        cnt++;
        if (cnt == shape_agent_num) {
            is_enough = true;
            break;
        }
    }
//    int max_height = int(0.6 * height);
    int max_height = height - 3;
    vector<vector<int>> incre_options;
    if (not is_enough) {
        for (int i = min_i + 1; i < width - 4; i++) {
            for (int j = 4; j <= max_height; j++) {
                if (j > min_j - 3) {
                    bool is_valid = true;
                    for (int px = i - 2; px <= i + 2 && px < width - 1; px++) {
                        for (int py = j - 2; py <= j + 2 && py < height - 1; py++) {
//                            cout<<px<<" "<<py<<endl;
                            if (grids[px][py] == 0) {
                                is_valid = is_valid && true;
                            } else {
                                is_valid = is_valid && false;
                                break;
                            }
                        }
                        if (not is_valid) {
                            break;
                        }
                    }
                    if (is_valid) {
//                        swarm[cnt].set_config(i,j,2*sqrt(2));
                        agent_init_pos_in_shapes[f].push_back({i, j});
                        cnt++;
                    } else {
                        break;
                    }
                } else {
//                    swarm[cnt].set_config(i,j,2*sqrt(2));
                    agent_init_pos_in_shapes[f].push_back({i, j});
                    cnt++;
                }
                if (cnt == shape_agent_num) {
                    is_enough = true;
                    break;
                }
                if (j == max_height) {
                    incre_options.push_back({i, j});
                }
            }
            if (is_enough) {
                break;
            }
        }
    }
    if (not is_enough) {
        for (int i = min_i - 1; i >= 4; i--) {
            for (int j = 4; j <= max_height; j++) {
                if (j > min_j - 3) {
                    bool is_valid = true;
                    for (int px = i - 2; px <= i + 3 && px < width - 1; px++) {
                        for (int py = j - 2; py <= j + 2 && py < height - 1; py++) {
//                            cout<<px<<" "<<py<<endl;
                            if (grids[px][py] == 0) {
                                is_valid = is_valid && true;
                            } else {
                                is_valid = is_valid && false;
                                break;
                            }
                        }
                        if (not is_valid) {
                            break;
                        }
                    }
                    if (is_valid) {
//                        swarm[cnt].set_config(i,j,2*sqrt(2));
                        agent_init_pos_in_shapes[f].push_back({i, j});
                        cnt++;
                    } else {
                        break;
                    }
                } else {
//                    swarm[cnt].set_config(i,j,2*sqrt(2));
                    agent_init_pos_in_shapes[f].push_back({i, j});
                    cnt++;
                }
                if (cnt == shape_agent_num) {
                    is_enough = true;
                    break;
                }
                if (j == max_height) {
                    incre_options.push_back({i, j});
                }
            }
            if (is_enough) {
                break;
            }
        }
    }
    if (not is_enough) {
        for (int o = 0; o < incre_options.size(); o++) {
            int i = incre_options[o][0];
            for (int j = incre_options[o][1] + 1; j < height - 2; j++) {
                if (j > min_j - 3) {
                    bool is_valid = true;
                    for (int px = i - 2; px <= i + 2 && px < width - 1; px++) {
                        for (int py = j - 2; py <= j + 2 && py < height - 1; py++) {
//                            cout<<px<<" "<<py<<endl;
                            if (grids[px][py] == 0) {
                                is_valid = is_valid && true;
                            } else {
                                is_valid = is_valid && false;
                                break;
                            }
                        }
                        if (not is_valid) {
                            break;
                        }
                    }
                    if (is_valid) {
//                        swarm[cnt].set_config(i,j,2*sqrt(2));
                        agent_init_pos_in_shapes[f].push_back({i, j});
                        cnt++;
                    } else {
                        break;
                    }
                } else {
//                    swarm[cnt].set_config(i,j,2*sqrt(2));
                    agent_init_pos_in_shapes[f].push_back({i, j});
                    cnt++;
                }
                if (cnt == shape_agent_num) {
                    is_enough = true;
                    break;
                }
            }
            if (is_enough) {
                break;
            }
        }
    }
    if (not is_enough) {
        cout << "There isn't enough space for initialization." << endl;

    }
    cout << shape_agent_num << " " << cnt << endl;
    return is_enough;
}

/*specific initialization*/
//int main() {
//    // Let's imagine you need to assign N people to N jobs.  Additionally, each person will make
//    // your company a certain amount of money at each job, but each person has different skills
//    // so they are better at some jobs and worse at others.  You would like to find the best way
//    // to assign people to these jobs.  In particular, you would like to maximize the amount of
//    // money the group makes as a whole.  This is an example of an assignment problem and is
//    // what is solved by the max_cost_assignment() routine.
//    //
//    // So in this example, let's imagine we have 3 people and 3 jobs.  We represent the amount of
//    // money each person will produce at each job with a cost matrix.  Each row corresponds to a
//    // person and each column corresponds to a job.  So for example, below we are saying that
//    // person 0 will make $1 at job 0, $2 at job 1, and $6 at job 2.
//    dlib::matrix<int> cost(3, 4);
//    cost = 1, 2, 6, 4,
//            5, 3, 6, 8,
//            4, 5, 0, 0;
//
//    // To find out the best assignment of people to jobs we just need to call this function.
//    std::vector<long> assignment = max_cost_assignment(cost);
//
//    // This prints optimal assignments:  [2, 0, 1] which indicates that we should assign
//    // the person from the first row of the cost matrix to job 2, the middle row person to
//    // job 0, and the bottom row person to job 1.
//    for (unsigned int i = 0; i < assignment.size(); i++)
//        cout << assignment[i] << std::endl;
//
//    // This prints optimal cost:  16.0
//    // which is correct since our optimal assignment is 6+5+5.
//    cout << "optimal cost: " << assignment_cost(cost, assignment) << endl;
//}

//int main() {
//    string dictionary, out_dict;
//    if (RUN_ENV == "WIN") {
//        dictionary =
//                "D:\\projects\\CLionProjects\\InitSettingGenerator\\tmp_edge_display\\" + to_string(width) + '_' +
//                to_string(height);
//        out_dict = "D:\\projects\\CLionProjects\\HungarianShape2\\exp";
//    } else if (RUN_ENV == "MAC") {
//        dictionary =
//                "/Users/chuwenjie/CLionProjects/InitSettingGenerator/tmp_edge_display/" + to_string(width) + '*' +
//                to_string(height);
//        out_dict = "/Users/chuwenjie/CLionProjects/HungarianShape2/exp";
//    } else {
//        dictionary = "./exp/" + to_string(width) + '*' + to_string(height);
//        out_dict = "./exp";
//    }
//
////    string work_root = argv[1];
////    string dictionary, out_dict;
////    if (RUN_ENV == "WIN") {
////        dictionary = work_root +
////                     "\\inputs\\"+ argv[2]+"\\"+ to_string(width) + '_' + to_string(height);
////        out_dict = work_root + "\\outputs\\" + argv[2] +"\\"+to_string(width) + '_' + to_string(height);
////    } else {
////        dictionary = work_root + "/inputs/" + argv[2]+ "/" + to_string(width) + '_' + to_string(height);
////        out_dict = work_root + "/outputs/" + argv[2]+ "/" + to_string(width) + '_' + to_string(height);
////    }
//
//    DIR *dir;
//    struct dirent *ptr;
//    vector<string> name_posts;
//    vector<int> shape_nums;
//    vector<string> filelist;
//    vector<int> a_num_s;
//    const char *p = dictionary.c_str();
//    if ((dir = opendir(p)) == NULL) {
//        perror("Open dir error...");
//        exit(1);
//    }
//    while ((ptr = readdir(dir)) != NULL) {
//        if (string(ptr->d_name).compare(0, 4, "grid") == 0)    //file
//        {
//            string temp;
//            if (RUN_ENV == "WIN") {
//                temp = dictionary + '\\' + ptr->d_name;
//            } else {
//                temp = dictionary + '/' + ptr->d_name;
//            }
//            filelist.push_back(temp);
//            string _post = ptr->d_name;
//            int shape_num = _post[5] - '0';
//            shape_nums.push_back(shape_num);
//            int a_num = atoi(_post.substr(7, _post.size() - 11).c_str());
//            a_num_s.push_back(a_num);
//            _post = _post.substr(4, _post.size() - 4);
//            name_posts.push_back(_post);
//        }
//    }
//    closedir(dir);
//
//    for (int f = 0; f < filelist.size(); f++) {
//        vector<vector<int>> grids;
//        for (int x = 0; x < width; x++) {
//            grids.emplace_back();
//            for (int y = 0; y < height; y++) {
//                grids[x].push_back(0);
//            }
//        }
//
//        int shape_agent_num;
//        for (int k = 0; k < a_num_s[f]; k++) {
//            swarm.push_back(Agent());
//            agent_tars.push_back({-1,-1});
//        }
//        srand(time(NULL));
//        int rand_more = 4;
//        if (a_mode == more_a) {
//            //如果是多,则为N+3~N+0.05N+3
//            for (int i = 0; i < rand_more + 3; i++) {
//                swarm.push_back(Agent());
//                agent_tars.push_back({-1,-1});
//            }
//        } else if (a_mode == equal_a) {
//            //如果是正好,由于seed robot中只有v3在目标内,所以至少为N+3
//            for (int i = 0; i < 3; i++) {
//                swarm.push_back(Agent());
//                agent_tars.push_back({-1,-1});
//            }
//        } else {
//            rand_more = rand_more - 3;
//            //如果是少,则为N+3-0.05N~N+3
//            if (rand_more > 0) {
//                for (int i = 0; i < rand_more; i++) {
//                    swarm.erase(swarm.end());
//                    agent_tars.erase(agent_tars.end());
//                }
//            } else if (rand_more < 0) {
//                for (int i = 0; i < -rand_more; i++) {
//                    swarm.push_back(Agent());
//                    agent_tars.push_back({-1,-1});
//                }
//            }
//        }
//        shape_agent_num = swarm.size();
//
//        //read the grid environment
//        fstream infile;
//        infile.open(filelist[f], ios::in);
//        if (!infile) {
//            cout << "open failed" << endl;
//            exit(1);
//        }
//        int i = 0;
//        int j = height - 1;
//        vector<vector<int>> goals;
//        int min_i = width;
//        int min_j = height;
//        while (!infile.eof() and j >= 0) {
//            i = 0;
//            while (!infile.eof() and i < width) {
//                infile >> grids[i][j];
//                if (grids[i][j] == 1) {
//                    goals.push_back({i, j});
//                    if (j <= min_j) {
//                        if (j < min_j) {
//                            min_j = j;
//                            min_i = i;
//                        } else {
//                            if (i < min_i) {
//                                min_i = i;
//                            }
//                        }
//                    }
//                }
//                i++;
//            }
//            j--;
//        }
//        infile.close();
//
//        vector<vector<int>> agent_init_pos_in_a_shape;
//        agent_init_pos_in_shapes.push_back(agent_init_pos_in_a_shape);
//        //initialize seed robot v0-v3, noting that only v3 is within the target shape,
//        //so the total number of agent for n-n shape is at least n+3 agent.
//        agent_init_pos_in_shapes[f].push_back({min_i, min_j - 1});
//        agent_init_pos_in_shapes[f].push_back({min_i + 1, min_j - 1});
//        agent_init_pos_in_shapes[f].push_back({min_i, min_j - 2});
//        agent_init_pos_in_shapes[f].push_back({min_i, min_j});
//        //初始化剩余agent位置
//        bool is_enough = initialize_no_seed_agent_positions(grids, f, min_i, min_j, shape_agent_num);
//        if (not is_enough) {
//            cout << "Experiment failed!" << endl;
//            continue;
//        }
//
//        int exp_num = 20;
//        const int THREAD_NUM = TEMP_THREAD_NUM > shape_agent_num ? shape_agent_num : TEMP_THREAD_NUM;
//
//        double exp_avg_iter = 0;
//        double exp_avg_iter_t = 0;
//        int valid_exp = 0;
//        int minor_valid_exp = 0;
//        double minor_exp_avg_iter = 0;
//        double minor_exp_avg_iter_t = 0;
//        double exp_avg_similarity = 0;
//        string out_arg;
//        if (RUN_ENV == "WIN") {
//            out_arg = out_dict + "\\args_" + to_string(width) + "_" +
//                      to_string(height) + name_posts[f];
//        } else {
//            out_arg = out_dict + "/args_" + to_string(width) + "_" +
//                      to_string(height) + name_posts[f];
//        }
//        ofstream outarg(out_arg, ios::app);
//
//        for (int e = 0; e < exp_num; e++) {
//            //在每轮实验开始初始化所有agent的位置
//            swarm[0].set_config(0,agent_init_pos_in_shapes[f][0][0], agent_init_pos_in_shapes[f][0][1], 2);
//            swarm[1].set_config(1,agent_init_pos_in_shapes[f][1][0], agent_init_pos_in_shapes[f][1][1], 2);
//            swarm[2].set_config(2, agent_init_pos_in_shapes[f][2][0], agent_init_pos_in_shapes[f][2][1], 2);
//            swarm[3].set_config(3, agent_init_pos_in_shapes[f][3][0], agent_init_pos_in_shapes[f][3][1], 2);
//            for (int s = 4; s < agent_init_pos_in_shapes[f].size(); s++) {
//                swarm[s].set_config(s, agent_init_pos_in_shapes[f][s][0], agent_init_pos_in_shapes[f][s][1], 2);
//            }
//
//            cout<<"Begin allocation"<<endl;
//            vector<vector<double>> matrix;
//            for (int m = 0; m < swarm.size(); m++) {
//                vector<double> task_costs;
//                for(int t = 0; t < goals.size(); t++) {
//                    double cost = sqrt(pow((swarm[m].pos_x-goals[t][0]),2)+pow((swarm[m].pos_y-goals[t][1]),2));
//                    task_costs.push_back(cost);
//                }
//                matrix.push_back(task_costs);
//            }
//            vector<vector<int>> assignment;
//            hungarian(matrix, assignment);
//            for(int m = 0; m < swarm.size(); m++){
//                for(int t = 0; t < goals.size(); t++){
//                    if(assignment[m][t] == 0)
//                        continue;
//                    else{
//                        swarm[m].set_tar(goals[t][0], goals[t][1]);
//                        agent_tars[m][0] = goals[t][0];
//                        agent_tars[m][1] = goals[t][1];
//                        break;
//                    }
//                }
//            }
//            cout<<"End allocation"<<endl;
//
//            //record the initialization
//            string out_name;
//            if (RUN_ENV == "WIN") {
//                out_name = out_dict + "\\poses_" + to_string(e) + "_" +
//                           to_string(width) + "_" + to_string(height) +
//                           name_posts[f];
//            } else {
//                out_name = out_dict + "/poses_" + to_string(e) + "_" +
//                           to_string(width) + "_" + to_string(height) +
//                           name_posts[f];
//            }
//            ofstream outfile(out_name, ios::app);
//            outfile << "arguments: " << width << ' ' << height << ' ' << shape_agent_num << ' ' << THREAD_NUM
//                    << endl;
//            outfile << "agent positions:" << endl;
//            for (int k = 0; k < shape_agent_num; k++) {
//                if (k < shape_agent_num - 1) {
//                    outfile << swarm[k].pos_x << ',' << swarm[k].pos_y << ' ';
//                } else {
//                    outfile << swarm[k].pos_x << ',' << swarm[k].pos_y;
//                }
//                agent_poses[swarm[k].pos_x][swarm[k].pos_y] = swarm[k].index;
//            }
//            outfile << endl;
//
//            vector<thread> threads;
//            int terminal = 0;
//            int minor_terminal = 0;
//            clock_t startT, endT;
//            vector<double> dec_times;
//
//            int ac_tar = shape_agent_num;
//            int out_agent = shape_agent_num-1;
//            vector<int> ac_tar_decay;
//            ac_tar_decay.push_back(ac_tar);
//            vector<int> out_a_decay;
//            out_a_decay.push_back(out_agent);
//            int re_assign_cnt = 0;
//            while (ac_tar > rand_more+3 && out_agent>rand_more+3 && terminal < 1000 and minor_terminal < 500) {
//                cout << ac_tar <<" "<<out_agent<< endl;
//                terminal += 1;
//                startT = clock();
//                int left = int(shape_agent_num % THREAD_NUM);
//                int alloc = 0;
//                int s_ids[THREAD_NUM];
//                int e_ids[THREAD_NUM];
//                srand(time(NULL));
//                for (int k = 0; k < THREAD_NUM; k++) {
//                    s_ids[k] = -1, e_ids[k] = -1;
//                    if (left > alloc) {
//                        s_ids[k] = k * (int(shape_agent_num / THREAD_NUM) + 1);
//                        e_ids[k] = s_ids[k] + int(shape_agent_num / THREAD_NUM);
//                        alloc++;
//                    } else {
//                        s_ids[k] = alloc * (int(shape_agent_num / THREAD_NUM) + 1) +
//                                   (k - alloc) * int(shape_agent_num / THREAD_NUM);
//                        e_ids[k] = s_ids[k] + int(shape_agent_num / THREAD_NUM) - 1;
//                    }
//                }
//                bool block_flag = true;
//                if (terminal > 1) {
//                    int cnt = 0;
//                    for (int a = 0; a < swarm.size(); a++) {
//                        if ((!swarm[a].is_running) || swarm[a].not_move) {
//                            cnt += 1;
//                            block_flag = block_flag && true;
//                        } else {
//                            block_flag = block_flag && false;
//                        }
//                    }
//                    cout << "blk: " << cnt << endl;
//                } else {
//                    block_flag = false;
//                }
//
//                for (int k = 0; k < THREAD_NUM; k++) {
//                    threads.emplace_back(parallel_swarms, k, s_ids[k], e_ids[k], block_flag);
//                }
//                // 等待其他线程join
//                for (int k = 0; k < THREAD_NUM; k++) {
//                    threads[k].join();
//                }
//                threads.clear();
//
//                if (block_flag ) {
//                    cout << "A block happens, reassignment......" << endl;
//                    re_assign_cnt += 1;
//                    //some agents are blocked
//                    vector<vector<double>> clear_matrix(0, vector<double>(0));
//                    matrix.swap(clear_matrix);
//                    for (int a = 0; a < swarm.size(); a++) {
//                        vector<double> task_costs;
//                        for (int g = 0; g < goals.size(); g++) {
//                            int dis = sqrt(pow((swarm[a].pos_x - goals[g][0]),2)+pow((swarm[a].pos_y - goals[g][1]),2));
//                            task_costs.push_back(dis);
//                        }
//                        matrix.push_back(task_costs);
//                    }
//                    vector<vector<int>> clear_assign(0, vector<int>(0));
//                    assignment.swap(clear_assign);
//                    hungarian(matrix, assignment);
//                    for (int m = 0; m < swarm.size(); m++) {
//                        for (int t = 0; t < goals.size(); t++) {
//                            if (assignment[m][t] == 0)
//                                continue;
//                            else {
//                                if(swarm[m].tar[0] != goals[t][0] || swarm[m].tar[1] != goals[t][1]) {
//                                    swarm[m].is_running = true;
//                                }else{
//                                    swarm[m].is_running = false;
//                                }
//                                swarm[m].set_tar(goals[t][0], goals[t][1]);
//                                agent_tars[m][0] = goals[t][0];
//                                agent_tars[m][1] = goals[t][1];
//                                break;
//                            }
//                        }
//                    }
//                }
//
//
//                ac_tar = 0;
//                out_agent = 0;
//                for(int a = 0; a < swarm.size(); a++){
//                    if(swarm[a].is_running){
//                        ac_tar++;
//                    }
//                    if(grids[swarm[a].pos_x][swarm[a].pos_y]==0){
//                        out_agent++;
//                    }
//                }
//
//                endT = clock();
//                dec_times.push_back((double) (endT - startT));
//
//                //record new positions for all agents
//                outfile << "agent positions:" << endl;
//                for (int k = 0; k < shape_agent_num; k++) {
//                    if (k < shape_agent_num - 1) {
//                        outfile << swarm[k].pos_x << ',' << swarm[k].pos_y << ' ';
//                    } else {
//                        outfile << swarm[k].pos_x << ',' << swarm[k].pos_y;
//                    }
//                }
//                outfile << endl;
//                ac_tar_decay.push_back(ac_tar);
//                out_a_decay.push_back(out_agent);
//                if (out_agent <= rand_more+6) {
//                    minor_terminal += 1;
//                }
//            }
//            double avg_t = 0;
//            for (int k = 0; k < dec_times.size(); k++) {
//                avg_t += dec_times[k];
//            }
//            double avg_iteration = avg_t / (dec_times.size() * CLOCKS_PER_SEC);
//            outfile.flush();
//            outfile.close();
//
//            vector<vector<int>> formed_shape;
//            for(int g=0;g<width;g++){
//                formed_shape.push_back(vector<int>());
//                for(int h=0;h<height;h++){
//                    if(agent_poses[g][h]>=0){
//                        formed_shape[g].push_back(1);
//                    }else{
//                        formed_shape[g].push_back(0);
//                    }
//                }
//            }
//            double mse_similarity = getSimilarity(grids,formed_shape);
//            exp_avg_similarity = exp_avg_similarity + mse_similarity;
//            outarg << "Experiment " << e << ":" << endl;
//            outarg << "The average decision time for each iteration is: " << avg_iteration << "s." << endl;
//            outarg << "Main: program exiting after " << terminal << " steps, and " << minor_terminal
//                   << " steps for the last 3 positions, the similarity is " <<mse_similarity<< endl;
//            outarg << "Decay ac_tar line:";
//            for (int k = 0; k < ac_tar_decay.size(); k++) {
//                outarg << ' ' << ac_tar_decay[k];
//            }
//            ac_tar_decay.clear();
//            outarg << endl;
//            outarg << "Decay out_tar line:";
//            for (int k = 0; k < out_a_decay.size(); k++) {
//                outarg << ' ' << out_a_decay[k];
//            }
//            out_a_decay.clear();
//            outarg << endl;
//
//            cout << "End an experiment! Clearing..." << endl;
//
//            vector<vector<int>> array(width, vector<int>(height, -1));
//            agent_poses.swap(array);
//            for (int k = 0; k < swarm.size(); k++) {
//                if (swarm[k].guard.owns_lock()) {
//                    swarm[k].guard.unlock();
//                }
//            }
//
//            if (ac_tar == rand_more+3 or out_agent == rand_more+3) {
//                valid_exp += 1;
//                exp_avg_iter += terminal;
//                exp_avg_iter_t += avg_iteration;
//            }
//            if (ac_tar <= rand_more+6 or out_agent <= rand_more+6) {
//                minor_valid_exp += 1;
//                minor_exp_avg_iter = minor_exp_avg_iter + (terminal - minor_terminal);
//                minor_exp_avg_iter_t += avg_iteration;
//            }
//        }
//        exp_avg_iter = exp_avg_iter / valid_exp;
//        exp_avg_iter_t = exp_avg_iter_t / valid_exp;
//        exp_avg_similarity = exp_avg_similarity / exp_num;
//        minor_exp_avg_iter = minor_exp_avg_iter / minor_valid_exp;
//        minor_exp_avg_iter_t = minor_exp_avg_iter_t / minor_valid_exp;
//        cout << f << ": exp_avg_iter=" << exp_avg_iter << ", exp_avg_iter_time=" << exp_avg_iter_t << endl;
//        outarg << endl << endl;
//        outarg << "After " << exp_num << " experiments, for shape " << f <<", avg_similarity="<<exp_avg_similarity
//               <<", "<<valid_exp<<" experiments success, exp_avg_iter=" << exp_avg_iter
//               << ", exp_avg_iter_time=" << exp_avg_iter_t << ", minor_exp_avg_iter=" << minor_exp_avg_iter
//               << ", minor_exp_avg_iter_time=" << minor_exp_avg_iter_t << endl;
//        outarg.flush();
//        outarg.close();
//        swarm.clear();
//    }
//
//    return 0;
//}

/*random initialization*/

double MinCostMatching(const vector<vector<double>> &cost, vector<int> &Lmate, vector<int> &Rmate) {
    int n = int(cost.size());

    // construct dual feasible solution
    vector<double> u(n);
    vector<double> v(n);
    for (int i = 0; i < n; i++) {
        u[i] = cost[i][0];
        for (int j = 1; j < n; j++) u[i] = min(u[i], cost[i][j]);
    }
    for (int j = 0; j < n; j++) {
        v[j] = cost[0][j] - u[0];
        for (int i = 1; i < n; i++) v[j] = min(v[j], cost[i][j] - u[i]);
    }

    // construct primal solution satisfying complementary slackness
    Lmate = vector<int>(n, -1);
    Rmate = vector<int>(n, -1);
    int mated = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (Rmate[j] != -1) continue;
            if (fabs(cost[i][j] - u[i] - v[j]) < 1e-10) {
                Lmate[i] = j;
                Rmate[j] = i;
                mated++;
                break;
            }
        }
    }

    vector<double> dist(n);
    vector<int> dad(n);
    vector<int> seen(n);

    // repeat until primal solution is feasible
    while (mated < n) {

        // find an unmatched left node
        int s = 0;
        while (Lmate[s] != -1) s++;

        // initialize Dijkstra
        fill(dad.begin(), dad.end(), -1);
        fill(seen.begin(), seen.end(), 0);
        for (int k = 0; k < n; k++)
            dist[k] = cost[s][k] - u[s] - v[k];

        int j = 0;
        while (true) {

            // find closest
            j = -1;
            for (int k = 0; k < n; k++) {
                if (seen[k]) continue;
                if (j == -1 || dist[k] < dist[j]) j = k;
            }
            seen[j] = 1;

            // termination condition
            if (Rmate[j] == -1) break;

            // relax neighbors
            const int i = Rmate[j];
            for (int k = 0; k < n; k++) {
                if (seen[k]) continue;
                const double new_dist = dist[j] + cost[i][k] - u[i] - v[k];
                if (dist[k] > new_dist) {
                    dist[k] = new_dist;
                    dad[k] = j;
                }
            }
        }

        // update dual variables
        for (int k = 0; k < n; k++) {
            if (k == j || !seen[k]) continue;
            const int i = Rmate[k];
            v[k] += dist[k] - dist[j];
            u[i] -= dist[k] - dist[j];
        }
        u[s] += dist[j];

        // augment along path
        while (dad[j] >= 0) {
            const int d = dad[j];
            Rmate[j] = Rmate[d];
            Lmate[Rmate[j]] = j;
            j = d;
        }
        Rmate[j] = s;
        Lmate[s] = j;

        mated++;
    }

    double value = 0;
    for (int i = 0; i < n; i++)
        value += cost[i][Lmate[i]];

    return value;
}

void split(const string & s, vector<string> & tokens, const string & delimiters = " "){
    string::size_type lastPos = s.find_first_not_of(delimiters,0);
    string::size_type pos = s.find_first_of(delimiters, lastPos);
    while(string::npos != pos || string::npos != lastPos){
        tokens.emplace_back(s.substr(lastPos, pos-lastPos));
        lastPos = s.find_first_not_of(delimiters, pos);
        pos = s.find_first_of(delimiters, lastPos);
    }
}

// Random Initialization
//int main(int argc, char ** argv) {
////    string dictionary, out_dict;
////    if (RUN_ENV == "WIN") {
////        dictionary =
////                "D:\\projects\\CLionProjects\\InitSettingGenerator\\display\\" + to_string(width) + '_' +
////                to_string(height);
////        out_dict = "D:\\projects\\CLionProjects\\HungarianShape2\\exp";
////    } else if (RUN_ENV == "MAC") {
////        dictionary =
////                "/Users/chuwenjie/CLionProjects/InitSettingGenerator/display/" + to_string(width) + '*' +
////                to_string(height);
////        out_dict = "/Users/chuwenjie/CLionProjects/HungarianShape2/exp";
////    } else {
////        dictionary = "./exp/" + to_string(width) + '*' + to_string(height);
////        out_dict = "./exp";
////    }
//
////    string work_root = argv[1];
////    string dictionary, out_dict;
////    if (RUN_ENV == "WIN") {
////        dictionary = work_root +
////                     "\\inputs\\"+ argv[2]+"\\"+ to_string(width) + '_' + to_string(height);
////        out_dict = work_root + "\\outputs\\" + argv[2] +"\\"+to_string(width) + '_' + to_string(height);
////    } else {
////        dictionary = work_root + "/inputs/" + argv[2]+ "/" + to_string(width) + '_' + to_string(height);
////        out_dict = work_root + "/outputs/" + argv[2]+ "/" + to_string(width) + '_' + to_string(height);
//
//    string work_root = argv[1];
//    string dictionary, out_dict;
//    if (RUN_ENV == "WIN") {
//        dictionary = work_root +
//                     "\\inputs\\" + argv[2] + "\\" + to_string(width) + '_' + to_string(height) + "\\" + argv[4] ;
//        out_dict = work_root + "\\outputs\\" + argv[3] + "\\" + to_string(width) + '_' + to_string(height) + "\\" + argv[4];
//    } else {
//        dictionary = work_root + "/inputs/" + argv[2] + "/" + to_string(width) + '_' + to_string(height) + '/' + argv[4];
//        out_dict = work_root + "/outputs/" + argv[3] + "/" + to_string(width) + '_' + to_string(height) + '/' + argv[4];
//    }
//
//    DIR *dir;
//    struct dirent *ptr;
//    vector<string> name_posts;
//    vector<int> shape_nums;
//    vector<string> filelist;
//    vector<int> a_num_s;
//    const char *p = dictionary.c_str();
//    if ((dir = opendir(p)) == NULL) {
//        perror("Open dir error...");
//        exit(1);
//    }
//    while ((ptr = readdir(dir)) != NULL) {
//        if (string(ptr->d_name).compare(0, 4, "grid") == 0)    //file
//        {
//            string temp;
//            if (RUN_ENV == "WIN") {
//                temp = dictionary + '\\' + ptr->d_name;
//            } else {
//                temp = dictionary + '/' + ptr->d_name;
//            }
//            filelist.push_back(temp);
//            string _post = ptr->d_name;
//            //Add
//            vector<string> split_list;
//            split(_post, split_list, "_");
//            string shape_num_str = split_list[split_list.size()-2];
//            string a_num_str = split_list[split_list.size()-1];
//            a_num_str = a_num_str.substr(0,a_num_str.find("."));
//            int shape_num = std::stoi(shape_num_str);
//            int a_num = std::stoi(a_num_str);
//            //Add End
////            int shape_num = _post[5] - '0';
////            int a_num = atoi(_post.substr(7, _post.size() - 11).c_str());
//            shape_nums.push_back(shape_num);
//            a_num_s.push_back(a_num);
//            _post = _post.substr(4, _post.size() - 4);
//            name_posts.push_back(_post);
//        }
//    }
//    closedir(dir);
//
//    for (int f = 0; f < filelist.size(); f++) {
//        vector<vector<int>> grids;
//        for (int x = 0; x < width; x++) {
//            grids.emplace_back();
//            for (int y = 0; y < height; y++) {
//                grids[x].push_back(0);
//            }
//        }
//
//        const int THREAD_NUM = TEMP_THREAD_NUM > a_num_s[f] ? a_num_s[f] : TEMP_THREAD_NUM;
//        for (int k = 0; k < a_num_s[f]; k++) {
//            swarm.push_back(Agent());
//            agent_tars.push_back({-1, -1});
//        }
//
//        int exp_num = 30;
//        double exp_avg_iter = 0;
//        double exp_avg_iter_t = 0;
//        int valid_exp = 0;
//        int minor_valid_exp = 0;
//        double minor_exp_avg_iter = 0;
//        double minor_exp_avg_iter_t = 0;
//        double exp_avg_similarity = 0;
//        string out_arg;
//        if (RUN_ENV == "WIN") {
//            out_arg = out_dict + "\\args_" + to_string(width) + "_" +
//                      to_string(height) + name_posts[f];
//        } else {
//            out_arg = out_dict + "/args_" + to_string(width) + "_" +
//                      to_string(height) + name_posts[f];
//        }
//        ofstream outarg(out_arg, ios::app);
//
//        //read the grid environment
//        fstream infile;
//        infile.open(filelist[f], ios::in);
//        if (!infile) {
//            cout << "open failed" << endl;
//            exit(1);
//        }
//        int i = 0;
//        int j = height - 1;
//        vector<int> valid_l;
//        vector<vector<int>> goals;
//        while (!infile.eof() and j >= 0) {
//            i = 0;
//            while (!infile.eof() and i < width) {
//                infile >> grids[i][j];
//                if (grids[i][j] == 0) {
//                    valid_l.push_back(i * height + j);
//                } else {
//                    goals.push_back({i, j});
//                }
//                i++;
//            }
//            j--;
//        }
//        infile.close();
//        double avg_alloc_clock = 0;
//        int alloc_cnt = 0;
//        for (int e = 0; e < exp_num; e++) {
//            vector<int> temp_valid_l = valid_l;
//            //initialize the positions of agents
//            int max_size = temp_valid_l.size();
//            srand((unsigned) time(NULL));
//            for (int r = 0; r < a_num_s[f]; r++) {
//                int rand_p = rand() % max_size;
//                int px = (int) temp_valid_l[rand_p] / height;
//                int py = (int) temp_valid_l[rand_p] % height;
//                swarm[r].set_config(r, px, py, 2);
//                int tmp = temp_valid_l[max_size - 1];
//                temp_valid_l[max_size - 1] = temp_valid_l[rand_p];
//                temp_valid_l[rand_p] = tmp;
//                max_size--;
//            }
//            cout << "Begin allocation" << endl;
//            vector<vector<double>> matrix;
//            for (int m = 0; m < swarm.size(); m++) {
//                vector<double> task_costs;
//                for (int t = 0; t < goals.size(); t++) {
//                    double cost = sqrt(pow((swarm[m].pos_x - goals[t][0]), 2) + pow((swarm[m].pos_y - goals[t][1]), 2));
//                    task_costs.push_back(-cost);
//                }
//                matrix.push_back(task_costs);
//            }
//            vector<int> Lmate, Rmate;
//            clock_t _alloc = clock();
//            MinCostMatching(matrix, Lmate, Rmate);
//            _alloc = clock()-_alloc;
//            avg_alloc_clock += double(_alloc);
//            alloc_cnt += 1;
//            for (int m = 0; m < swarm.size(); m++) {
//                swarm[m].set_tar(goals[Lmate[m]][0], goals[Lmate[m]][1]);
//                agent_tars[m][0] = goals[Lmate[m]][0];
//                agent_tars[m][1] = goals[Lmate[m]][1];
//            }
//            cout << "End allocation" << endl;
//
//            //record the initialization
//            string out_name;
//            if (RUN_ENV == "WIN") {
//                out_name = out_dict + "\\poses_" + to_string(e) + "_" +
//                           to_string(width) + "_" + to_string(height) +
//                           name_posts[f];
//            } else {
//                out_name = out_dict + "/poses_" + to_string(e) + "_" +
//                           to_string(width) + "_" + to_string(height) +
//                           name_posts[f];
//            }
//            ofstream outfile(out_name, ios::app);
//            outfile << "arguments: " << width << ' ' << height << ' ' << a_num_s[f] << ' ' << THREAD_NUM
//                    << endl;
//            outfile << "agent positions:" << endl;
//            for (int k = 0; k < a_num_s[f]; k++) {
//                if (k < a_num_s[f] - 1) {
//                    outfile << swarm[k].pos_x << ',' << swarm[k].pos_y << ' ';
//                } else {
//                    outfile << swarm[k].pos_x << ',' << swarm[k].pos_y;
//                }
//                agent_poses[swarm[k].pos_x][swarm[k].pos_y] = swarm[k].index;
//            }
//            outfile << endl;
//
//            vector<thread> threads;
//            int terminal = 0;
//            int minor_terminal = 0;
//            clock_t startT, endT;
//            vector<double> dec_times;
//
//            int ac_tar = a_num_s[f];
//            int out_agent = a_num_s[f];
//            vector<int> ac_tar_decay;
//            ac_tar_decay.push_back(ac_tar);
//            vector<int> out_a_decay;
//            out_a_decay.push_back(out_agent);
//            int re_assign_cnt = 0;
//            while (ac_tar > 0 && out_agent > 0 && terminal < 1000 and minor_terminal < 500) {
//                cout << ac_tar << " " << out_agent << endl;
//                terminal += 1;
//                startT = clock();
//                int left = int(a_num_s[f] % THREAD_NUM);
//                int alloc = 0;
//                int s_ids[THREAD_NUM];
//                int e_ids[THREAD_NUM];
//                srand(time(NULL));
//                for (int k = 0; k < THREAD_NUM; k++) {
//                    s_ids[k] = -1, e_ids[k] = -1;
//                    if (left > alloc) {
//                        s_ids[k] = k * (int(a_num_s[f] / THREAD_NUM) + 1);
//                        e_ids[k] = s_ids[k] + int(a_num_s[f] / THREAD_NUM);
//                        alloc++;
//                    } else {
//                        s_ids[k] = alloc * (int(a_num_s[f] / THREAD_NUM) + 1) +
//                                   (k - alloc) * int(a_num_s[f] / THREAD_NUM);
//                        e_ids[k] = s_ids[k] + int(a_num_s[f] / THREAD_NUM) - 1;
//                    }
//                }
//                bool block_flag = true;
//                if (terminal > 1) {
//                    int cnt = 0;
//                    for (int a = 0; a < swarm.size(); a++) {
//                        if ((!swarm[a].is_running) || swarm[a].not_move) {
//                            cnt += 1;
//                            block_flag = block_flag && true;
//                        } else {
//                            block_flag = block_flag && false;
//                        }
//                    }
//                    cout << "blk: " << cnt << endl;
//                } else {
//                    block_flag = false;
//                }
//
//                for (int k = 0; k < THREAD_NUM; k++) {
//                    threads.emplace_back(parallel_swarms, k, s_ids[k], e_ids[k], block_flag);
//                }
//                // 等待其他线程join
//                for (int k = 0; k < THREAD_NUM; k++) {
//                    threads[k].join();
//                }
//                threads.clear();
//
//                if (block_flag) {
//                    cout << "A block happens, reassignment......" << endl;
//                    re_assign_cnt += 1;
//
//                    vector<vector<double>> clear_matrix(0, vector<double>(0));
//                    matrix.swap(clear_matrix);
//                    for (int a = 0; a < swarm.size(); a++) {
//                        vector<double> task_costs;
//                        for (int g = 0; g < goals.size(); g++) {
//                            int dis = sqrt(
//                                    pow((swarm[a].pos_x - goals[g][0]), 2) + pow((swarm[a].pos_y - goals[g][1]), 2));
//                            task_costs.push_back(dis);
//                        }
//                        matrix.push_back(task_costs);
//                    }
//                    Lmate.clear();
//                    Rmate.clear();
//                    _alloc = clock();
//                    MinCostMatching(matrix, Lmate, Rmate);
//                    _alloc = clock()-_alloc;
//                    avg_alloc_clock += double(_alloc);
//                    alloc_cnt += 1;
//                    for (int m = 0; m < swarm.size(); m++) {
//                        swarm[m].is_running =
//                                (swarm[m].tar[0] != goals[Lmate[m]][0] || swarm[m].tar[1] != goals[Lmate[m]][1]);
//                        swarm[m].set_tar(goals[Lmate[m]][0], goals[Lmate[m]][1]);
//                        agent_tars[m][0] = goals[Lmate[m]][0];
//                        agent_tars[m][1] = goals[Lmate[m]][1];
//                    }
//                }
//
//
//                ac_tar = 0;
//                out_agent = 0;
//                for (int a = 0; a < swarm.size(); a++) {
//                    if (swarm[a].is_running) {
//                        ac_tar++;
//                    }
//                    if (grids[swarm[a].pos_x][swarm[a].pos_y] == 0) {
//                        out_agent++;
//                    }
//                }
//
//                endT = clock();
//                dec_times.push_back((double) (endT - startT));
//
//                //record new positions for all agents
//                outfile << "agent positions:" << endl;
//                for (int k = 0; k < a_num_s[f]; k++) {
//                    if (k < a_num_s[f] - 1) {
//                        outfile << swarm[k].pos_x << ',' << swarm[k].pos_y << ' ';
//                    } else {
//                        outfile << swarm[k].pos_x << ',' << swarm[k].pos_y;
//                    }
//                }
//                outfile << endl;
//                ac_tar_decay.push_back(ac_tar);
//                out_a_decay.push_back(out_agent);
//                if (out_agent <= 3) {
//                    minor_terminal += 1;
//                }
//            }
//            double avg_t = 0;
//            for (int k = 0; k < dec_times.size(); k++) {
//                avg_t += dec_times[k];
//            }
//            double avg_iteration = avg_t / (dec_times.size() * CLOCKS_PER_SEC);
//            outfile.flush();
//            outfile.close();
//
//            vector<vector<int>> formed_shape;
//            for (int g = 0; g < width; g++) {
//                formed_shape.push_back(vector<int>());
//                for (int h = 0; h < height; h++) {
//                    if (agent_poses[g][h] >= 0) {
//                        formed_shape[g].push_back(1);
//                    } else {
//                        formed_shape[g].push_back(0);
//                    }
//                }
//            }
//            double mse_similarity = getSimilarity(grids, formed_shape);
//            exp_avg_similarity = exp_avg_similarity + mse_similarity;
//            outarg << "Experiment " << e << ":" << endl;
//            outarg << "The average decision time for each iteration is: " << avg_iteration << "s." << endl;
//            outarg << "Main: program exiting after " << terminal << " steps, and " << minor_terminal
//                   << " steps for the last 3 positions, the similarity is " << mse_similarity << endl;
//            outarg << "Decay ac_tar line:";
//            for (int k = 0; k < ac_tar_decay.size(); k++) {
//                outarg << ' ' << ac_tar_decay[k];
//            }
//            ac_tar_decay.clear();
//            outarg << endl;
//            outarg << "Decay out_tar line:";
//            for (int k = 0; k < out_a_decay.size(); k++) {
//                outarg << ' ' << out_a_decay[k];
//            }
//            out_a_decay.clear();
//            outarg << endl;
//
//            cout << "End an experiment! Clearing..." << endl;
//
//            vector<vector<int>> array(width, vector<int>(height, -1));
//            agent_poses.swap(array);
//            for (int k = 0; k < swarm.size(); k++) {
//                if (swarm[k].guard.owns_lock()) {
//                    swarm[k].guard.unlock();
//                }
//            }
//
//            if (ac_tar == 0 or out_agent == 0) {
//                valid_exp += 1;
//                exp_avg_iter += terminal;
//                exp_avg_iter_t += avg_iteration;
//            }
//            if (ac_tar <= 3 or out_agent <= 3) {
//                minor_valid_exp += 1;
//                minor_exp_avg_iter = minor_exp_avg_iter + (terminal - minor_terminal);
//                minor_exp_avg_iter_t += avg_iteration;
//            }
//        }
//        avg_alloc_clock = avg_alloc_clock / (alloc_cnt * CLOCKS_PER_SEC);
//        cout<<"The average allocation time is: "<<avg_alloc_clock<<endl;
//        exp_avg_iter = exp_avg_iter / valid_exp;
//        exp_avg_iter_t = exp_avg_iter_t / valid_exp;
//        exp_avg_similarity = exp_avg_similarity / exp_num;
//        minor_exp_avg_iter = minor_exp_avg_iter / minor_valid_exp;
//        minor_exp_avg_iter_t = minor_exp_avg_iter_t / minor_valid_exp;
//        cout << f << ": exp_avg_iter=" << exp_avg_iter << ", exp_avg_iter_time=" << exp_avg_iter_t << endl;
//        outarg << endl << endl;
//        outarg << "After " << exp_num << " experiments, for shape " << f << ", avg_similarity=" << exp_avg_similarity
//               << ", " << valid_exp << " experiments success, exp_avg_iter=" << exp_avg_iter
//               << ", exp_avg_iter_time=" << exp_avg_iter_t << ", minor_exp_avg_iter=" << minor_exp_avg_iter
//               << ", minor_exp_avg_iter_time=" << minor_exp_avg_iter_t << endl;
//        outarg.flush();
//        outarg.close();
//        swarm.clear();
//    }
//
//    return 0;
//}


/* Specific Initialization */
int main(int argc, char ** argv) {
//    string work_root = argv[1];
//    string dictionary, out_dict;
//    if (RUN_ENV == "WIN") {
//        dictionary = work_root +
//                     "\\inputs\\"+ argv[2]+"\\"+ to_string(width) + '_' + to_string(height);
//        out_dict = work_root + "\\outputs\\" + argv[2] +"\\"+to_string(width) + '_' + to_string(height);
//    } else {
//        dictionary = work_root + "/inputs/" + argv[2]+ "/" + to_string(width) + '_' + to_string(height);
//        out_dict = work_root + "/outputs/" + argv[2]+ "/" + to_string(width) + '_' + to_string(height);
//    }
    string work_root = argv[1];
    string dictionary, out_dict;
    if (RUN_ENV == "WIN") {
        dictionary = work_root +
                     "\\inputs\\" + argv[2] + "\\" + to_string(width) + '_' + to_string(height) + "\\" + argv[4] ;
        out_dict = work_root + "\\outputs\\" + argv[3] + "\\" + to_string(width) + '_' + to_string(height) + "\\" + argv[4];
    } else {
        dictionary = work_root + "/inputs/" + argv[2] + "/" + to_string(width) + '_' + to_string(height) + '/' + argv[4];
        out_dict = work_root + "/outputs/" + argv[3] + "/" + to_string(width) + '_' + to_string(height) + '/' + argv[4];
    }

    DIR *dir;
    struct dirent *ptr;
    vector<string> name_posts;
    vector<int> shape_nums;
    vector<string> filelist;
    vector<int> a_num_s;
    const char *p = dictionary.c_str();
    if ((dir = opendir(p)) == NULL) {
        perror("Open dir error...");
        exit(1);
    }
    while ((ptr = readdir(dir)) != NULL) {
        if (string(ptr->d_name).compare(0, 4, "grid") == 0)    //file
        {
            string temp;
            if (RUN_ENV == "WIN") {
                temp = dictionary + '\\' + ptr->d_name;
            } else {
                temp = dictionary + '/' + ptr->d_name;
            }
            filelist.push_back(temp);
            string _post = ptr->d_name;
            //Add
            vector<string> split_list;
            split(_post, split_list, "_");
            string shape_num_str = split_list[split_list.size()-2];
            string a_num_str = split_list[split_list.size()-1];
            a_num_str = a_num_str.substr(0,a_num_str.find("."));
            int shape_num = std::stoi(shape_num_str);
            int a_num = std::stoi(a_num_str);
            //Add End
//            int shape_num = _post[5] - '0';
//            int a_num = atoi(_post.substr(7, _post.size() - 11).c_str());
            shape_nums.push_back(shape_num);
            a_num_s.push_back(a_num);
            _post = _post.substr(4, _post.size() - 4);
            name_posts.push_back(_post);
        }
    }
    closedir(dir);

    for (int f = 0; f < filelist.size(); f++) {
        vector<vector<int>> grids;
        for (int x = 0; x < width; x++) {
            grids.emplace_back();
            for (int y = 0; y < height; y++) {
                grids[x].push_back(0);
            }
        }

        int shape_agent_num;
        for (int k = 0; k < a_num_s[f]; k++) {
            swarm.push_back(Agent());
            agent_tars.push_back({-1,-1});
        }
        srand(time(NULL));
        int rand_more = 0;
        if (a_mode == more_a) {
            //如果是多,则为N+3~N+0.05N+3
            for (int i = 0; i < rand_more + 3; i++) {
                swarm.push_back(Agent());
                agent_tars.push_back({-1,-1});
            }
        } else if (a_mode == equal_a) {
            //如果是正好,由于seed robot中只有v3在目标内,所以至少为N+3
            for (int i = 0; i < 3; i++) {
                swarm.push_back(Agent());
                agent_tars.push_back({-1,-1});
            }
        } else {
            rand_more = rand_more - 3;
            //如果是少,则为N+3-0.05N~N+3
            if (rand_more > 0) {
                for (int i = 0; i < rand_more; i++) {
                    swarm.erase(swarm.end());
                    agent_tars.erase(agent_tars.end());
                }
            } else if (rand_more < 0) {
                for (int i = 0; i < -rand_more; i++) {
                    swarm.push_back(Agent());
                    agent_tars.push_back({-1,-1});
                }
            }
        }
        shape_agent_num = swarm.size();

        //read the grid environment
        fstream infile;
        infile.open(filelist[f], ios::in);
        if (!infile) {
            cout << "open failed" << endl;
            exit(1);
        }
        int i = 0;
        int j = height - 1;
        vector<vector<int>> goals;
        int min_i = width;
        int min_j = height;
        while (!infile.eof() and j >= 0) {
            i = 0;
            while (!infile.eof() and i < width) {
                infile >> grids[i][j];
                if (grids[i][j] == 1) {
                    goals.push_back({i, j});
                    if (j <= min_j) {
                        if (j < min_j) {
                            min_j = j;
                            min_i = i;
                        } else {
                            if (i < min_i) {
                                min_i = i;
                            }
                        }
                    }
                }
                i++;
            }
            j--;
        }
        infile.close();

        vector<vector<int>> agent_init_pos_in_a_shape;
        agent_init_pos_in_shapes.push_back(agent_init_pos_in_a_shape);
        //initialize seed robot v0-v3, noting that only v3 is within the target shape,
        //so the total number of agent for n-n shape is at least n+3 agent.
        agent_init_pos_in_shapes[f].push_back({min_i, min_j - 1});
        agent_init_pos_in_shapes[f].push_back({min_i + 1, min_j - 1});
        agent_init_pos_in_shapes[f].push_back({min_i, min_j - 2});
        agent_init_pos_in_shapes[f].push_back({min_i, min_j});
        //初始化剩余agent位置
        bool is_enough = initialize_no_seed_agent_positions(grids, f, min_i, min_j, shape_agent_num);
        if (not is_enough) {
            cout << "Experiment failed!" << endl;
            continue;
        }

        grids[min_i][min_j-1] = 1;
        grids[min_i+1][min_j-1] = 1;
        grids[min_i][min_j-2] = 1;
        grids[min_i][min_j] = 1;

        int exp_num = 14;
        const int THREAD_NUM = TEMP_THREAD_NUM > shape_agent_num ? shape_agent_num : TEMP_THREAD_NUM;

        double exp_avg_iter = 0;
        double exp_avg_iter_t = 0;
        int valid_exp = 0;
        int minor_valid_exp = 0;
        double minor_exp_avg_iter = 0;
        double minor_exp_avg_iter_t = 0;
        double exp_avg_similarity = 0;
        string out_arg;
        if (RUN_ENV == "WIN") {
            out_arg = out_dict + "\\args_" + to_string(width) + "_" +
                      to_string(height) + name_posts[f];
        } else {
            out_arg = out_dict + "/args_" + to_string(width) + "_" +
                      to_string(height) + name_posts[f];
        }
        ofstream outarg(out_arg, ios::app);

        for (int e = 0; e < exp_num; e++) {
            //在每轮实验开始初始化所有agent的位置
            swarm[0].set_config(0,agent_init_pos_in_shapes[f][0][0], agent_init_pos_in_shapes[f][0][1], 2);
            swarm[1].set_config(1,agent_init_pos_in_shapes[f][1][0], agent_init_pos_in_shapes[f][1][1], 2);
            swarm[2].set_config(2, agent_init_pos_in_shapes[f][2][0], agent_init_pos_in_shapes[f][2][1], 2);
            swarm[3].set_config(3, agent_init_pos_in_shapes[f][3][0], agent_init_pos_in_shapes[f][3][1], 2);
            for (int s = 4; s < agent_init_pos_in_shapes[f].size(); s++) {
                swarm[s].set_config(s, agent_init_pos_in_shapes[f][s][0], agent_init_pos_in_shapes[f][s][1], 2);
            }


            cout << "Begin allocation" << endl;
            vector<vector<double>> matrix;
            for (int m = 0; m < swarm.size(); m++) {
                vector<double> task_costs;
                for (int t = 0; t < goals.size(); t++) {
                    double cost = sqrt(pow((swarm[m].pos_x - goals[t][0]), 2) + pow((swarm[m].pos_y - goals[t][1]), 2));
                    task_costs.push_back(-cost);
                }
                matrix.push_back(task_costs);
            }
            vector<int> Lmate, Rmate;
            MinCostMatching(matrix, Lmate, Rmate);
            for (int m = 0; m < Lmate.size(); m++) {
                if(Lmate[m]>=goals.size()){
                    continue;
                }else {
                    swarm[m].set_tar(goals[Lmate[m]][0], goals[Lmate[m]][1]);
                    agent_tars[m][0] = goals[Lmate[m]][0];
                    agent_tars[m][1] = goals[Lmate[m]][1];
                }
            }
            cout << "End allocation" << endl;

            //record the initialization
            string out_name;
            if (RUN_ENV == "WIN") {
                out_name = out_dict + "\\poses_" + to_string(e) + "_" +
                           to_string(width) + "_" + to_string(height) +
                           name_posts[f];
            } else {
                out_name = out_dict + "/poses_" + to_string(e) + "_" +
                           to_string(width) + "_" + to_string(height) +
                           name_posts[f];
            }
            ofstream outfile(out_name, ios::app);
            outfile << "arguments: " << width << ' ' << height << ' ' << shape_agent_num << ' ' << THREAD_NUM
                    << endl;
            outfile << "agent positions:" << endl;
            for (int k = 0; k < shape_agent_num; k++) {
                if (k < shape_agent_num - 1) {
                    outfile << swarm[k].pos_x << ',' << swarm[k].pos_y << ' ';
                } else {
                    outfile << swarm[k].pos_x << ',' << swarm[k].pos_y;
                }
                agent_poses[swarm[k].pos_x][swarm[k].pos_y] = swarm[k].index;
            }
            outfile << endl;

            vector<thread> threads;
            int terminal = 0;
            int minor_terminal = 0;
            clock_t startT, endT;
            vector<double> dec_times;

            int ac_tar = shape_agent_num;
            int out_agent = shape_agent_num-4;
            vector<int> ac_tar_decay;
            ac_tar_decay.push_back(ac_tar);
            vector<int> out_a_decay;
            out_a_decay.push_back(out_agent);
            int re_assign_cnt = 0;
//            while (ac_tar > rand_more+3 && out_agent>rand_more+3 && terminal < 1000 and minor_terminal < 500) {
            while (ac_tar > 0 && out_agent>0 && terminal < 1000 and minor_terminal < 500) {
                cout << ac_tar <<" "<<out_agent<< endl;
                terminal += 1;
                startT = clock();
                int left = int(shape_agent_num % THREAD_NUM);
                int alloc = 0;
                int s_ids[THREAD_NUM];
                int e_ids[THREAD_NUM];
                srand(time(NULL));
                for (int k = 0; k < THREAD_NUM; k++) {
                    s_ids[k] = -1, e_ids[k] = -1;
                    if (left > alloc) {
                        s_ids[k] = k * (int(shape_agent_num / THREAD_NUM) + 1);
                        e_ids[k] = s_ids[k] + int(shape_agent_num / THREAD_NUM);
                        alloc++;
                    } else {
                        s_ids[k] = alloc * (int(shape_agent_num / THREAD_NUM) + 1) +
                                   (k - alloc) * int(shape_agent_num / THREAD_NUM);
                        e_ids[k] = s_ids[k] + int(shape_agent_num / THREAD_NUM) - 1;
                    }
                }
                bool block_flag = true;
                if (terminal > 1) {
                    int cnt = 0;
                    for (int a = 0; a < swarm.size(); a++) {
                        if ((!swarm[a].is_running) || swarm[a].not_move) {
                            cnt += 1;
                            block_flag = block_flag && true;
                        } else {
                            block_flag = block_flag && false;
                        }
                    }
                    cout << "blk: " << cnt << endl;
                } else {
                    block_flag = false;
                }

                for (int k = 0; k < THREAD_NUM; k++) {
                    threads.emplace_back(parallel_swarms, k, s_ids[k], e_ids[k], block_flag);
                }
                // 等待其他线程join
                for (int k = 0; k < THREAD_NUM; k++) {
                    threads[k].join();
                }
                threads.clear();

                if (block_flag ) {
                    cout << "A block happens, reassignment......" << endl;
                    re_assign_cnt += 1;
                    //some agents are blocked
                    vector<vector<double>> clear_matrix(0, vector<double>(0));
                    matrix.swap(clear_matrix);
                    for (int a = 0; a < swarm.size(); a++) {
                        vector<double> task_costs;
                        for (int g = 0; g < goals.size(); g++) {
                            int dis = sqrt(
                                    pow((swarm[a].pos_x - goals[g][0]), 2) + pow((swarm[a].pos_y - goals[g][1]), 2));
                            task_costs.push_back(dis);
                        }
                        matrix.push_back(task_costs);
                    }
                    Lmate.clear();
                    Rmate.clear();
                    MinCostMatching(matrix, Lmate, Rmate);
                    for (int m = 0; m < Lmate.size(); m++) {
                        if(Lmate[m]<goals.size()) {
                            swarm[m].is_running =
                                    (swarm[m].tar[0] != goals[Lmate[m]][0] || swarm[m].tar[1] != goals[Lmate[m]][1]);
                            swarm[m].set_tar(goals[Lmate[m]][0], goals[Lmate[m]][1]);
                            agent_tars[m][0] = goals[Lmate[m]][0];
                            agent_tars[m][1] = goals[Lmate[m]][1];
                        }else{
                            swarm[m].is_running =
                                    (swarm[m].tar[0] != -1 || swarm[m].tar[1] != -1);
                            swarm[m].set_tar(-1, -1);
                            agent_tars[m][0] = -1;
                            agent_tars[m][1] = -1;
                        }
                    }
                }


                ac_tar = 0;
                out_agent = 0;
                for(int a = 0; a < swarm.size(); a++){
                    if(swarm[a].is_running){
                        ac_tar++;
                    }
                    if(grids[swarm[a].pos_x][swarm[a].pos_y]==0){
                        out_agent++;
                    }
                }

                endT = clock();
                dec_times.push_back((double) (endT - startT));

                //record new positions for all agents
                outfile << "agent positions:" << endl;
                for (int k = 0; k < shape_agent_num; k++) {
                    if (k < shape_agent_num - 1) {
                        outfile << swarm[k].pos_x << ',' << swarm[k].pos_y << ' ';
                    } else {
                        outfile << swarm[k].pos_x << ',' << swarm[k].pos_y;
                    }
                }
                outfile << endl;
                ac_tar_decay.push_back(ac_tar);
                out_a_decay.push_back(out_agent);
                if (out_agent <= rand_more+6) {
                    minor_terminal += 1;
                }
            }
            double avg_t = 0;
            for (int k = 0; k < dec_times.size(); k++) {
                avg_t += dec_times[k];
            }
            double avg_iteration = avg_t / (dec_times.size() * CLOCKS_PER_SEC);
            outfile.flush();
            outfile.close();

            vector<vector<int>> formed_shape;
            for(int g=0;g<width;g++){
                formed_shape.push_back(vector<int>());
                for(int h=0;h<height;h++){
                    if(agent_poses[g][h]>=0){
                        formed_shape[g].push_back(1);
                    }else{
                        formed_shape[g].push_back(0);
                    }
                }
            }
            double mse_similarity = getSimilarity(grids,formed_shape);
            exp_avg_similarity = exp_avg_similarity + mse_similarity;
            outarg << "Experiment " << e << ":" << endl;
            outarg << "The average decision time for each iteration is: " << avg_iteration << "s." << endl;
            outarg << "Main: program exiting after " << terminal << " steps, and " << minor_terminal
                   << " steps for the last 3 positions, the similarity is " <<mse_similarity<< endl;
            outarg << "Decay ac_tar line:";
            for (int k = 0; k < ac_tar_decay.size(); k++) {
                outarg << ' ' << ac_tar_decay[k];
            }
            ac_tar_decay.clear();
            outarg << endl;
            outarg << "Decay out_tar line:";
            for (int k = 0; k < out_a_decay.size(); k++) {
                outarg << ' ' << out_a_decay[k];
            }
            out_a_decay.clear();
            outarg << endl;

            cout << "End an experiment! Clearing..." << endl;

            vector<vector<int>> array(width, vector<int>(height, -1));
            agent_poses.swap(array);
            for (int k = 0; k < swarm.size(); k++) {
                if (swarm[k].guard.owns_lock()) {
                    swarm[k].guard.unlock();
                }
            }

            if (ac_tar == rand_more+3 or out_agent == rand_more+3) {
                valid_exp += 1;
                exp_avg_iter += terminal;
                exp_avg_iter_t += avg_iteration;
            }
            if (ac_tar <= rand_more+6 or out_agent <= rand_more+6) {
                minor_valid_exp += 1;
                minor_exp_avg_iter = minor_exp_avg_iter + (terminal - minor_terminal);
                minor_exp_avg_iter_t += avg_iteration;
            }
        }
        exp_avg_iter = exp_avg_iter / valid_exp;
        exp_avg_iter_t = exp_avg_iter_t / valid_exp;
        exp_avg_similarity = exp_avg_similarity / exp_num;
        minor_exp_avg_iter = minor_exp_avg_iter / minor_valid_exp;
        minor_exp_avg_iter_t = minor_exp_avg_iter_t / minor_valid_exp;
        cout << f << ": exp_avg_iter=" << exp_avg_iter << ", exp_avg_iter_time=" << exp_avg_iter_t << endl;
        outarg << endl << endl;
        outarg << "After " << exp_num << " experiments, for shape " << f <<", avg_similarity="<<exp_avg_similarity
               <<", "<<valid_exp<<" experiments success, exp_avg_iter=" << exp_avg_iter
               << ", exp_avg_iter_time=" << exp_avg_iter_t << ", minor_exp_avg_iter=" << minor_exp_avg_iter
               << ", minor_exp_avg_iter_time=" << minor_exp_avg_iter_t << endl;
        outarg.flush();
        outarg.close();
        swarm.clear();
    }

    return 0;
}