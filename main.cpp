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

#include <algorithm>
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

namespace Munkres {

/* Utility function to print Matrix */
    template<template<typename, typename...> class Container,
            typename T,
            typename... Args>
//disable for string, which is std::basic_string<char>, a container itself
    typename std::enable_if<!std::is_convertible<Container<T, Args...>, std::string>::value &&
                            !std::is_constructible<Container<T, Args...>, std::string>::value,
            std::ostream &>::type
    operator<<(std::ostream &os, const Container<T, Args...> &con) {
        os << " ";
        for (auto &elem: con)
            os << elem << " ";

        os << "\n";
        return os;
    }

/* Handle negative elements if present. If allowed = true, add abs(minval) to
 * every element to create one zero. Else throw an exception */
    template<typename T>
    void handle_negatives(std::vector<std::vector<T>> &matrix,
                          bool allowed = true) {
        T minval = std::numeric_limits<T>::max();

        for (auto &elem: matrix)
            for (auto &num: elem)
                minval = std::min(minval, num);

        if (minval < 0) {
            if (!allowed) { //throw
                throw std::runtime_error("Only non-negative values allowed");
            } else { // add abs(minval) to every element to create one zero
                minval = abs(minval);

                for (auto &elem: matrix)
                    for (auto &num: elem)
                        num += minval;
            }
        }
    }

/* Ensure that the matrix is square by the addition of dummy rows/columns if necessary */
    template<typename T>
    void pad_matrix(std::vector<std::vector<T>> &matrix) {
        std::size_t i_size = matrix.size();
        std::size_t j_size = matrix[0].size();

        if (i_size > j_size) {
            for (auto &vec: matrix)
                vec.resize(i_size, std::numeric_limits<T>::max());
        } else if (i_size < j_size) {
            while (matrix.size() < j_size)
                matrix.push_back(std::vector<T>(j_size, std::numeric_limits<T>::max()));
        }
    }

/* For each row of the matrix, find the smallest element and subtract it from every
 * element in its row.
 * For each col of the matrix, find the smallest element and subtract it from every
 * element in its col. Go to Step 2. */
    template<typename T>
    void step1(std::vector<std::vector<T>> &matrix,
               int &step) {
        // process rows
        for (auto &row: matrix) {
            auto smallest = *std::min_element(begin(row), end(row));
            if (smallest > 0)
                for (auto &n: row)
                    n -= smallest;
        }

        // process cols
        int sz = matrix.size(); // square matrix is granted
        for (int j = 0; j < sz; ++j) {
            T minval = std::numeric_limits<T>::max();
            for (int i = 0; i < sz; ++i) {
                minval = std::min(minval, matrix[i][j]);
            }

            if (minval > 0) {
                for (int i = 0; i < sz; ++i) {
                    matrix[i][j] -= minval;
                }
            }
        }

        step = 2;
    }

/* helper to clear the temporary vectors */
    inline void clear_covers(std::vector<int> &cover) {
        for (auto &n: cover) n = 0;
    }

/* Find a zero (Z) in the resulting matrix.  If there is no starred zero in its row or
 * column, star Z. Repeat for each element in the matrix. Go to Step 3.  In this step,
 * we introduce the mask matrix M, which in the same dimensions as the cost matrix and
 * is used to star and prime zeros of the cost matrix.  If M(i,j)=1 then C(i,j) is a
 * starred zero,  If M(i,j)=2 then C(i,j) is a primed zero.  We also define two vectors
 * RowCover and ColCover that are used to "cover" the rows and columns of the cost matrix.
 * In the nested loop (over indices i and j) we check to see if C(i,j) is a zero value
 * and if its column or row is not already covered.  If not then we star this zero
 * (i.e. set M(i,j)=1) and cover its row and column (i.e. set R_cov(i)=1 and C_cov(j)=1).
 * Before we go on to Step 3, we uncover all rows and columns so that we can use the
 * cover vectors to help us count the number of starred zeros. */
    template<typename T>
    void step2(const std::vector<std::vector<T>> &matrix,
               std::vector<std::vector<int>> &M,
               std::vector<int> &RowCover,
               std::vector<int> &ColCover,
               int &step) {
        int sz = matrix.size();

        for (int r = 0; r < sz; ++r)
            for (int c = 0; c < sz; ++c)
                if (matrix[r][c] == 0)
                    if (RowCover[r] == 0 && ColCover[c] == 0) {
                        M[r][c] = 1;
                        RowCover[r] = 1;
                        ColCover[c] = 1;
                    }

        clear_covers(RowCover); // reset vectors for posterior using
        clear_covers(ColCover);

        step = 3;
    }


/* Cover each column containing a starred zero.  If K columns are covered, the starred
 * zeros describe a complete set of unique assignments.  In this case, Go to DONE,
 * otherwise, Go to Step 4. Once we have searched the entire cost matrix, we count the
 * number of independent zeros found.  If we have found (and starred) K independent zeros
 * then we are done.  If not we procede to Step 4.*/
    void step3(const std::vector<std::vector<int>> &M,
               std::vector<int> &ColCover,
               int &step) {
        int sz = M.size();
        int colcount = 0;

        for (int r = 0; r < sz; ++r)
            for (int c = 0; c < sz; ++c)
                if (M[r][c] == 1)
                    ColCover[c] = 1;

        for (auto &n: ColCover)
            if (n == 1)
                colcount++;

        if (colcount >= sz) {
            step = 7; // solution found
        } else {
            step = 4;
        }
    }

// Following functions to support step 4
    template<typename T>
    void find_a_zero(int &row,
                     int &col,
                     const std::vector<std::vector<T>> &matrix,
                     const std::vector<int> &RowCover,
                     const std::vector<int> &ColCover) {
        int r = 0;
        int c = 0;
        int sz = matrix.size();
        bool done = false;
        row = -1;
        col = -1;

        while (!done) {
            c = 0;
            while (true) {
                if (matrix[r][c] == 0 && RowCover[r] == 0 && ColCover[c] == 0) {
                    row = r;
                    col = c;
                    done = true;
                }
                c += 1;
                if (c >= sz || done)
                    break;
            }
            r += 1;
            if (r >= sz)
                done = true;
        }
    }

    bool star_in_row(int row,
                     const std::vector<std::vector<int>> &M) {
        bool tmp = false;
        for (unsigned c = 0; c < M.size(); c++)
            if (M[row][c] == 1)
                tmp = true;

        return tmp;
    }


    void find_star_in_row(int row,
                          int &col,
                          const std::vector<std::vector<int>> &M) {
        col = -1;
        for (unsigned c = 0; c < M.size(); c++)
            if (M[row][c] == 1)
                col = c;
    }


/* Find a noncovered zero and prime it.  If there is no starred zero in the row containing
 * this primed zero, Go to Step 5.  Otherwise, cover this row and uncover the column
 * containing the starred zero. Continue in this manner until there are no uncovered zeros
 * left. Save the smallest uncovered value and Go to Step 6. */
    template<typename T>
    void step4(const std::vector<std::vector<T>> &matrix,
               std::vector<std::vector<int>> &M,
               std::vector<int> &RowCover,
               std::vector<int> &ColCover,
               int &path_row_0,
               int &path_col_0,
               int &step) {
        int row = -1;
        int col = -1;
        bool done = false;

        while (!done) {
            find_a_zero(row, col, matrix, RowCover, ColCover);

            if (row == -1) {
                done = true;
                step = 6;
            } else {
                M[row][col] = 2;
                if (star_in_row(row, M)) {
                    find_star_in_row(row, col, M);
                    RowCover[row] = 1;
                    ColCover[col] = 0;
                } else {
                    done = true;
                    step = 5;
                    path_row_0 = row;
                    path_col_0 = col;
                }
            }
        }
    }

// Following functions to support step 5
    void find_star_in_col(int c,
                          int &r,
                          const std::vector<std::vector<int>> &M) {
        r = -1;
        for (unsigned i = 0; i < M.size(); i++)
            if (M[i][c] == 1)
                r = i;
    }

    void find_prime_in_row(int r,
                           int &c,
                           const std::vector<std::vector<int>> &M) {
        for (unsigned j = 0; j < M.size(); j++)
            if (M[r][j] == 2)
                c = j;
    }

    void augment_path(std::vector<std::vector<int>> &path,
                      int path_count,
                      std::vector<std::vector<int>> &M) {
        for (int p = 0; p < path_count; p++)
            if (M[path[p][0]][path[p][1]] == 1)
                M[path[p][0]][path[p][1]] = 0;
            else
                M[path[p][0]][path[p][1]] = 1;
    }

    void erase_primes(std::vector<std::vector<int>> &M) {
        for (auto &row: M)
            for (auto &val: row)
                if (val == 2)
                    val = 0;
    }


/* Construct a series of alternating primed and starred zeros as follows.
 * Let Z0 represent the uncovered primed zero found in Step 4.  Let Z1 denote the
 * starred zero in the column of Z0 (if any). Let Z2 denote the primed zero in the
 * row of Z1 (there will always be one).  Continue until the series terminates at a
 * primed zero that has no starred zero in its column.  Unstar each starred zero of
 * the series, star each primed zero of the series, erase all primes and uncover every
 * line in the matrix.  Return to Step 3.  You may notice that Step 5 seems vaguely
 * familiar.  It is a verbal description of the augmenting path algorithm (for solving
 * the maximal matching problem). */
    void step5(std::vector<std::vector<int>> &path,
               int path_row_0,
               int path_col_0,
               std::vector<std::vector<int>> &M,
               std::vector<int> &RowCover,
               std::vector<int> &ColCover,
               int &step) {
        int r = -1;
        int c = -1;
        int path_count = 1;

        path[path_count - 1][0] = path_row_0;
        path[path_count - 1][1] = path_col_0;

        bool done = false;
        while (!done) {
            find_star_in_col(path[path_count - 1][1], r, M);
            if (r > -1) {
                path_count += 1;
                path[path_count - 1][0] = r;
                path[path_count - 1][1] = path[path_count - 2][1];
            } else { done = true; }

            if (!done) {
                find_prime_in_row(path[path_count - 1][0], c, M);
                path_count += 1;
                path[path_count - 1][0] = path[path_count - 2][0];
                path[path_count - 1][1] = c;
            }
        }

        augment_path(path, path_count, M);
        clear_covers(RowCover);
        clear_covers(ColCover);
        erase_primes(M);

        step = 3;
    }

// methods to support step 6
    template<typename T>
    void find_smallest(T &minval,
                       const std::vector<std::vector<T>> &matrix,
                       const std::vector<int> &RowCover,
                       const std::vector<int> &ColCover) {
        for (unsigned r = 0; r < matrix.size(); r++)
            for (unsigned c = 0; c < matrix.size(); c++)
                if (RowCover[r] == 0 && ColCover[c] == 0)
                    if (minval > matrix[r][c])
                        minval = matrix[r][c];
    }

/* Add the value found in Step 4 to every element of each covered row, and subtract it
 * from every element of each uncovered column.  Return to Step 4 without altering any
 * stars, primes, or covered lines. Notice that this step uses the smallest uncovered
 * value in the cost matrix to modify the matrix.  Even though this step refers to the
 * value being found in Step 4 it is more convenient to wait until you reach Step 6
 * before searching for this value.  It may seem that since the values in the cost
 * matrix are being altered, we would lose sight of the original problem.
 * However, we are only changing certain values that have already been tested and
 * found not to be elements of the minimal assignment.  Also we are only changing the
 * values by an amount equal to the smallest value in the cost matrix, so we will not
 * jump over the optimal (i.e. minimal assignment) with this change. */
    template<typename T>
    void step6(std::vector<std::vector<T>> &matrix,
               const std::vector<int> &RowCover,
               const std::vector<int> &ColCover,
               int &step) {
        T minval = std::numeric_limits<T>::max();
        find_smallest(minval, matrix, RowCover, ColCover);

        int sz = matrix.size();
        for (int r = 0; r < sz; r++)
            for (int c = 0; c < sz; c++) {
                if (RowCover[r] == 1)
                    matrix[r][c] += minval;
                if (ColCover[c] == 0)
                    matrix[r][c] -= minval;
            }

        step = 4;
    }

/* Calculates the optimal cost from mask matrix */
    template<template<typename, typename...> class Container,
            typename T,
            typename... Args>
    T output_solution(const Container<Container<T, Args...>> &original,
                      const std::vector<std::vector<int>> &M) {
        T res = 0;

        for (unsigned j = 0; j < original.begin()->size(); ++j)
            for (unsigned i = 0; i < original.size(); ++i)
                if (M[i][j]) {
                    auto it1 = original.begin();
                    std::advance(it1, i);
                    auto it2 = it1->begin();
                    std::advance(it2, j);
                    res += *it2;
                    continue;
                }

        return res;
    }


/* Main function of the algorithm */
    template<template<typename, typename...> class Container,
            typename T,
            typename... Args>
//    typename std::enable_if<std::is_integral<T>::value, T>::type // Work only on integral types
    void
    hungarian(const Container<Container<T, Args...>> &original, std::vector<std::vector<int>> &assignment,
              bool allow_negatives = true) {
        /* Initialize data structures */

        // Work on a vector copy to preserve original matrix
        // Didn't passed by value cause needed to access both
        std::vector<std::vector<T>> matrix(original.size(),
                                           std::vector<T>(original.begin()->size()));

        auto it = original.begin();
        for (auto &vec: matrix) {
            std::copy(it->begin(), it->end(), vec.begin());
            it = std::next(it);
        }

        // handle negative values -> pass true if allowed or false otherwise
        // if it is an unsigned type just skip this step
        if (!std::is_unsigned<T>::value) {
            handle_negatives(matrix, allow_negatives);
        }


        // make square matrix
        pad_matrix(matrix);
        std::size_t sz = matrix.size();

        /* The masked matrix M.  If M(i,j)=1 then C(i,j) is a starred zero,
         * If M(i,j)=2 then C(i,j) is a primed zero. */
        std::vector<std::vector<int>> M(sz, std::vector<int>(sz, 0));

        /* We also define two vectors RowCover and ColCover that are used to "cover"
         *the rows and columns of the cost matrix C*/
        std::vector<int> RowCover(sz, 0);
        std::vector<int> ColCover(sz, 0);

        int path_row_0, path_col_0; //temporary to hold the smallest uncovered value

        // Array for the augmenting path algorithm
        std::vector<std::vector<int>> path(sz + 1, std::vector<int>(2, 0));

        /* Now Work The Steps */
        bool done = false;
        int step = 1;
        while (!done) {
            switch (step) {
                case 1:
                    step1(matrix, step);
                    break;
                case 2:
                    step2(matrix, M, RowCover, ColCover, step);
                    break;
                case 3:
                    step3(M, ColCover, step);
                    break;
                case 4:
                    step4(matrix, M, RowCover, ColCover, path_row_0, path_col_0, step);
                    break;
                case 5:
                    step5(path, path_row_0, path_col_0, M, RowCover, ColCover, step);
                    break;
                case 6:
                    step6(matrix, RowCover, ColCover, step);
                    break;
                case 7:
                    for (auto &vec: M) { vec.resize(original.begin()->size()); }
                    M.resize(original.size());
                    done = true;
                    break;
                default:
                    done = true;
                    break;
            }
        }

        //Printing part (optional)
//        std::cout << "Cost Matrix: \n" << original << std::endl
//                  << "Optimal assignment: \n" << M;
        assignment = M;

//        return output_solution(original, M);
    }


} // end of namespace

using namespace Munkres;
using namespace std;

const string RUN_ENV = "MAC";
const int TEMP_THREAD_NUM = 7;
#define inf 0x7f7f7f7f
//vector<Agent> swarm;
//const int width = 80;
//const int height = 80;
mutex locks[width][height];
vector<vector<int>> agent_poses(width,vector<int>(height,-1));
vector<vector<int>> agent_tars;
vector<Agent> swarm;

struct node{
    int v,w;
    node(){}
    node(int V,int W){v=V;w=W;}
    bool operator < (const node &rhs) const {
        return rhs.w>w;
    }
};


const int N_NUM = width*height; //节点总数
vector<node> G[N_NUM]; //邻接矩阵,G[i]为一个node列表,表示出发点为i,到每个node.v的距离为node.w


int dijkstra(int s, int t) {
    int dis[N_NUM];
    for(int i=1;i<=N_NUM;i++)
        dis[i]=inf;
    dis[s]=0;
    priority_queue<node> q;
    q.push( node(s,0) );
    while(!q.empty()) {
        node x=q.top();q.pop();
        for(int i=0;i<G[x.v].size();i++) {
            node y=G[x.v][i];
            if(dis[y.v]>x.w+y.w) {
                dis[y.v]=x.w+y.w;
                q.push( node(y.v,dis[y.v]) );
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

int main() {
    string dictionary, out_dict;
    if (RUN_ENV == "WIN") {
        dictionary =
                "D:\\projects\\CLionProjects\\InitSettingGenerator\\display\\" + to_string(width) + '_' +
                to_string(height);
        out_dict = "D:\\projects\\CLionProjects\\HungarianShape2\\exp";
    } else if (RUN_ENV == "MAC") {
        dictionary =
                "/Users/chuwenjie/CLionProjects/InitSettingGenerator/display/" + to_string(width) + '*' +
                to_string(height);
        out_dict = "/Users/chuwenjie/CLionProjects/HungarianShape2/exp";
    } else {
        dictionary = "./exp/" + to_string(width) + '*' + to_string(height);
        out_dict = "./exp";
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
            int shape_num = _post[5] - '0';
            shape_nums.push_back(shape_num);
            int a_num = atoi(_post.substr(7, _post.size() - 11).c_str());
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

        const int THREAD_NUM = TEMP_THREAD_NUM > a_num_s[f] ? a_num_s[f] : TEMP_THREAD_NUM;
        for (int k = 0; k < a_num_s[f]; k++) {
            swarm.push_back(Agent());
            agent_tars.push_back({-1,-1});
        }

        int exp_num = 10;
        double exp_avg_iter = 0;
        double exp_avg_iter_t = 0;
        int valid_exp = 0;
        int minor_valid_exp = 0;
        double minor_exp_avg_iter = 0;
        double minor_exp_avg_iter_t = 0;
        string out_arg;
        if (RUN_ENV == "WIN") {
            out_arg = out_dict + "\\args_" + to_string(width) + "_" +
                      to_string(height) + name_posts[f];
        } else {
            out_arg = out_dict + "/args_" + to_string(width) + "_" +
                      to_string(height) + name_posts[f];
        }
        ofstream outarg(out_arg, ios::app);

        //read the grid environment
        fstream infile;
        infile.open(filelist[f], ios::in);
        if (!infile) {
            cout << "open failed" << endl;
            exit(1);
        }
        int i = 0;
        int j = height - 1;
        vector<int> valid_l;
        vector<vector<int>> goals;
        while (!infile.eof() and j >= 0) {
            i = 0;
            while (!infile.eof() and i < width) {
                infile >> grids[i][j];
                if (grids[i][j] == 0) {
                    valid_l.push_back(i * height + j);
                } else {
                    goals.push_back({i, j});
                }
                i++;
            }
            j--;
        }
        infile.close();

        for (int e = 0; e < exp_num; e++) {
            vector<int> temp_valid_l = valid_l;
            //initialize the positions of agents
            int max_size = temp_valid_l.size();
            srand((unsigned) time(NULL));
            for (int r = 0; r < a_num_s[f]; r++) {
                int rand_p = rand() % max_size;
                int px = (int) temp_valid_l[rand_p] / height;
                int py = (int) temp_valid_l[rand_p] % height;
                swarm[r].set_config(r, px, py, 2);
                int tmp = temp_valid_l[max_size - 1];
                temp_valid_l[max_size - 1] = temp_valid_l[rand_p];
                temp_valid_l[rand_p] = tmp;
                max_size--;
            }
            //perform the goal assignment for each agent
            // work on multiple containers of the STL
            // |matrix| = M*N, M denotes the number of workers, while N denotes the number of tasks
//            vector<vector<int>> matrix;
//            for (int m = 0; m < swarm.size(); m++) {
//                vector<int> task_costs;
//                for(int t = 0; t < goals.size(); t++) {
//                    int cost = max(abs(swarm[m].pos_x-goals[t][0]),abs(swarm[m].pos_y-goals[t][1]));
//                    task_costs.push_back(cost);
//                }
//                matrix.push_back(task_costs);
//            }
            cout<<"Begin allocation"<<endl;
            vector<vector<double>> matrix;
            for (int m = 0; m < swarm.size(); m++) {
                vector<double> task_costs;
                for(int t = 0; t < goals.size(); t++) {
                    double cost = sqrt(pow((swarm[m].pos_x-goals[t][0]),2)+pow((swarm[m].pos_y-goals[t][1]),2));
                    task_costs.push_back(cost);
                }
                matrix.push_back(task_costs);
            }
            vector<vector<int>> assignment;
//            auto res = hungarian(matrix, assignment);
            hungarian(matrix, assignment);
            for(int m = 0; m < swarm.size(); m++){
                for(int t = 0; t < goals.size(); t++){
                    if(assignment[m][t] == 0)
                        continue;
                    else{
                        swarm[m].set_tar(goals[t][0], goals[t][1]);
                        agent_tars[m][0] = goals[t][0];
                        agent_tars[m][1] = goals[t][1];
                        break;
                    }
                }
            }
            cout<<"End allocation"<<endl;
//            cout << "Assignment: " << assignment << endl;
//            std::cout << "Optimal cost: " << res << std::endl;

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
            outfile << "arguments: " << width << ' ' << height << ' ' << a_num_s[f] << ' ' << THREAD_NUM
                    << endl;
            outfile << "agent positions:" << endl;
            for (int k = 0; k < a_num_s[f]; k++) {
                if (k < a_num_s[f] - 1) {
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

            int ac_tar = a_num_s[f];
            vector<int> ac_tar_decay;
            ac_tar_decay.push_back(ac_tar);
            int re_assign_cnt = 0;
            while (ac_tar > 0 && terminal < 500 and minor_terminal < 500) {
                cout << ac_tar << endl;
                terminal += 1;
                startT = clock();
                int left = int(a_num_s[f] % THREAD_NUM);
                int alloc = 0;
                int s_ids[THREAD_NUM];
                int e_ids[THREAD_NUM];
                srand(time(NULL));
                for (int k = 0; k < THREAD_NUM; k++) {
                    s_ids[k] = -1, e_ids[k] = -1;
                    if (left > alloc) {
                        s_ids[k] = k * (int(a_num_s[f] / THREAD_NUM) + 1);
                        e_ids[k] = s_ids[k] + int(a_num_s[f] / THREAD_NUM);
                        alloc++;
                    } else {
                        s_ids[k] = alloc * (int(a_num_s[f] / THREAD_NUM) + 1) +
                                   (k - alloc) * int(a_num_s[f] / THREAD_NUM);
                        e_ids[k] = s_ids[k] + int(a_num_s[f] / THREAD_NUM) - 1;
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


//                bool block_flag = true;
//                int cnt = 0;
//                for(int a = 0; a < swarm.size(); a++) {
//                    if ((!swarm[a].is_running) || swarm[a].not_move) {
//                        cnt+=1;
//                        block_flag = block_flag && true;
//                    }else{
//                        block_flag = block_flag && false;
//                    }
//                }
//                cout<<"blk: "<<cnt<<endl;
//                block_flag = false;
//                if (block_flag && re_assign_cnt < 20) {
                if (block_flag ) {
                    cout << "A block happens, reassignment......" << endl;
                    re_assign_cnt += 1;
                    //some agents are blocked
                    //需要根据当前所有agent的位置,求取最短距离,并再次分配
//                    for(int i=0;i<N_NUM;i++){
//                        G[i].clear();
//                    }
//                    for(int x=0; x<width; x++) {
//                        for(int y=0; y<height; y++) {
//                            int from = x * height + y;
//                            for (int i = -1; i <= 1; i++) {
//                                for(int j = -1; j<1; j++) {
//                                    if((!((i==0)&&(j==0))) && x+i>=0 && x+i<width && y+j>=0 && y+j<height) {
//                                        int to = (x+i) * height + (y+j);
//                                        if(agent_poses[x+i][y+j]!=-1)
//                                            G[from].push_back(node(to, inf));
//                                        else
//                                            G[from].push_back(node(to, 1));
//                                        if(agent_poses[x][y] != -1)
//                                            G[to].push_back(node(from, inf));
//                                        else
//                                            G[to].push_back(node(from, 1));
//                                    }
//                                }
//                            }
//                        }
//                    }
//                    vector<vector<int>> clear_matrix(0,vector<int>(0));
//                    matrix.swap(clear_matrix);
//                    for(int a=0; a<swarm.size(); a++) {
//                        vector<int> task_costs;
//                        int s=swarm[a].pos_x * height +swarm[a].pos_y;
//                        for(int g=0; g<goals.size(); g++) {
//                            int t = goals[g][0]*height + goals[g][1];
//                            int dis = 0;
//                            if(s!=t)
//                                dis = dijkstra(s, t);
//                            task_costs.push_back(dis);
//                        }
//                        matrix.push_back(task_costs);
//                    }

//                    vector<vector<int>> clear_matrix(0, vector<int>(0));
//                    matrix.swap(clear_matrix);
//                    for (int a = 0; a < swarm.size(); a++) {
//                        vector<int> task_costs;
//                        for (int g = 0; g < goals.size(); g++) {
//                            int dis = max(abs(swarm[a].pos_x - goals[g][0]), abs(swarm[a].pos_y - goals[g][1]));
//                            task_costs.push_back(dis);
//                        }
//                        matrix.push_back(task_costs);
//                    }
                    vector<vector<double>> clear_matrix(0, vector<double>(0));
                    matrix.swap(clear_matrix);
                    for (int a = 0; a < swarm.size(); a++) {
                        vector<double> task_costs;
                        for (int g = 0; g < goals.size(); g++) {
                            int dis = sqrt(pow((swarm[a].pos_x - goals[g][0]),2)+pow((swarm[a].pos_y - goals[g][1]),2));
                            task_costs.push_back(dis);
                        }
                        matrix.push_back(task_costs);
                    }
                    vector<vector<int>> clear_assign(0, vector<int>(0));
                    assignment.swap(clear_assign);
//                    auto res = hungarian(matrix, assignment);
                    hungarian(matrix, assignment);
                    for (int m = 0; m < swarm.size(); m++) {
                        for (int t = 0; t < goals.size(); t++) {
                            if (assignment[m][t] == 0)
                                continue;
                            else {
//                                if (!swarm[m].is_running &&
//                                    (swarm[m].tar[0] != goals[t][0] || swarm[m].tar[1] != goals[t][1])) {
//                                    swarm[m].is_running = true;
//                                } else if (!swarm[m].is_running && swarm[m].tar[0] == goals[t][0] &&
//                                           swarm[m].tar[1] == goals[t][1]) {
//                                    swarm[m].is_running = false;
//                                }
                                if(swarm[m].tar[0] != goals[t][0] || swarm[m].tar[1] != goals[t][1]) {
                                    swarm[m].is_running = true;
                                }else{
                                    swarm[m].is_running = false;
                                }
                                swarm[m].set_tar(goals[t][0], goals[t][1]);
                                agent_tars[m][0] = goals[t][0];
                                agent_tars[m][1] = goals[t][1];
                                break;
                            }
                        }
                    }
                }


                ac_tar = 0;
                for(int a = 0; a < swarm.size(); a++){
                    if(swarm[a].is_running){
                        ac_tar++;
                    }
                }

                endT = clock();
                dec_times.push_back((double) (endT - startT));

                //record new positions for all agents
                outfile << "agent positions:" << endl;
                for (int k = 0; k < a_num_s[f]; k++) {
                    if (k < a_num_s[f] - 1) {
                        outfile << swarm[k].pos_x << ',' << swarm[k].pos_y << ' ';
                    } else {
                        outfile << swarm[k].pos_x << ',' << swarm[k].pos_y;
                    }
                }
                outfile << endl;
                ac_tar_decay.push_back(ac_tar);
                if (ac_tar <= 3) {
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

            outarg << "Experiment " << e << ":" << endl;
            outarg << "The average decision time for each iteration is: " << avg_iteration << "s." << endl;
            outarg << "Main: program exiting after " << terminal << " steps, and " << minor_terminal
                   << " steps for the last 3 positions." << endl;
            outarg << "Decay line:";
            for (int k = 0; k < ac_tar_decay.size(); k++) {
                outarg << ' ' << ac_tar_decay[k];
            }
            ac_tar_decay.clear();
            outarg << endl;

            cout << "End an experiment! Clearing..." << endl;

            vector<vector<int>> array(width, vector<int>(height, -1));
            agent_poses.swap(array);
            for (int k = 0; k < swarm.size(); k++) {
                if (swarm[k].guard.owns_lock()) {
                    swarm[k].guard.unlock();
                }
            }

            if (ac_tar == 0) {
                valid_exp += 1;
                exp_avg_iter += terminal;
                exp_avg_iter_t += avg_iteration;
            }
            if (ac_tar <= 3) {
                minor_valid_exp += 1;
                minor_exp_avg_iter = minor_exp_avg_iter + (terminal - minor_terminal);
                minor_exp_avg_iter_t += avg_iteration;
            }
        }
        exp_avg_iter = exp_avg_iter / valid_exp;
        exp_avg_iter_t = exp_avg_iter_t / valid_exp;
        minor_exp_avg_iter = minor_exp_avg_iter / minor_valid_exp;
        minor_exp_avg_iter_t = minor_exp_avg_iter_t / minor_valid_exp;
        cout << f << ": exp_avg_iter=" << exp_avg_iter << ", exp_avg_iter_time=" << exp_avg_iter_t << endl;
        outarg << endl << endl;
        outarg << "After " << exp_num << " experiments, for shape " << f << ", exp_avg_iter=" << exp_avg_iter
               << ", exp_avg_iter_time=" << exp_avg_iter_t << ", minor_exp_avg_iter=" << minor_exp_avg_iter
               << ", minor_exp_avg_iter_time=" << minor_exp_avg_iter_t << endl;
        outarg.flush();
        outarg.close();
        swarm.clear();
    }

    return 0;
}