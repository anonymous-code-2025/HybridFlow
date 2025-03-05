#include "graphoperations.h"
#include <memory.h>
#include <queue>

void GraphOperations::getKCore(const Graph *graph, int *core_table) {
    int vertices_count = graph->getVerticesCount();
    int max_degree = graph->getGraphMaxDegree();

    int* vertices = new int[vertices_count];          // Vertices sorted by degree.
    int* position = new int[vertices_count];          // The position of vertices in vertices array.
    int* degree_bin = new int[max_degree + 1];      // Degree from 0 to max_degree.
    int* offset = new int[max_degree + 1];          // The offset in vertices array according to degree.

    std::fill(degree_bin, degree_bin + (max_degree + 1), 0);

    for (int i = 0; i < vertices_count; ++i) {
        int degree = graph->getVertexDegree(i);
        core_table[i] = degree;
        degree_bin[degree] += 1;
    }

    int start = 0;
    for (int i = 0; i < max_degree + 1; ++i) {
        offset[i] = start;
        start += degree_bin[i];
    }

    for (int i = 0; i < vertices_count; ++i) {
        int degree = graph->getVertexDegree(i);
        position[i] = offset[degree];
        vertices[position[i]] = i;
        offset[degree] += 1;
    }

    for (int i = max_degree; i > 0; --i) {
        offset[i] = offset[i - 1];
    }
    offset[0] = 0;

    for (int i = 0; i < vertices_count; ++i) {
        int v = vertices[i];

        ui count;
        const VertexID * neighbors = graph->getVertexNeighbors(v, count);

        for(int j = 0; j < static_cast<int>(count); ++j) {
            int u = neighbors[j];

            if (core_table[u] > core_table[v]) {

                // Get the position and vertex which is with the same degree
                // and at the start position of vertices array.
                int cur_degree_u = core_table[u];
                int position_u = position[u];
                int position_w = offset[cur_degree_u];
                int w = vertices[position_w];

                if (u != w) {
                    // Swap u and w.
                    position[u] = position_w;
                    position[w] = position_u;
                    vertices[position_u] = w;
                    vertices[position_w] = u;
                }

                offset[cur_degree_u] += 1;
                core_table[u] -= 1;
            }
        }
    }

    delete[] vertices;
    delete[] position;
    delete[] degree_bin;
    delete[] offset;
}


void GraphOperations::dfs(TreeNode *tree, VertexID cur_vertex, VertexID *dfs_order, ui &count) {
    dfs_order[count++] = cur_vertex;

    for (ui i = 0; i < tree[cur_vertex].children_count_; ++i) {
        dfs(tree, tree[cur_vertex].children_[i], dfs_order, count);
    }
}

void GraphOperations::compute_degeneracy_order(const Graph *graph, uint32_t *degeneracy_order) {
    int vertices_count = graph->getVerticesCount();
    int max_degree = graph->getGraphMaxDegree();

    int* core_table = new int[vertices_count];        // core values.
    int* vertices = new int[vertices_count];          // Vertices sorted by degree.
    int* position = new int[vertices_count];          // The position of vertices in vertices array.
    int* degree_bin = new int[max_degree + 1];      // Degree from 0 to max_degree.
    int* offset = new int[max_degree + 1];          // The offset in vertices array according to degree.

    std::fill(degree_bin, degree_bin + (max_degree + 1), 0);

    for (int i = 0; i < vertices_count; ++i) {
        int degree = graph->getVertexDegree(i);
        core_table[i] = degree;
        degree_bin[degree] += 1;
    }

    int start = 0;
    for (int i = 0; i < max_degree + 1; ++i) {
        offset[i] = start;
        start += degree_bin[i];
    }

    for (int i = 0; i < vertices_count; ++i) {
        int degree = graph->getVertexDegree(i);
        position[i] = offset[degree];
        vertices[position[i]] = i;
        offset[degree] += 1;
    }

    for (int i = max_degree; i > 0; --i) {
        offset[i] = offset[i - 1];
    }
    offset[0] = 0;

    for (int i = 0; i < vertices_count; ++i) {
        int v = vertices[i];

        ui count;
        const VertexID * neighbors = graph->getVertexNeighbors(v, count);

        for(int j = 0; j < static_cast<int>(count); ++j) {
            int u = neighbors[j];

            if (core_table[u] > core_table[v]) {

                // Get the position and vertex which is with the same degree
                // and at the start position of vertices array.
                int cur_degree_u = core_table[u];
                int position_u = position[u];
                int position_w = offset[cur_degree_u];
                int w = vertices[position_w];

                if (u != w) {
                    // Swap u and w.
                    position[u] = position_w;
                    position[w] = position_u;
                    vertices[position_u] = w;
                    vertices[position_w] = u;
                }

                offset[cur_degree_u] += 1;
                core_table[u] -= 1;
            }
        }

        degeneracy_order[i] = v;
    }

    delete[] core_table;
    delete[] vertices;
    delete[] position;
    delete[] degree_bin;
    delete[] offset;
}

// 函数，计算自同构结构，第一个参数表示查询图，第二个参数用来存自同构结构
void GraphOperations::compute_automorphism(const Graph *graph, std::vector<std::vector<uint32_t>> &embeddings) {
    // Note that this method is working on small graphs (tens of vertices) only.
    // Initialize resource.
    uint32_t n = graph->getVerticesCount(); // 获取图的顶点数量
    std::vector<bool> visited(n, false);    // 标志访问状态，自同构映射过程中是否已经访问过了
    std::vector<uint32_t> idx(n);           // 存储索引
    std::vector<uint32_t> mapping(n);       // 存储映射关系
    std::vector<std::vector<uint32_t>> local_candidates(n); // 存放局部候选集
    std::vector<std::vector<uint32_t>> global_candidates(n);// 存放全局候选集
    std::vector<std::vector<uint32_t>> backward_neighbors(n);// 存储反向邻居，用于计算local_candidates

    // Initialize global candidates.
    // 初始化 全局候选集，遍历查询图中的所有点，对于每一个u，将具有相同标签且度数大于等于u的顶点v, 添加到全局候选集
    // 这里本质为自己向自己映射，有点像Q-D的映射
    for (uint32_t u = 0; u < n; ++u) {
        uint32_t u_label = graph->getVertexLabel(u);
        uint32_t u_degree = graph->getVertexDegree(u);

        for (uint32_t v = 0; v < n; ++v) {
            uint32_t v_label = graph->getVertexLabel(v);
            uint32_t v_degree = graph->getVertexDegree(v);

            if (v_label == u_label && v_degree >= u_degree)
                global_candidates[u].push_back(v);
        }
    }

    // Generate a matching order.
    // 生成一个匹配顺序

    // 选择一个顶点作为匹配顺序的起点，选择的顶点是具有最少候选集的顶点，，减小搜索空间
    std::vector<uint32_t> matching_order;
    uint32_t selected_vertex = 0;
    uint32_t selected_vertex_selectivity = global_candidates[selected_vertex].size();
    for (uint32_t u = 1; u < n; ++u) {
        if (global_candidates[u].size() < selected_vertex_selectivity){
            selected_vertex = u;
            selected_vertex_selectivity = global_candidates[u].size();
        }
    }
    // 把起始的第一个点 先放入匹配顺序中，同时标记为已经访问
    matching_order.push_back(selected_vertex);
    visited[selected_vertex] = true; 

    // 遍历剩余的n-1个点，为每一点都确定匹配顺序
    for (uint32_t i = 1; i < n; ++i) {
        // 每次初始化成 一个比任何可能的候选集都大的值
        selected_vertex_selectivity = n + 1; 
        // 再遍历这n个点，尝试找到下一个匹配的点
        for (uint32_t u = 0; u < n; ++u) {
            if (!visited[u]) { // 如果u没被访问过，即它还没有被选为匹配顺序的一部分
                bool is_feasible = false; // 标记u是否可以成为匹配顺序的下一个点

                uint32_t u_nbr_count; // 获取它的邻居
                auto u_nbr = graph->getVertexNeighbors(u, u_nbr_count);
                for (uint32_t j = 0; j < u_nbr_count; ++j) { // 如果它有一个邻居在匹配顺序里边了，才有可能
                    uint32_t uu = u_nbr[j];

                    if (visited[uu]) {
                        is_feasible = true;
                        break;
                    }
                }
                // 如果它有一个被标记的邻居，并且它的候选集小于当前迭代的候选集，即u成为下一个匹配的点
                if (is_feasible && global_candidates[u].size() < selected_vertex_selectivity) {
                    selected_vertex = u;
                    selected_vertex_selectivity = global_candidates[u].size();
                }
            }
        }
        matching_order.push_back(selected_vertex);
        visited[selected_vertex] = true;
    }
    // 匹配顺序至此生成完毕
    
    
    // 将标志位全部清空
    std::fill(visited.begin(), visited.end(), false);

    // *设置反向邻居，它存放的是 在匹配顺序中，排在当前顶点之前的顶点，并且与当前顶点之间相连（一跳前驱）
    // Set backward neighbors to compute local candidates.
    for (uint32_t i = 1; i < n; ++i) {
        uint32_t u = matching_order[i];
        for (uint32_t j = 0; j < i; ++j) {
            uint32_t uu = matching_order[j];

            if (graph->checkEdgeExistence(uu, u)) {
                backward_neighbors[u].push_back(uu);
            }
        }
    }


    // *Recursive search along the matching order.
    // 沿着制定后的匹配顺序递归搜索
    
    int cur_level = 0; // 记录当前的层级
    // 将全局候选集 global_candidates 中与匹配顺序中的第一个顶点对应的候选集赋值给 local_candidates 中的当前层级
    local_candidates[cur_level] = global_candidates[matching_order[0]];
    
    // 开始递归搜索所有可能的自同构结构
    while (true) {
        while (idx[cur_level] < local_candidates[cur_level].size()) {//遍历当前层级（当前点）的所有局部候选点，idx[]用来标记当前的搜索位置
            uint32_t u = matching_order[cur_level]; // u是匹配顺序中的当前顶点
            uint32_t v = local_candidates[cur_level][idx[cur_level]]; // v是u的第一个候选顶点
            idx[cur_level] += 1;

            if (cur_level == n - 1) { // 到达最后一层
                // Find an embedding
                mapping[u] = v;
                embeddings.push_back(mapping); // 找到一个mapping自同构结构，存到embedding中
            } else { // 没到最后一层
                mapping[u] = v;   // 将当前顶点u映射到v
                visited[v] = true;// v标记为已访问
                cur_level += 1;   // 遍历到下一层级（下一个顶点）
                idx[cur_level] = 0;// 下一层索引 重置到0号

                {
                    // Compute local candidates.
                    // 计算当前层级顶点的局部候选集
                    u = matching_order[cur_level]; // 沿着当前的匹配顺序，现在顶点为u
                    for (auto temp_v: global_candidates[u]) { // 遍历u的候选点
                        if (!visited[temp_v]) { // 如果v没有被访问过
                            bool is_feasible = true; // 标记temp_v是否能成为局部候选点
                            for (auto uu: backward_neighbors[u]) { // 找到u的所有前一跳邻居
                                uint32_t edge_label = graph->getEdgeLabelByVertex(u, uu); // 获取边标签
                                // 如果 temp_v 与所有反向邻居都兼容，那么它是可行的。
                                uint32_t temp_vv = mapping[uu];
                                // TODO: check edge label, 如果边不存在或者边边标签不符
                                if (!graph->checkEdgeExistence(temp_v, temp_vv)
                                    || edge_label != graph->getEdgeLabelByVertex(temp_v, temp_vv)) {
                                    is_feasible = false;
                                    continue;
                                }
                            }
                            if (is_feasible)
                                local_candidates[cur_level].push_back(temp_v);
                        }
                    }
                }
            }
        }
        
        local_candidates[cur_level].clear();
        cur_level -= 1;
        if (cur_level < 0) {
            break;
        }
        visited[mapping[matching_order[cur_level]]] = false;
    }
}

