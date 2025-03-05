#include <queue>
#include <computesetintersection.h>
#include "order_manager.h"
#include "graphoperations.h"

void OrderManager::initialize(const Graph *query_graph) {}

void OrderManager::release() {
    automorphisms_.clear();
    orders_.clear();
    edge_orders_mapping_.clear();
    label_edge_mapping_.clear();
    label_automorphism_mapping_.clear();
    automorphism_edges_.clear();
}

// 函数：获取 特定查询边的 顺序信息
OrdersPerEdge *OrderManager::get_orders(Edge edge) {
    auto it = edge_orders_mapping_.find(edge);
    if (it != edge_orders_mapping_.end()) {
        return &orders_[it->second];
    }

    return nullptr;
}
// 函数：获取指定标签组的自同构集合
std::vector<uint32_t> *OrderManager::get_mapped_automorphism(LabelTriple label_triple) {
    auto it = label_automorphism_mapping_.find(label_triple);
    if (it != label_automorphism_mapping_.end()) {
        return &it->second;
    }

    return nullptr;
}
// 函数：获取指定标签组 的 边集合
std::vector<Edge> *OrderManager::get_mapped_edges(LabelTriple label_triple) {
    auto it = label_edge_mapping_.find(label_triple);
    if (it != label_edge_mapping_.end()) {
        return &it->second;
    }
    return nullptr;
}

// 函数：获取与特定自同构ID相关的元数据， 包含 对应的边，以及边的顺序信息指针， 和自同构的数量
std::tuple<Edge, OrdersPerEdge*, uint32_t> OrderManager::get_automorphism_meta(uint32_t id) {
    return std::make_tuple(automorphism_edges_[id][0], &orders_[id], (uint32_t)(automorphism_edges_[id].size()));
}

// 函数：检测自同构边，由create_order()函数调用，即在图的自同构变换下保持不变的边。自同构边完成了分组
/*
[
  [ {u0, u1}, {u1, u2}, {u2, u3}, {u3, u0} ],  // 所有边在恒等自同构下
  [ {u0, u1}, {u2, u3} ]               // {u0, u1} 和 {u2, u3} 在非恒等自同构下
]
*/
void OrderManager::detect_automorphism_edges(const Graph *query_graph) {
    // Divide the vertex into disjoint sets based on automorphisms.
    // 调用自定义的图操作工具包，计算自同构结构
    GraphOperations::compute_automorphism(query_graph, automorphisms_);// automorphism中的每一个元素都是一个自同构，表示为顶点索引的排列

    uint32_t n = query_graph->getVerticesCount();
    spp::sparse_hash_set<Edge> selected;// 存储已经选择的边
    // test by gz
    // printf("\n starting test automorphisms\n");
    /*printf("\n************\n");
    for(auto embedding: automorphisms_)
    {
        for(int i = 0; i < embedding.size(); i++)
        {
            printf("%d ", embedding[i]);
        }
        printf("\n");
    }
    printf("\n************\n");*/
    
    // 循环，遍历查询图中的所有点
    for (uint32_t u = 0; u < n; ++u) {
        uint32_t u_nbr_count;
        auto u_nbr = query_graph->getVertexNeighbors(u, u_nbr_count);
        // 对于当前节点u，遍历它的每一个邻居
        for (uint32_t i = 0; i < u_nbr_count; ++i) {
            uint32_t uu = u_nbr[i];
            Edge e = {u, uu};
            if (!selected.contains(e)) { // 如果当前边u--uu还没有被选择过
                selected.insert(e); // 标记选择，
                automorphism_edges_.push_back({e}); //初始边只包含e，{}变向量
                for (auto &embedding: automorphisms_) { // 遍历自同构结构集合
                    Edge mapped_e = {embedding[u], embedding[uu]}; // 对于每个自同构结构，创建一个映射边（顶点u和顶点uu在自同构下的映射）
                    if (!selected.contains(mapped_e)) { // 如果该映射边没有被选择过
                        selected.insert(mapped_e); 
                        automorphism_edges_.back().push_back(mapped_e); // 并将该映射边添加到automorphism_edges的最后一个自同构边集合
                    }
                }
            }
        }
    }
    /*printf("\n 测试自同构边集\n");
    printf("\n************\n");
    for(auto embedding: automorphism_edges_)
    {
        for(int i = 0; i < embedding.size(); i++)
        {
            printf("%d %d   ", embedding[i].first, embedding[i].second);
        }
        printf("\n");
    }
    printf("\n************\n");*/
    
}

// ***函数 为缩减后的查询图创建索引顺序 （这里引用传递是会改变原本主函数中order的东西的），手推完毕
void OrderManager::create_indexing_order_for_reduced_graph(const Graph *query_graph, Edge edge, OrdersPerEdge &order) {
    uint32_t n = query_graph->getVerticesCount(); // 获取查询图顶点的数量
    // 获取order对象中存储索引顺序和反向邻居信息的成员变量的引用
    auto& indexing_order = order.indexing_order_; 
    auto& indexing_order_bn = order.indexing_order_bn_;
    auto& indexing_order_bn_offset = order.indexing_order_bn_offset_;

    
    indexing_order = {edge.first, edge.second}; // 初始化索引顺序，开始时只包含给定边的两个顶点
    indexing_order_bn = {edge.first}; // 初始化反向邻居列表，开始时只包含给定边的一个端点
    indexing_order_bn_offset = {0, 0, 1}; //初始化反向邻居列表的偏移量，开始时只有一个元素，值为1，表示反向邻居列表的大小

    std::vector<bool> visited(n);  // 设置标志位，给定边的两个端点已访问
    visited[edge.first] = true;     
    visited[edge.second] = true;

    // 初始化三个向量，用于存储给定边的两个端点相邻的顶点
    std::vector<uint32_t> adjacent_to_both; // 和给定边的两个点 都是邻居
    std::vector<uint32_t> adjacent_to_one;  // 和给定边的其中一个点 是邻居
    std::vector<uint32_t> adjacent_to_none; // 和给定边的两个点 都不是邻居
    // 遍历Q中的每个点
    for (uint32_t u = 0; u < n; ++u) {
        if (u == edge.first || u == edge.second)
            continue;
        // 检查当前点u是否与给定边的两个端点存在边
        bool f1 = query_graph->checkEdgeExistence(u, edge.first);
        bool f2 = query_graph->checkEdgeExistence(u, edge.second);

        if (f1 && f2) {
            adjacent_to_both.push_back(u);
        }
        else if (!f1 && !f2) {
            adjacent_to_none.push_back(u);
        }
        else {
            adjacent_to_one.push_back(u);
        }
    }
    // 定义一个lambda函数
    auto generate_function = [query_graph, &visited, &indexing_order, &indexing_order_bn, &indexing_order_bn_offset]
            (std::vector<uint32_t>& target_vertex) {
        // 循环遍历target_vertex向量中的所有顶点
        for (uint32_t i = 0; i < target_vertex.size(); ++i) {  
            std::vector<uint32_t> current_vertex_bn; // 当前点的反向邻居
            std::vector<uint32_t> selected_vertex_bn;// 选定顶点的反向邻居
            uint32_t selected_vertex;// 当前选定点

            // 找到
            for (auto u : target_vertex) { // 遍历target中的所有顶点
                if (!visited[u]) { // 如果没被访问过
                    current_vertex_bn.clear(); // 当前点的反向邻居 清空
                    // Get backward neighbors.
                    // 获取反向邻居
                    for (auto uu : indexing_order) { // 遍历indexing_order向量
                        if (query_graph->checkEdgeExistence(u, uu)) {
                            current_vertex_bn.push_back(uu); // 如果存在边，就
                        }
                    }
                    if (current_vertex_bn.size() > selected_vertex_bn.size()) {
                        current_vertex_bn.swap(selected_vertex_bn);
                        selected_vertex = u;
                    }
                }
            }

            indexing_order.push_back(selected_vertex);
            indexing_order_bn.insert(indexing_order_bn.end(), selected_vertex_bn.begin(), selected_vertex_bn.end());
            indexing_order_bn_offset.push_back(indexing_order_bn.size());

            visited[selected_vertex] = true;
        }
    };
    order.triangle_end_index_ = adjacent_to_both.size() + 2;
    order.adjacent_end_index_ = adjacent_to_both.size() + adjacent_to_one.size() + 2;

    generate_function(adjacent_to_both);
    generate_function(adjacent_to_one);
    generate_function(adjacent_to_none);
    // tested by gz

}
// 为什么要设置强连通分量，因为在为一条边创建匹配顺序时，是不需要考虑这连个点的，这条边可能会把整个查询图划分为 好几个连通分量
// ***函数：为缩减后的查询图的每一条边 都创建一个匹配顺序（获取所有的连通分量 + 为每一个连通分量创建匹配顺序 + 合并这些连通分量）
void OrderManager::create_matching_order_for_reduced_graph(const Graph *query_graph, Edge edge, OrdersPerEdge &order) {
    /**
     * 1. Get all connected components.
     * 2. Generate the matching order for each connected component.
     * 3. Merge the connected component.
     */
    // 获取查询图的顶点数量，并把这条边的两个顶点的 标志位设为true
     uint32_t n = query_graph->getVerticesCount();
     std::vector<bool> visited(n, false);
     visited[edge.first] = true;
     visited[edge.second] = true;

     std::vector<std::vector<uint32_t>> connected_components;
     std::queue<uint32_t> q;  // 设置一个队列q，用于广度优先搜索去生成连通分量

     // 1. Generate connected components by conducting a BFS from each vertex.
     // 广度优先搜索BFS 用来 生成匹配顺序
     for (uint32_t u = 0; u < n; ++u) {
         if (!visited[u]) { // 如果u没被访问过，说明找到了一个新的连通分量
             std::vector<uint32_t> component;
             component.push_back(u);
             q.push(u);
             visited[u] = true; // 把u压入队列， 标志位设置为true

             while (!q.empty()) { // 只要队列不空，基操
                 uint32_t uu = q.front();
                 q.pop(); // 取出元素 + 弹队列
                 uint32_t uu_nbrs_count;
                 auto uu_nbrs = query_graph->getVertexNeighbors(uu, uu_nbrs_count);
                // 遍历u的一跳邻居
                 for (uint32_t i = 0; i < uu_nbrs_count; ++i) {
                     uint32_t uuu = uu_nbrs[i];
                     if (!visited[uuu]) {
                         q.push(uuu);
                         component.push_back(uuu);
                         visited[uuu] = true;
                     }
                 }
             }

             connected_components.emplace_back(component);
         }
     }

     // 2. Generate the matching order for each connected component.
     // 2. 为每一个连通分量，生成一个匹配顺序
     std::vector<std::vector<uint32_t>> matching_orders;

     // map the old vertex id to new vertex id
     std::vector<uint32_t> reverse_mapping(n, 0);

     std::vector<Graph*> graphs;// 用于存储每个连通分量
     // 遍历所有连通分量，为每个连通分量计算匹配顺序
     for (auto& component : connected_components) {
         if (component.size() == 1 || component.size() == 2) { // 如果其中一个连通分量的大小为1或2，之间将这两个点作为匹配顺序
             matching_orders.push_back(component);
             graphs.push_back(nullptr);
         }
         else {
             for (int i = 0; i < component.size(); ++i) {
                reverse_mapping[component[i]] = i; // ****重新编号的目的***，在为当前连通分量生成匹配顺序时候，只需要考虑当前分量的相关点和边，将旧的顶点ID 映射到新的顶点ID
             }

             // Create vertex list (vertex_id, vertex_label_id)  and edge list (begin_vertex_id, end_vertex_id, edge_label_id)
             std::vector<std::pair<uint32_t, uint32_t>> vertex_list;
             std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> edge_list;
             // 遍历该连通分量重的每一个顶点
             for (auto u: component) {
                 uint32_t u_id = reverse_mapping[u]; // 获取u的新ID
                 uint32_t u_label = query_graph->getVertexLabel(u); // 获取u的点标签
                 vertex_list.emplace_back(u_id, u_label); // 存入vertex_list

                 uint32_t u_nbrs_count;
                 auto u_nbrs = query_graph->getVertexNeighbors(u, u_nbrs_count);
                 // 遍历u的每一个邻居
                 for (uint32_t i = 0; i < u_nbrs_count; ++i) {
                     uint32_t uu = u_nbrs[i];
                     if (uu != edge.first && uu != edge.second) { // 如果u不是当前的两个点
                         uint32_t uu_id = reverse_mapping[uu]; // 获取uu的新ID

                         if (u_id < uu_id) { 
                             // 如果顶点 u 的 ID 小于邻居 uu 的 ID，将边 (u_id, uu_id) 及其标签添加到边列表中
                             edge_list.emplace_back(u_id, uu_id, query_graph->getEdgeLabelByVertex(u_id, uu_id));
                         }
                     }
                 }
             }
             
             // 创建一个新图来表示当前连通分量，并从顶点列表和边列表中加载图数据
             Graph *graph = new Graph(false);
             graph->is_edge_labeled = true;
             graph->loadGraphFromMemory(vertex_list, edge_list);
             // 为当前连通分量生成匹配顺序
             std::vector<uint32_t> matching_order;
             generate_matching_order_with_RI(graph, matching_order); /*****************************/
             matching_orders.push_back(matching_order);
             graphs.push_back(graph);
         }
     }
    
     // 3. Merge the matching orders into one based on 1) the 2-core size; and 2) the component size.
     // 基于 两核大小和连通分量块的大小 合并 匹配顺序
     auto& updated_matching_order = order.matching_order_;
     auto& updated_matching_order_bn = order.matching_order_bn_;  // 记录前驱邻居
     auto& updated_matching_order_bn_offset = order.matching_order_bn_offset_; // 前驱邻居的偏移量
     auto& updated_matching_order_view_mapping = order.matching_order_view_mappings_; // 匹配顺序和视图的映射
     auto& updated_matching_order_edge_type = order.matching_order_edge_type_; // 引用边的类型

     updated_matching_order_bn_offset = {0};

     std::vector<uint32_t> core_size; // 存储每个子图的两核 大小
     std::vector<uint32_t> graph_size; // 存储每个子图的大小

     // Initialize the core size and graph size for each subgraph.
     // 初始化子图的两核大小和图大小
     for (uint32_t i = 0; i < graphs.size(); ++i) {
         Graph* graph = graphs[i];
         if (graph != nullptr) { 
             uint32_t size = 0;
             std::vector<int> core(graph->getVerticesCount(), 0);
             GraphOperations::getKCore(graph, core.data());/****************************/

             for (auto core_value: core) {
                 if (core_value >= 2) { // 如果顶点的两核值大于等于2
                     size += 1;// 两核心大小计数+1
                 }
             }
             core_size.push_back(size);
             graph_size.push_back(graph->getVerticesCount());
         }
         else {
             core_size.push_back(0);
             graph_size.push_back(connected_components[i].size());
         }
     }
     // 重置访问标记数组
     std::fill(visited.begin(), visited.end(), false);
     // 遍历所有子图
     for (uint32_t i = 0; i < graphs.size(); ++i) {
         // Pick the graph.
         uint32_t selected_core_value = 0;
         uint32_t selected_graph_size = 0;
         uint32_t selected_graph = 0;
         // 选择一个子图（选择标准：核的大小越大 或者 核值相同但图更大）
         for (uint32_t j = 0; j < graphs.size(); ++j) {
             if (!visited[j]) {
                 if ((core_size[j] > selected_core_value)
                     || (core_size[j] == selected_core_value && graph_size[j] > selected_graph_size)) {
                     selected_graph = j;
                     selected_core_value = core_size[j];
                     selected_graph_size = graph_size[j];
                 }
             }
         }

         visited[selected_graph] = true; // 设置图已经被访问

         // Update matching order and backward neighbors.
         // 更新 匹配顺序和前驱邻居
         if (graph_size[selected_graph] == 1) { // 如果选择的子图只有一个顶点
             uint32_t u = connected_components[selected_graph][0]; // 获取该顶点
             updated_matching_order.push_back(u); // 将该点放入匹配顺序
             if (query_graph->getVertexDegree(u) == 1) { // 如果顶点度为1，则判断把边的哪个顶点加入前驱邻居列表
                 if (query_graph->checkEdgeExistence(u, edge.first)) {
                     updated_matching_order_bn.push_back(edge.first);
                 } else {
                     updated_matching_order_bn.push_back(edge.second);
                 }
             }
             updated_matching_order_edge_type.push_back(RelationEdgeType::REGULAR); // 设置边的类型为Regular
             updated_matching_order_bn_offset.push_back(updated_matching_order_bn.size()); // 更新前驱邻居的偏移量
         } else if (graph_size[selected_graph] == 2) { // 如果该子图（连通分量）有两个顶点，则选择度更大的那个点，放在前面
             if (query_graph->getVertexDegree(connected_components[selected_graph][0]) < query_graph->getVertexDegree(connected_components[selected_graph][1])) {
                 std::swap(connected_components[selected_graph][0], connected_components[selected_graph][1]);
             }
             // 将连通分量的两个顶点加入匹配顺序中
             updated_matching_order.insert(updated_matching_order.end(), connected_components[selected_graph].begin(), connected_components[selected_graph].end());
             // Insert the first vertex.
             updated_matching_order_bn_offset.push_back(updated_matching_order_bn.size()); // 放入偏移量

             // Insert the second vertex.
             updated_matching_order_bn.push_back(connected_components[selected_graph][0]);
             updated_matching_order_edge_type.push_back(RelationEdgeType::REGULAR);
             updated_matching_order_bn_offset.push_back(updated_matching_order_bn.size());
         } else { // 如果该子图(连通分量)有多个顶点
             auto graph = graphs[selected_graph];
             std::vector<bool> vertex_visited(matching_orders[selected_graph].size(), false);
             for (auto u: matching_orders[selected_graph]) { // 从该子图(连通分量)确定好的匹配顺序中，遍历每个点
                 uint32_t u_nbrs_count;
                 auto u_nbrs = graph->getVertexNeighbors(u, u_nbrs_count);

                 for (uint32_t j = 0; j < u_nbrs_count; ++j) { // 遍历该点的邻居
                     uint32_t uu = u_nbrs[j];
                     if (vertex_visited[uu]) {
                         updated_matching_order_bn.push_back(connected_components[selected_graph][uu]);
                         updated_matching_order_edge_type.push_back(RelationEdgeType::REGULAR);
                     }
                 }

                 updated_matching_order.push_back(connected_components[selected_graph][u]);
                 updated_matching_order_bn_offset.push_back(updated_matching_order_bn.size());
                 vertex_visited[u] = true;
             }
         }
     }

    updated_matching_order_view_mapping.resize(updated_matching_order_bn.size());

     // Release graph.
     for (auto graph : graphs) {
         delete graph;
     }
}

// ***函数：为每一个每一个连通分量（这里用graph来存储）生成匹配顺序，存放到matching_order中
void OrderManager::generate_matching_order_with_RI(const Graph *graph, std::vector<uint32_t> &matching_order) {
    
    // 获取当前连通分量的点的数量以及设置标志位 
    uint32_t n = graph->getVerticesCount();
    std::vector<bool> visited(n, false);
    // Select the vertex with the maximum degree as the start vertex.
    // 第一步是选择具有最大度的顶点
    uint32_t selected_vertex = 0;
    uint32_t selected_vertex_selectivity = graph->getVertexDegree(selected_vertex);
    for (ui u = 1; u < n; ++u) {
        uint32_t u_selectivity = graph->getVertexDegree(u);
        if (u_selectivity > selected_vertex_selectivity) {
            selected_vertex = u;
            selected_vertex_selectivity = u_selectivity;
        }
    }
    
    matching_order.push_back(selected_vertex); // 整个连通分量重最大度的点，先压入matchingorder，然后标志位设为true
    visited[selected_vertex] = true;

    // Order vertices.
    std::vector<uint32_t> tie_vertices; // 这两个辅助数组是用来存储具有相同选择度的顶点
    std::vector<uint32_t> temp;

    // 遍历该连通分量中 剩下的n-1个点，总是选择 最多（已选择）前驱结点的点 作为下一个匹配顺序点
    for (uint32_t i = 1; i < n; ++i) {
        // Select the vertices with the maximum number of backward neighbors.
        selected_vertex_selectivity = 0;
        for (uint32_t u = 0; u < n; ++u) {
            if (!visited[u]) {
                // Compute the number of backward neighbors of u.
                uint32_t u_selectivity = 0;
                for (auto uu : matching_order) {
                    if (graph->checkEdgeExistence(u, uu)) {
                        u_selectivity += 1;
                    }
                }

                // Update the vertices under consideration.
                if (u_selectivity > selected_vertex_selectivity) {
                    selected_vertex_selectivity = u_selectivity;
                    tie_vertices.clear();
                    tie_vertices.push_back(u);
                } else if (u_selectivity == selected_vertex_selectivity) {
                    tie_vertices.push_back(u);
                }
            }
        }

        // 当没有唯一选择度最高的顶点时的情况。具体来说，它处理了当存在多个顶点具有相同数量的前驱邻居时，如何选择下一个顶点加入到匹配顺序中。
        if (tie_vertices.size() != 1) {
            temp.swap(tie_vertices);
            tie_vertices.clear();

            uint32_t count = 0;
            std::vector<uint32_t> u_fn; 
            for (auto u : temp) { // 遍历这些具有 相同数量的前驱邻居 的顶点
                // Compute the number of vertices in the matching order that has at least one vertex
                // not in the matching order && connected with u.

                // Get the neighbors of u that are not in the matching order.
                uint32_t un_count;
                auto un = graph->getVertexNeighbors(u, un_count); // 获取与当前u的邻居且这个邻居不在matching_order中
                for (uint32_t j = 0; j < un_count; ++j) {
                    if (!visited[un[j]]) {
                        u_fn.push_back(un[j]);
                    }
                }

                // Compute the valid number of vertices.
                // 计算有效顶点的数量
                uint32_t cur_count = 0;
                for (auto uu : matching_order) { // 遍历匹配顺序中的点，并找他们的邻居
                    uint32_t uun_count;
                    auto uun = graph->getVertexNeighbors(uu, uun_count);
                    uint32_t common_neighbor_count = 0;
                    // 第一次计算交集，这里是复用了交集代码
                    ComputeSetIntersection::ComputeCandidates(uun, uun_count, u_fn.data(), (uint32_t)u_fn.size(), common_neighbor_count);
                    if (common_neighbor_count > 0) { // 如果共同邻居非空，则计算一个有效点
                        cur_count += 1;
                    }
                }

                u_fn.clear();

                // Update the vertices under consideration.
                if (cur_count > count) {
                    count = cur_count;
                    tie_vertices.clear();
                    tie_vertices.push_back(u);
                }
                else if (cur_count == count){
                    tie_vertices.push_back(u);
                }
            }
        }
        // 如果上边仍然一致，
        if (tie_vertices.size() != 1) {
            temp.swap(tie_vertices);
            tie_vertices.clear();

            uint32_t count = 0;
            std::vector<uint32_t> u_fn;
            for (auto u : temp) {
                // Compute the number of vertices not in the matching order && not the neighbor of vertices in the
                // matching order, but is connected with u.

                // Get the neighbors of u that are not in the matching order.
                uint32_t un_count;
                auto un = graph->getVertexNeighbors(u, un_count);
                for (uint32_t j = 0; j < un_count; ++j) { // 
                    if (!visited[un[j]]) {
                        u_fn.push_back(un[j]);
                    }
                }
                // 计算有效点的数量
                // Compute the valid number of vertices.
                uint32_t cur_count = 0;
                for (auto uu : u_fn) {
                    bool valid = true;
                    
                    for (auto uuu : matching_order) { // 遍历匹配顺序中点，只要u的邻居uu和匹配顺序中点有一个 不存在边，就立马退出
                        if (graph->checkEdgeExistence(uu, uuu)) { 
                            valid = false;
                            break;
                        }
                    }

                    if (valid) {
                        cur_count += 1;
                    }
                }

                u_fn.clear();

                // Update the vertices under consideration.
                if (cur_count > count) {
                    count = cur_count;
                    tie_vertices.clear();
                    tie_vertices.push_back(u);
                }
                else if (cur_count == count){
                    tie_vertices.push_back(u);
                }
            }
        }

        // 实在没办法了，就把第一个点假如匹配顺序吧
        matching_order.push_back(tie_vertices[0]);
        visited[tie_vertices[0]] = true;
        tie_vertices.clear();
        temp.clear();
    }
}

// ***函数：开始函数，根据查询图 创建匹配顺序（检测自同构信息 + 创建匹配顺序并 对每条边完成对匹配顺序的映射）
void OrderManager::create_orders(const Graph *query_graph) {
    // 1、检测自同构的边
    detect_automorphism_edges(query_graph); 

    // 2. Create orders and map each edge to them.
    uint32_t count = 0; // 标记 创建的顺序信息的数量
    // 遍历每一个自同构边的集合
    for (auto& automorphism : automorphism_edges_) {
        Edge edge = automorphism.front(); // 当前自同构边集合中的第一对边
        orders_.emplace_back(OrdersPerEdge()); // 存储与当前边 相关的顺序信息

        // 为简化后的查询图创建索引顺序和匹配顺序，并将结果存储在orders向量的最后一个元素
        create_indexing_order_for_reduced_graph(query_graph, edge, orders_.back());
        create_matching_order_for_reduced_graph(query_graph, edge, orders_.back());

        // 定义标签组 = 自同构第一条边，构建 标签组--边集 的映射
        LabelTriple label_triple = { query_graph->getVertexLabel(edge.first),
                                     query_graph->getEdgeLabelByVertex(edge.first, edge.second), query_graph->getVertexLabel(edge.second)};
        {
            auto it = label_edge_mapping_.find(label_triple); // 标签映射到 边集合
            if (it == label_edge_mapping_.end()) { // 如果不存在，则创建一个新的条目
                auto temp_it = label_edge_mapping_.emplace(label_triple, std::vector<Edge>());
                it = temp_it.first;
            }

            for (auto e: automorphism) { // 遍历当前自同构边集合的所有边，
                edge_orders_mapping_.insert({e, orders_.size() - 1});// 一个自同构里边的边 映射到 一个ID 上
                it->second.push_back(e); // 
            }
        }
        {
            auto it = label_automorphism_mapping_.find(label_triple);
            if (it == label_automorphism_mapping_.end()) {
                auto temp_it = label_automorphism_mapping_.emplace(label_triple, std::vector<uint32_t>());
                it = temp_it.first;
            }
            it->second.push_back(orders_.size() - 1);
            count += 1;
        }
    }

//    printf("\n%zu\n", label_automorphism_mapping_.size());
//    printf("\n%zu\n", label_edge_mapping_.size());
//    printf("\n%u\n", count);
//    count = 0;
//    for (auto& x : label_edge_mapping_) {
//        count += x.second.size();
//    }
//    printf("\n%u\n", count);
//    count = 0;
//    for (auto& x : label_automorphism_mapping_) {
//        count += x.second.size();
//    }
//    printf("\n%u\n", count);
//
//    for (auto& kv : label_edge_mapping_) {
//        auto &test = kv.second;
//        printf("%u, %u, %u: ", kv.first.src_label_, kv.first.edge_label_, kv.first.dst_label_);
//        for (auto &x: test) {
//            printf("(%u, %u)\n", x.first, x.second);
//        }
//        printf("\n");
//    }
//
//    printf("%u, %u\n", query_graph->getEdgeLabelByVertex(4, 5), query_graph->getEdgeLabelByVertex(5, 4));
//    exit(-1);
}

void OrderManager::print_info() {
    /*printf("Order Info:\n");
    printf("The number of automorphisms: %zu\n", automorphisms_.size());

    for (auto& automorphism : automorphisms_) {
        for (auto u : automorphism) {
            printf("%u ", u);
        }
        printf("\n");
    }

    printf("-----\n");
    printf("The number of disjoint set: %zu\n", automorphism_edges_.size());
    for (uint32_t i = 0; i < automorphism_edges_.size(); ++i) {
        printf("Edge set: ");
        for (auto e : automorphism_edges_[i]) {
            printf("(%u, %u), ", e.first, e.second);
        }
        printf("\n");
        printf("Index order: ");
        for (auto u : orders_[i].indexing_order_) {
            printf("%u ", u);
        }
        printf("\n");

        printf("Index order bn offset: ");
        for (auto u : orders_[i].indexing_order_bn_offset_) {
            printf("%u ", u);
        }
        printf("\n");

        printf("Index order bn: ");
        for (auto u : orders_[i].indexing_order_bn_) {
            printf("%u ", u);
        }
        printf("\n");

        printf("Matching order: ");
        for (auto u : orders_[i].matching_order_) {
            printf("%u ", u);
        }
        printf("\n");

        printf("Matching order bn offset: ");
        for (auto u : orders_[i].matching_order_bn_offset_) {
            printf("%u ", u);
        }
        printf("\n");

        printf("Matching order bn: ");
        for (auto u : orders_[i].matching_order_bn_) {
            printf("%u ", u);
        }
        printf("\n");
        printf("-----\n");
    }*/
}
