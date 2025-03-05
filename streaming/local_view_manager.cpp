#include <chrono>
#include <cmath>
#include "local_view_manager.h"
#include "streaming_config.h"
#include "computesetintersection.h"
void LocalViewManager::initialize(const Graph *query_graph, const Graph *data_graph) {
    ui n = query_graph->getVerticesCount();
    candidates_store_.resize(n);
    for (auto& candidates_set : candidates_store_) { // 为每个查询顶点的候选集预分配1024个元素的空间
        candidates_set.reserve(1024);
    }
    
    // ...............................................................................................
    views_.resize(query_graph->getEdgesCount() * 2); // 分配两倍查询边的空间，正向和反向
    buffer_pool_.reserve(1024 * 1024); // 分配1M的空间，存储 构建视图过程中 产生的中间数据
    edge_view_mapping_.reserve(256);
    flag_array_.resize(data_graph->getVerticesCount(), 0); // 记录数据图中的顶点状态，大小与数据图中的顶点数量相同
    updated_flag_.resize(data_graph->getVerticesCount());  // 
    si_buffer_.resize(data_graph->getVerticesCount());  // 与 交集运算相关
    candidate_set_size_.resize(n);
    candidate_set_pointer_.resize(n);

    u_back.resize(n);
    u_pruning.resize(n);
    u_pruning_class.resize(n);
    leaf_ID.resize(n);
}

void LocalViewManager::release() {
    candidates_store_.clear();
    views_.clear();
    buffer_pool_.clear();
    edge_view_mapping_.clear();
    flag_array_.clear();
    updated_flag_.clear();
    candidate_set_size_.clear();
    candidate_set_pointer_.clear();
    si_buffer_.clear();

    u_v_ID_map.clear();
    u_pruning.clear();
    u_pruning_class.clear();
    u_back.clear();
    leaf_ID.clear();
}

void LocalViewManager::destroy_view() {
    for (auto& view : views_) {
        view.clear();
    }

    for (auto& candidate_set : candidates_store_) {
        candidate_set.clear();
    }

    buffer_pool_.clear();
    edge_view_mapping_.clear();

    u_v_ID_map.clear();
}

LocalEdgeView *LocalViewManager::get_view(uint32_t view_id) { 
    assert(view_id < views_.size());
    return &views_[view_id];
} 
LocalEdgeView *LocalViewManager::get_view(Edge query_edge) {
    auto it = edge_view_mapping_.find(query_edge);
    if (it != edge_view_mapping_.end()) {
        return &views_[it->second];
    }
    return nullptr;
}
uint32_t LocalViewManager::get_view_id(Edge query_edge) {
    auto it = edge_view_mapping_.find(query_edge);
    if (it != edge_view_mapping_.end()) {
        return it->second;
    }
    // Local view does not exist.
    return edge_view_mapping_.size();
}
uint32_t *LocalViewManager::get_candidate_set(uint32_t u) {
    return candidates_store_[u].data();
}
uint32_t LocalViewManager::get_candidate_set_size(uint32_t u) {
    return candidates_store_[u].size();
}

// ***函数： 创建全局视图，返回的局部视图创建是否成功
bool
LocalViewManager::create_view(const Graph *query_graph, OrdersPerEdge &orders, GlobalViewManager &gvm, Edge data_edge) {

    build_visited_neighbor_count_  = 0; // 已经访问了多少邻居
    generate_visited_neighbor_count_ = 0; // 记录生成候选集时，访问的邻居数量
    first_vertex_neighbor_ = 0; // 记录生成局部候选集中，第一个顶点的邻居数量
    bool is_valid = optimized_generate_local_candidates_v2(query_graph, orders, gvm, data_edge);
    if (!is_valid)  // 
        return false;
    deal_is_prunning(query_graph, orders);
    optimized_build_local_view(query_graph, orders, gvm, data_edge);
    // print_local_information(query_graph);
    return true;
}

ui LocalViewManager::get_ID(ui u, ui v){
    return u_v_ID_map[u][v];
}

bool LocalViewManager::optimized_generate_local_candidates(const Graph *query_graph, OrdersPerEdge &orders,
                                                           GlobalViewManager &gvm, Edge exclude_data_edge) {
    // Generate the candidate set for each vertex along the indexing order.
    auto& indexing_order = orders.indexing_order_;
    auto& indexing_order_bn = orders.indexing_order_bn_;
    auto& indexing_order_bn_offset = orders.indexing_order_bn_offset_;

    uint32_t u0 = indexing_order[0];
    uint32_t u1 = indexing_order[1];
    uint32_t v0 = exclude_data_edge.first;
    uint32_t v1 = exclude_data_edge.second;

    if (!gvm.nlf_check(u0, v0)
        || !gvm.nlf_check(u1, v1)) {
        return false;
    }

    if (query_graph->getVertexDegree(u0) > 1) {
        uint32_t global_encoded_v0 = gvm.get_encoded_id(u0, v0);
        candidates_store_[u0].push_back(global_encoded_v0);
    }
    else {
        candidates_store_[u0].push_back(v0);
    }

    if (query_graph->getVertexDegree(u1) > 1) {
        uint32_t global_encoded_v1 = gvm.get_encoded_id(u1, v1);
        candidates_store_[u1].push_back(global_encoded_v1);
    }
    else {
        candidates_store_[u1].push_back(v1);
    }

    for (uint32_t i = 2; i < indexing_order.size(); ++i) {
        uint32_t u = indexing_order[i];

        // Skip leaf node.
        if (query_graph->getVertexDegree(u) == 1)
            continue;

        uint32_t begin = indexing_order_bn_offset[i];
        uint32_t end = indexing_order_bn_offset[i + 1];
        updated_count_ = 0;

        uint32_t flag_value = 0;
        for (uint32_t j = begin; j < end; ++j) {
            uint32_t uu = indexing_order_bn[j];
            auto gv = gvm.get_nlf_view({uu, u});
            auto reverse_gv = gvm.get_nlf_view({u, uu});

            bool flag = false;
            if (end - begin >= 2) {
                flag = (indexing_order_bn[begin] == u0 && indexing_order_bn[begin + 1] == u1);
            }
            for (auto vv : candidates_store_[uu]) {
                uint32_t vv_nbrs_count;
                auto vv_nbrs = gv->get_neighbor(vv, vv_nbrs_count);
                generate_visited_neighbor_count_ += vv_nbrs_count;
                if ((uu == indexing_order[0] || uu == indexing_order[1]) && flag)
                    first_vertex_neighbor_ += vv_nbrs_count;

                // If it is the first bn or the cost of binary search is lower than that of scan, then use the first approach.
                if ((j == begin) || vv_nbrs_count < 1024 || vv_nbrs_count < updated_count_ * 32) {
                    for (uint32_t k = 0; k < vv_nbrs_count; ++k) {
                        uint32_t v = vv_nbrs[k];

                        if (flag_array_[v] == flag_value) {
                            flag_array_[v] += 1;
                            if (flag_value == 0) {
                                updated_flag_[updated_count_++] = v;
                            }
                        }
                    }
                }
                else {
                    for (uint32_t k = 0; k < updated_count_; ++k) {
                        uint32_t v = updated_flag_[k];
                        if (flag_array_[v] == flag_value) {
                            uint32_t v_nbrs_count;
                            auto v_nbrs = reverse_gv->get_neighbor(v, v_nbrs_count);
                            if (vv_nbrs_count < v_nbrs_count) {
                                auto it = std::lower_bound(vv_nbrs, vv_nbrs + vv_nbrs_count, v);
                                if (it != vv_nbrs + vv_nbrs_count && *it == v) {
                                    flag_array_[v] += 1;
                                }
                            }
                            else {
                                auto it = std::lower_bound(v_nbrs, v_nbrs + v_nbrs_count, vv);
                                if (it != v_nbrs + v_nbrs_count && *it == vv) {
                                    flag_array_[v] += 1;
                                }
                            }
                        }
                    }
                }
            }

            flag_value += 1;
        }
        uint32_t local_encoded_v0 = std::numeric_limits<uint32_t>::max();
        if (gvm.candidate_check(u, v0)) {
            local_encoded_v0 = gvm.get_encoded_id(u, v0);
        }

        uint32_t local_encoded_v1 = std::numeric_limits<uint32_t>::max();
        if (gvm.candidate_check(u, v1)) {
            local_encoded_v1 = gvm.get_encoded_id(u, v1);
        }

        for (uint32_t j = 0; j < updated_count_; ++j) {
            uint32_t v = updated_flag_[j];
            if (flag_array_[v] == flag_value && v != local_encoded_v0 && v != local_encoded_v1) {
                candidates_store_[u].push_back(v);
            }

            flag_array_[v] = 0;
        }

        if (candidates_store_[u].empty())
            return false;
    }

    for (auto& candidate_set : candidates_store_) {
        std::sort(candidate_set.begin(), candidate_set.end());
    }

    return true;
}
void LocalViewManager::deal_is_prunning(const Graph *query_graph, OrdersPerEdge &orders){
    auto& matching_order = orders.matching_order_;
    auto& index_order = orders.indexing_order_;
    ui u0 = index_order[0];
    ui u1 = index_order[1];
    LabelID label0 = query_graph->getVertexLabel(u0);
    LabelID label1 = query_graph->getVertexLabel(u1);
    std::fill(u_pruning.begin(),u_pruning.end(),false);
    for(ui i = 0; i < matching_order.size(); i++){
        ui u = matching_order[i];
        LabelID label = query_graph->getVertexLabel(u);
        ui count =  query_graph->getLabelsFrequency(label);
        if(label == label0 || label == label1){
            u_pruning[u] = count == 2 ? true : false;
            continue;
        }
        u_pruning[u] = count == 1 ? true : false; 
    }
}

void LocalViewManager::optimized_build_local_view(const Graph *query_graph, OrdersPerEdge &orders, GlobalViewManager &gvm, Edge data_edge) {
    auto& matching_order = orders.matching_order_;
    auto& matching_order_bn = orders.matching_order_bn_;
    auto& matching_order_bn_offset = orders.matching_order_bn_offset_;
    auto& matching_order_view_mappings = orders.matching_order_view_mappings_; // very important
    auto& matching_order_relation_edge_type = orders.matching_order_edge_type_;
    
    uint32_t local_view_id = 0;
      
    // 对u_v_vector2初始化
    std::vector<std::vector<std::vector<ui>>> u_v_vector2;
    ui N = query_graph->getVerticesCount();
    u_v_vector2.resize(N);
    for(ui u = 0; u < N; u++){
        if(candidate_set_size_[u] == 0)
            continue;
        u_v_vector2[u].resize(candidate_set_size_[u]);
    }
    std::fill(u_back.begin(),u_back.end(),false);
    // 遍历匹配顺序中的查询点
    for (uint32_t i = 0; i < matching_order.size(); ++i) {
        uint32_t u = matching_order[i];
        updated_count_ = 0; // 记录的是当前查询点u对应的候选点 有几个更新了
        // 这里是遍历当前查询点u的 前驱邻居uu
        for (uint32_t j = matching_order_bn_offset[i]; j < matching_order_bn_offset[i + 1]; ++j) {
            uint32_t uu = matching_order_bn[j];
            u_back[uu] = true;
            RelationEdgeType type = matching_order_relation_edge_type[j];
            Edge query_edge = {uu, u};
            Edge reverse_query_edge = {u, uu};
            // 如果当前查询点的度 为1，为叶子结点，直接获取它的视图ID
            if (query_graph->getVertexDegree(u) == 1) {
                //强行处理叶子节点
                auto gv = gvm.get_view(query_edge); // 获取uu-u的视图
                // printf("发现叶子节点u%u和它前驱结点u%u, 前驱结点u%u的候选点规模为%u", u, uu, uu, candidates_store_[uu].size());
                for(ui k = 0; k < candidates_store_[uu].size(); k++){

                    ui vv = candidates_store_[uu][k];
                    ui true_vv = gvm.get_id(uu,vv);
                    // printf("v%u: ",true_vv);
                    ui lc_count = gv->get_neighbor_num(true_vv);
                    if(u_pruning[uu]){
                        u_v_vector2[uu][k].push_back(lc_count);
                    }
                }
                matching_order_view_mappings[j]= gvm.get_view_id(query_edge);
                continue;
            }
            // 非叶子：更新匹配顺序-view-ID映射, 边视图映射
            matching_order_view_mappings[j] = local_view_id;
            edge_view_mapping_[query_edge] = local_view_id;
            
            // 先取local_view_id的值，再 + 1
            auto &lv = views_[local_view_id++]; // 获取其view,下面做填充和修改
            // printf("LocalView[%u]:现在开始处理u%u和它的前驱u%u,即(u%u-u%u)视图  \n",local_view_id-1,u,uu,uu,u);
            if (type == RelationEdgeType::REGULAR) {
                // Set the flag array.
                if (updated_count_ == 0) { // 如果还没有更新的候选点
                    for (uint32_t k = 0; k < candidates_store_[u].size(); ++k) {
                        // modified by gz********************************************************************
                        flag_array_[candidates_store_[u][k]] = k + 1;
                        // flag_array_[candidates_store_[u][k]] = candidates_store_[u][k] + 1;
                        updated_flag_[updated_count_++] = candidates_store_[u][k];
                    }
                }
                // 生成局部视图
                // Generate the local view.
                auto gv = gvm.get_nlf_view(query_edge); // 获取正向和反向的全局NLF视图
                auto reverse_gv = gvm.get_nlf_view(reverse_query_edge);
                // 遍历前驱结点uu的候选集
                
                // for (auto vv : candidates_store_[uu]) {
                for (ui l = 0; l < candidate_set_size_[uu]; l++){
                    ui vv = candidate_set_pointer_[uu][l];
                    // ui cnt = 0;
                    uint32_t begin_pos = buffer_pool_.size(); // 开始下标为buffer_pool的大小

                    uint32_t vv_nbrs_count;
                    auto vv_nbrs = gv->get_neighbor(vv, vv_nbrs_count);
                    build_visited_neighbor_count_ += vv_nbrs_count;
                    // 如果vv的邻接点规模不是太大
                    // printf("(u%u--u%u):",uu,u);
                    if (vv_nbrs_count < 1024 || vv_nbrs_count < updated_count_ * 32) {
                        for (uint32_t k = 0; k < vv_nbrs_count; ++k) {
                            uint32_t v = vv_nbrs[k];
                            // printf("k=%u: flag_array_[%u] = %u\n",k, v, flag_array_[v]);
                            if (flag_array_[v] > 0) {
                                buffer_pool_.push_back(flag_array_[v] - 1);
                                // cnt++;
                                if(u_pruning[uu]){
                                    // printf("u%u的候选点(%u--v%u)加入(%u--v%u)\n",uu, vv, gvm.get_id(uu,vv), flag_array_[v]-1, gvm.get_id(u,flag_array_[v]-1));
                                    // printf("u%u的候选点(%u--v%u)加入(%u--v%u)\n",uu, vv, gvm.get_id(uu,vv), candidates_store_[u][k], gvm.get_id(u,candidates_store_[u][k]));
                                    //printf("u%u的候选点(%u--v%u)加入(%u--v%u)\n",uu, vv, gvm.get_id(uu,vv),v , gvm.get_id(u,v));
                                    u_v_vector2[uu][l].push_back(v); // 这里想想办法
                                }
                            }
                        }
                    }
                    else {
 //                       printf("B %u, %u, %d, %d\n", candidates_store_[u].size(), vv_nbrs_count, updated_count_);
                        for (auto v : candidates_store_[u]) {
                            uint32_t v_nbrs_count;
                            auto v_nbrs = reverse_gv->get_neighbor(v, v_nbrs_count);
                            if (vv_nbrs_count < v_nbrs_count) {
                                auto it = std::lower_bound(vv_nbrs, vv_nbrs + vv_nbrs_count, v);
                                if (it != vv_nbrs + vv_nbrs_count && *it == v) {
                                    buffer_pool_.push_back(flag_array_[v] - 1);
                                    // cnt++;
                                    if(u_pruning[uu])
                                        u_v_vector2[uu][l].push_back(v);
                                }
                            }
                            else {
                                auto it = std::lower_bound(v_nbrs, v_nbrs + v_nbrs_count, vv);
                                if (it != v_nbrs + v_nbrs_count && *it == vv) {
                                    buffer_pool_.push_back(flag_array_[v] - 1);
                                    // cnt++;
                                    if(u_pruning[uu])
                                        u_v_vector2[uu][l].push_back(v);
                                }
                            }
                        }
                    }
                    //if(u_pruning[uu]){   
                        //u_v_vector2[uu][l].insert(u_v_vector2[uu][l].end(),buffer_pool_.begin() + begin_pos, buffer_pool_.end());
                        //u_v_vector2[uu][l].push_back(buffer_pool_[begin_pos]);
                        //u_v_vector2[uu][l].push_back(cnt);
                    //}
                    // 
                    // printf("LocalView[%u]信息",local_view_id);
                    // modified by gz*******************************************************************
                    lv.trie_.emplace_back(begin_pos, buffer_pool_.size());
                    lv.cardinality_ += buffer_pool_.size() - begin_pos;
                }
                
            }

        }
        
        for (uint32_t j = 0; j < updated_count_; ++j) {
            uint32_t v = updated_flag_[j];
            flag_array_[v] = 0;
        }
    }
    for (uint32_t i = 0; i < local_view_id; ++i) {
        views_[i].data_ = buffer_pool_.data();
    }
    ui ID = 1;
    #pragma omp parallel for
    for(ui i = 0; i < matching_order.size(); i++){
        ui u = matching_order[i];
        if(!u_pruning[u]){
            u_pruning_class[u] = 0;       
            continue; 
        } 
        if(query_graph->getVertexDegree(u) == 1){ // 如果度为1，则统一设置一个ID即可
            u_pruning_class[u] = 1;
            leaf_ID[u] = ID;
            ID++;
            continue;
        }
        if(!u_back[u]){ // 如果没有以它为前驱的结点
            u_pruning_class[u] = 1;
            leaf_ID[u] = ID;
            ID++;
            continue;
        }
        u_pruning_class[u] = 2;
        VectorMap vector_map;
        for(ui k = 0; k < candidate_set_size_[u]; k++){
            ui encode_v = candidate_set_pointer_[u][k];
            std::vector<ui> &tmp = u_v_vector2[u][k];
            vector_map[tmp].push_back(encode_v);
        }

        for(auto it = vector_map.begin(); it != vector_map.end(); it++){
            #pragma omp parallel for
            for(ui i = 0; i < it->second.size(); i++){
                ui true_v = gvm.get_id(u, it->second[i]);
                u_v_ID_map[u][true_v] = ID; 
            }
            ID++;
        }
    }
    /*for(ui i = 0; i < matching_order.size(); i++){
        ui u = matching_order[i];
        if(u_pruning[u]){
            printf("u%u的ID信息: ",u);
            for(ui j = 0; j < candidates_store_[u].size(); j++){
                ui v = candidates_store_[u][j];
                ui true_v = gvm.get_id(u,v);
                printf("v%u ID:%u      ", true_v, get_ID(u,true_v));
            }
        }
        else{
            printf("u%u层不剪枝",u);
        }
        printf("\n");
    }*/
    for (uint32_t u = 0; u < query_graph->getVerticesCount(); ++u) {
        if (query_graph->getVertexDegree(u) > 1) {
            for (uint32_t i = 0; i < candidates_store_[u].size(); ++i) {
                candidates_store_[u][i] = gvm.get_id(u, candidates_store_[u][i]);
            }
        }
    }
    
}

bool LocalViewManager::prune_local_candidates(GlobalViewManager &gvm, uint32_t u, std::vector<uint32_t> &bn) {
    uint32_t flag_value = 1; //标志位，用于在flag_array中标记有效的候选顶点，不会和v0 v1 产生匹配冲突的标志
    uint32_t current_valid_count = updated_count_;// update_count是更新的候选点的数量
    // 遍历u的各个前驱结点uu
    for (auto uu : bn) {
        // 获取uu-u和u-uu的NLF视图（正向 + 反向）
        auto gv = gvm.get_nlf_view({uu, u});
        auto reverse_gv = gvm.get_nlf_view({u, uu});
        // 
        uint32_t target_valid_count = current_valid_count;
        current_valid_count = 0; // 记录 u结点 当前的有效候选数量 
        // 遍历uu的候选集vv，并计算vv的邻居
        for (uint32_t j = 0; j < candidate_set_size_[uu]; ++j) {
            uint32_t vv = candidate_set_pointer_[uu][j];
            uint32_t vv_nbrs_count;
            auto vv_nbrs = gv->get_neighbor(vv, vv_nbrs_count);

            // If the cost of binary search is lower than that of scan, then use the first approach.
            // 如果vv的邻居规模较小，则基于扫描的方法即可
            if (vv_nbrs_count < 1024 || vv_nbrs_count < updated_count_ * 32) {
                for (uint32_t k = 0; k < vv_nbrs_count; ++k) { // 遍历vv的邻居v,
                    uint32_t v = vv_nbrs[k];
                    // 如果v的标志位等于当前的标志位，则增加其标志位并增加有效候选顶点的数量
                    if (flag_array_[v] == flag_value) {
                        flag_array_[v] += 1;
                        current_valid_count += 1;
                    }
                }
            }
            // 基于二分
            else {
                // 遍历之前更新的候选顶点，
                for (uint32_t k = 0; k < updated_count_; ++k) {
                    uint32_t v = updated_flag_[k];
                    // 如果标志v的标志位 等于 当前的标志。
                    if (flag_array_[v] == flag_value) {
                        uint32_t v_nbrs_count;
                        auto v_nbrs = reverse_gv->get_neighbor(v, v_nbrs_count);
                        if (vv_nbrs_count < v_nbrs_count) { // 如果vv的邻居数量 小于v的邻居数量
                            auto it = std::lower_bound(vv_nbrs, vv_nbrs + vv_nbrs_count, v);
                            if (it != vv_nbrs + vv_nbrs_count && *it == v) { //如果找到v,则更新标志位
                                flag_array_[v] += 1;
                                current_valid_count += 1;
                            }
                        }
                        else {
                            auto it = std::lower_bound(v_nbrs, v_nbrs + v_nbrs_count, vv);
                            if (it != v_nbrs + v_nbrs_count && *it == vv) {
                                flag_array_[v] += 1;
                                current_valid_count += 1;
                            }
                        }
                    }
                }
            }
            // 有效候选顶点的数量达到了目标数量，则停止迭代
            if (current_valid_count == target_valid_count)
                break;
        }

        flag_value += 1;// 增加标志值，为迭代下一个前驱结点做准备
    }
    
    bool push = candidates_store_[u].empty(); // 检查u的候选点是否 为空
    // 下面更新候选集
    uint32_t local_pos = 0;
    for (uint32_t j = 0; j < updated_count_; ++j) { // 遍历之前更新的候选点数量
        uint32_t v = updated_flag_[j]; // 获取之前更新的一个候选顶点v
        if (flag_array_[v] == flag_value) {// **********成功用flag_array[]数组，达到了削减 候选集数量的 目的
            if (push) // 如果为空
                candidates_store_[u].push_back(v);
            else // 如果不为空，则添加即可
                candidates_store_[u][local_pos++] = v;
        }

        flag_array_[v] = 0;  // 标志位清零，为下一个 u 做准备
    }

    if (!push)  // 如果非空，重新改变candidates_store_[]的大小
        candidates_store_[u].resize(local_pos);

    candidate_set_pointer_[u] = candidates_store_[u].data();
    candidate_set_size_[u] = candidates_store_[u].size();

    if (candidate_set_size_[u] == 0)
        return false;

    return true;
}
bool LocalViewManager::optimized_generate_local_candidates_v2(const Graph *query_graph, OrdersPerEdge &orders,
                                                              GlobalViewManager &gvm, Edge exclude_data_edge) {
    // Generate the candidate set for each vertex along the indexing order.
    // 沿着索引顺序，为每一个点生成候选集。记录该边记录的索引顺序 + 前驱邻居信息 + 前驱邻居列表的迁移量
    auto& indexing_order = orders.indexing_order_;
    auto& indexing_order_bn = orders.indexing_order_bn_;
    auto& indexing_order_bn_offset = orders.indexing_order_bn_offset_;
    // 记录索引顺序的前两个点，和自同构边集对应第一条边的两个点，唯一且单一匹配
    uint32_t u0 = indexing_order[0];
    uint32_t u1 = indexing_order[1];
    uint32_t v0 = exclude_data_edge.first;      // 数据更新边的第一点
    uint32_t v1 = exclude_data_edge.second;     // 数据更新边的第二点

    if (!gvm.nlf_check(u0, v0) // NLF检查不通过，直接返回false
        || !gvm.nlf_check(u1, v1)) {
        return false;
    }
    // 如果u0的度大于1, 不是叶子顶点，取候选集u0v0的编码ID，存入候选集。如果是叶子结点，则直接添加v0
    if (query_graph->getVertexDegree(u0) > 1) { 
        uint32_t global_encoded_v0 = gvm.get_encoded_id(u0, v0);
        candidates_store_[u0].push_back(global_encoded_v0);
    }
    else {
        candidates_store_[u0].push_back(v0);
    }
    // 
    if (query_graph->getVertexDegree(u1) > 1) {
        uint32_t global_encoded_v1 = gvm.get_encoded_id(u1, v1);
        candidates_store_[u1].push_back(global_encoded_v1);
    }
    else {
        candidates_store_[u1].push_back(v1);
    }
    // 这里的u0 u1的候选集只能是v0,v1，候选集的规模只能是1,1（因为是CSM的性质）
    candidate_set_size_[u0] = 1;
    candidate_set_pointer_[u0] = candidates_store_[u0].data();
    candidate_set_size_[u1] = 1;
    candidate_set_pointer_[u1] = candidates_store_[u1].data();

    // First, compute the triangle.
    // 首先计算三角形，和前面的u0，u1构成一个三角形。这里的triangle_end_index = both_two_vertex + 2
    for (uint32_t i = 2; i < orders.triangle_end_index_; ++i) {
        // Skip the first two vertex.
        uint32_t u = indexing_order[i]; // 取出第一个与u0u1构成三角形的点

        // 获取与{u0,u}边相对应的NLF视图
        // 这里的candidates_set_pointer[u0][0]：u0对应候选集的第一个候选点，也是唯一一个点
        auto gv0 = gvm.get_nlf_view({u0, u});
        uint32_t v0_nbr_count;
        auto v0_nbr = gv0->get_neighbor(candidate_set_pointer_[u0][0], v0_nbr_count);
        
        auto gv1 = gvm.get_nlf_view({u1, u});
        uint32_t v1_nbr_count;
        auto v1_nbr = gv1->get_neighbor(candidate_set_pointer_[u1][0], v1_nbr_count);

        // 这里计算一次v0与v1的一次邻居,匹配一次三角环结构
        uint32_t lc = 0;
        ComputeSetIntersection::ComputeCandidates(v0_nbr, v0_nbr_count, v1_nbr, v1_nbr_count, updated_flag_.data(), lc);
        // 找不到共同邻居，无法匹配三角环，直接返回失败
        if (lc == 0)
            return false;

        // 将找到的交集放入 顶点u的候选集中，u的邻居就应该很多了，因为不限制数据图的第三个点
        candidates_store_[u].insert(candidates_store_[u].end(), updated_flag_.begin(), updated_flag_.begin() + lc);
        candidate_set_pointer_[u] = candidates_store_[u].data();
        candidate_set_size_[u] = candidates_store_[u].size();
    }
    //printf("\n不用测试，这里真的很简单。只需要注意这里的candidates_store_[][]存放的是编码ID↑\n");
    //printf(%)
    
    // Second, set the candidate for vertex adjacent to update.
    // 第二步，为与{u0，u1}一个顶点的点设置候选集
    for (uint32_t i = orders.triangle_end_index_; i < orders.adjacent_end_index_; ++i) {
        uint32_t u = indexing_order[i];

        // 跳过叶子结点
        // Skip leaf node.
        if (query_graph->getVertexDegree(u) == 1)
            continue;
        // 找到u的前驱结点uu，这里的uu一定是已经被访问过的了
        uint32_t uu = indexing_order_bn[indexing_order_bn_offset[i]];
        auto gv = gvm.get_nlf_view({uu, u}); // 获取uu--u的NLF视图
        // 这里是找到vv的邻居，candidates_set_pointer[uu][0]是指：uu候选集的第一个候选点，它也only只有一个
        candidate_set_pointer_[u] = gv->get_neighbor(candidate_set_pointer_[uu][0], candidate_set_size_[u]);

        if (candidate_set_size_[u] < 1024) {  // 说明候选集不是很大？
            uint32_t local_encoded_v0 = std::numeric_limits<uint32_t>::max(); // 假定为最大值
            if (gvm.candidate_check(u, v0)) { // 检查v0或v1是否为u的候选集
                local_encoded_v0 = gvm.get_encoded_id(u, v0);
            }

            uint32_t local_encoded_v1 = std::numeric_limits<uint32_t>::max();
            if (gvm.candidate_check(u, v1)) {
                local_encoded_v1 = gvm.get_encoded_id(u, v1);
            }
            // 循环遍历u的候选集，候选点用v表示，填充候选集
            for (uint32_t j = 0; j < candidate_set_size_[u]; ++j) {
                uint32_t v = candidate_set_pointer_[u][j];
                if (v != local_encoded_v0 && v != local_encoded_v1) { // 粗略理解：如果u的候选点不等于v0 v1本身 
                    candidates_store_[u].push_back(v); // 这样的话，就可以把这些候选点压入 候选集，
                }
            }

            candidate_set_pointer_[u] = candidates_store_[u].data();
            candidate_set_size_[u] = candidates_store_[u].size();
        }

        if (candidate_set_size_[u] == 0) // 在索引到u的候选点时，发现为空，已经不可能完全匹配了。
            return false;
    }

        /*
         * Regenerate indexing order for adjacent vertex based on candidate set size.
         */
        // 基于候选集的大小 为adjacent点重新生成 索引顺序

        // 创建一个新的索引顺序，从索引顺序的第三个点开始，直到相邻点的最后一个点结束。并把候选集规模较小的点排在前面
        std::vector<uint32_t> optimized_indexing_order(indexing_order.begin() + 2,
                                                       indexing_order.begin() + orders.adjacent_end_index_);
        std::sort(optimized_indexing_order.begin(), optimized_indexing_order.end(), [this](uint32_t a, uint32_t b) {
            return candidate_set_size_[a] < candidate_set_size_[b];
        });
        // 创建索引顺序的前驱邻居和偏移量
        std::vector<bool> visited(query_graph->getVerticesCount(), false);
        std::vector<uint32_t> optimized_indexing_order_bn;
        std::vector<uint32_t> optimized_indexing_order_bn_offset;

        optimized_indexing_order_bn_offset.push_back(0);
        // 遍历这些单个相邻顶点的索引顺序，目的是填写好optimized_indexing_order_bn和optimized_indexing_order_bn_offset数组
        for (auto u : optimized_indexing_order) { 
            uint32_t u_nbrs_count;
            auto u_nbrs = query_graph->getVertexNeighbors(u, u_nbrs_count); 
            for (uint32_t i = 0; i < u_nbrs_count; ++i) {
                uint32_t uu = u_nbrs[i];
                if (visited[uu]) {
                    optimized_indexing_order_bn.push_back(uu);
                }
            }
            optimized_indexing_order_bn_offset.push_back(optimized_indexing_order_bn.size());
            visited[u] = true;
        }

    // 再次遍历index_order，这里再次遍历的目的是要 更新候选点
    for (uint32_t i = 2; i < indexing_order.size(); ++i) {
        uint32_t u;
        std::vector<uint32_t> bn; // 存储u的前驱邻居
        if (i < orders.adjacent_end_index_) { // 如果还在adjacent_end_index范围之内，就使用新的优化后的索引顺序
            uint32_t index = i - 2;
            u = optimized_indexing_order[index]; 
            uint32_t begin = optimized_indexing_order_bn_offset[index];
            uint32_t end = optimized_indexing_order_bn_offset[index + 1];
            bn.insert(bn.end(), optimized_indexing_order_bn.begin() + begin, optimized_indexing_order_bn.begin() + end);
        }
        else { // 如果出了范围， 就还是使用原本的索引顺序
            u = indexing_order[i];
            uint32_t begin = indexing_order_bn_offset[i];
            uint32_t end = indexing_order_bn_offset[i + 1];
            bn.insert(bn.end(), indexing_order_bn.begin() + begin, indexing_order_bn.begin() + end);
        }
        /*printf("*u%u的前驱为: ", u);
        for(ui j = 0; j < bn.size(); j++)
            printf("u%u ",bn[j]);
        printf("\n");*/
        if (query_graph->getVertexDegree(u) == 1 || bn.empty()) // 如果当前结点u为叶子结点或者它的前驱邻居数量为空
            continue;


        uint32_t local_encoded_v0 = std::numeric_limits<uint32_t>::max();
        if (gvm.candidate_check(u, v0)) { // 这里是不是防止产生自冲突
            local_encoded_v0 = gvm.get_encoded_id(u, v0);
        }
        uint32_t local_encoded_v1 = std::numeric_limits<uint32_t>::max();
        if (gvm.candidate_check(u, v1)) {
            local_encoded_v1 = gvm.get_encoded_id(u, v1);
        }
        
        /*printf("u%u的前驱为: ", u);
        for(ui j = 0; j < bn.size(); j++)
            printf("u%u ",bn[j]);
        printf("\n");*/
        // 如果是adjacent_end_index范围内，设置set_adjacent_update_candidates_flag，否则设置set_update_candidate_flag
        if (i < orders.adjacent_end_index_) {
            set_adjacent_update_candidates_flag(gvm, u, bn, local_encoded_v0, local_encoded_v1);
        }
        else {
            set_update_candidates_flag(gvm, u, bn, local_encoded_v0, local_encoded_v1);
        }

        if (!prune_local_candidates(gvm, u, bn)) 
            return false;
    }
    // 遍历索引顺序，填充候选集，并且排序
    // printf("剪枝完");
    // print_local_information();
    for (auto u : indexing_order) {
        if (query_graph->getVertexDegree(u) > 1) {
            if (candidates_store_[u].empty()) {  // 对于|C(u)|大于1024的u，这里作填充
                for (uint32_t i = 0; i < candidate_set_size_[u]; ++i) {  // 无效的
                    candidates_store_[u].push_back(candidate_set_pointer_[u][i]);
                }
            }
            std::sort(candidates_store_[u].begin(), candidates_store_[u].end());
            candidate_set_pointer_[u] = candidates_store_[u].data();
            candidate_set_size_[u] = candidates_store_[u].size();
        }
    }
    return true;
}
// 函数：为adjacent的u查询点设置更新候选集 标志，参数（全局视图，当前查询点u，u的前驱邻居，更新边两个点的候选集编码（记录u与 前两个点 匹配冲突））
void LocalViewManager::set_adjacent_update_candidates_flag(GlobalViewManager &gvm, uint32_t u, std::vector<uint32_t> &bn,
                                                      uint32_t encoded_v0, uint32_t encoded_v1) {
    // 设置is_falg_based标志位，标志是否启动设置候选更新标志。如果u的候选集规模 小于1024 或者小于任何一个前驱结点候选集的32倍
    bool is_flag_based = false;
    if (candidate_set_size_[u] < 1024) {
        is_flag_based = true;
    }
    else {
        for (auto uu : bn) {
            if (candidate_set_size_[u] < 32 * candidate_set_size_[uu]) {
                is_flag_based = true;
                break;
            }
        }
    }
    // 如果启动了使用 基于标志的方法
    if (is_flag_based) {
        updated_count_ = 0; // 记录更新后符合条件的候选点数量
        for (uint32_t i = 0; i < candidate_set_size_[u]; ++i) {  // 遍历u的候选集
            uint32_t v = candidate_set_pointer_[u][i];
            if (v != encoded_v0 && v != encoded_v1) {
                flag_array_[v] = 1;
                updated_flag_[updated_count_++] = v;
            }
        }
    }

    // 如果不使用 基于标志的方法，即可候选集规模较大，即必须上交集来缩减候选集的规模了，重点关注一下怎么缩减规模的
    else {
        uint32_t uu = select_bn_with_minimum_degree_sum(gvm, u, bn); // 选择一个前驱节点，这个节点拥有最小的度数
        auto gv = gvm.get_nlf_view({uu, u}); // 获取uu-u的全局NLF视图

        updated_count_ = 0; // 
        // 遍历uu的所有候选顶点
        for (uint32_t i = 0; i < candidate_set_size_[uu]; ++i) {
            // 获取uu的候选顶点vv，获取vv的邻居（匹配uu）
            uint32_t vv = candidate_set_pointer_[uu][i];
            uint32_t vv_nbrs_count;
            auto vv_nbrs = gv->get_neighbor(vv, vv_nbrs_count);

            // 获取u的候选集 与 uu邻居vv集合的 交集，结果存放在si_buffed中
            uint32_t lc = 0;
            ComputeSetIntersection::ComputeCandidates(candidate_set_pointer_[u], candidate_set_size_[u],
                                                      vv_nbrs, vv_nbrs_count, si_buffer_.data(), lc); // 前驱匹配点和和u的候选点之间的交集

            // 遍历这些交集顶点
            for (uint32_t j = 0; j < lc; ++j) {
                uint32_t v = si_buffer_[j];
                // 如果flag[v]为0，则更新为1                
                if (flag_array_[v] == 0 && v != encoded_v0 && v != encoded_v1) {
                    flag_array_[v] = 1;
                    updated_flag_[updated_count_++] = v;
                }
            }
            // 如果原始的候选集大小 和 更新的候选集数量 相同，则退出，意味着所有候选点都已检查
            if (candidate_set_size_[u] == updated_count_)
                break;
        }
    }
}

// 为非 邻接的点 更新标志味
void LocalViewManager::set_update_candidates_flag(GlobalViewManager &gvm, uint32_t u, std::vector<uint32_t> &bn,
                                                  uint32_t encoded_v0, uint32_t encoded_v1) {
    // 选择具有最小度数和的 u的前驱结点
    uint32_t selected_bn = select_bn_with_minimum_degree_sum(gvm, u, bn);
    auto gv = gvm.get_nlf_view({selected_bn, u});
    updated_count_ = 0;
    for (uint32_t i = 0; i < candidate_set_size_[selected_bn]; ++i) {
        uint32_t vv = candidate_set_pointer_[selected_bn][i];
        uint32_t vv_nbrs_count;
        auto vv_nbrs = gv->get_neighbor(vv, vv_nbrs_count);

        for (uint32_t j = 0; j < vv_nbrs_count; ++j) {
            uint32_t v = vv_nbrs[j];
            if (flag_array_[v] == 0 && v != encoded_v0 && v != encoded_v1) {
                flag_array_[v] = 1;
                updated_flag_[updated_count_++] = v;
            }
        }
    }
}

uint32_t
LocalViewManager::select_bn_with_minimum_degree_sum(GlobalViewManager &gvm, uint32_t u, std::vector<uint32_t> &bn) {
    std::sort(bn.begin(), bn.end(), [this](const uint32_t a, const uint32_t b) -> bool {
        return candidate_set_size_[a] < candidate_set_size_[b];
    });

    uint64_t min_degree_sum = std::numeric_limits<uint64_t>::max();
    uint32_t selected_bn = 0;
    uint32_t selected_idx = 0;
    for (uint32_t i = 0; i < bn.size(); ++i) {
        uint32_t uu = bn[i];
        auto gv = gvm.get_nlf_view({uu, u});
        uint64_t local_degree_sum = 0;
        if (candidate_set_size_[uu] < min_degree_sum) {
            for (uint32_t j = 0; j < candidate_set_size_[uu]; ++j) {
                uint32_t vv = candidate_set_pointer_[uu][j];
                local_degree_sum += gv->get_neighbor_num(vv);
            }

            if (local_degree_sum < min_degree_sum) {
                min_degree_sum = local_degree_sum;
                selected_bn = uu;
                selected_idx = i;
            }
        }
    }

    bn.erase(bn.begin() + selected_idx);
    return selected_bn;
}



void LocalViewManager::print_local_information(const Graph* query_graph)
{
    printf("局部NLF视图信息\n");
    for(ui u = 0; u < query_graph->getVerticesCount(); u++){
        if(query_graph->getVertexDegree(u) == 1){
            printf("叶子点u%u的候选集: \n",u);
            continue;
        }
        else{
            printf("非叶点u%u的候选集: ",u);
            for(ui j = 0; j < candidates_store_[u].size(); j++)
                printf("v%u ", candidates_store_[u][j]);
            printf("\n");
        }
    }
}