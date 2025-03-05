#include <fstream>
#include <chrono>
#include "search_engine.h"
#include "computesetintersection.h"
#include "streaming_config.h"

void SearchEngine::initialize(const Graph *query_graph, const Graph *data_graph)
{
    uint32_t n = query_graph->getVerticesCount(); // 查询图点的数量
    uint32_t N = data_graph->getVerticesCount();  // 数据图点的数量
    // 初始化一下信息，分配空间
    embedding_.resize(n);
    encoded_embedding_.resize(n);
    local_idx_.resize(n);
    local_candidates_store_.resize(n);
    encoded_local_candidates_store_.resize(n);

    local_candidates_buffer1_.resize(n);
    for (auto &buffer : local_candidates_buffer1_)
    {
        buffer.resize(N);
    }

    local_candidates_buffer2_.resize(n);
    for (auto &buffer : local_candidates_buffer2_)
    {
        buffer.resize(N);
    }

    visited_ = new bool[data_graph->getVerticesCount()];
    std::fill(visited_, visited_ + data_graph->getVerticesCount(), false);

    
    level_ans.resize(n-2);
    current_solved.resize(n-2);
    current_ID_flag.resize(n-2);
    /*for(ui i = 0; i < n - 2; i++){
        current_encodev_ans[i].resize(N);
    }*/
        
    // ID_ans.resize(N);
    
    reset_performance_counters();
}

void SearchEngine::release()
{
    embedding_.clear();
    encoded_embedding_.clear();
    local_idx_.clear();
    local_candidates_store_.clear();
    encoded_local_candidates_store_.clear();
    local_candidates_buffer1_.clear();
    local_candidates_buffer2_.clear();
    delete[] visited_;
    
    level_ans.clear();
    current_solved.clear();
    // ID_ans.clear();
    current_ID_flag.clear();
}

// 优化搜索空间
uint64_t SearchEngine::optimized_search_on_reduced_query(const Graph *query_graph, OrdersPerEdge &orders, LocalViewManager &lvm,
                                                         GlobalViewManager &gvm)
{
    
    // 重置第一层的结果 
    level_ans[0] = 0;
    current_solved[0] = false;
    current_ID_flag[0].clear();
    spp::sparse_hash_map<ui,int64_t> ID_ans;
    auto &u_pruning_class = lvm.u_pruning_class;
    auto &u_v_ID = lvm.u_v_ID_map;
    auto &leaf_ID = lvm.leaf_ID;
    
    auto &order = orders.matching_order_;
    auto &bn_offset = orders.matching_order_bn_offset_;
    auto &bn = orders.matching_order_bn_;
    auto &view_mapping = orders.matching_order_view_mappings_;
    
    uint32_t target_depth = order.size(); // 目标搜索的深度，也就是节点数
    uint32_t start_vertex = order[0];     // 记录开始结点和最后一个结点
    // 这里就是取的更新边，把这两个点标记为已访问，并完成这两个点的匹配
    Edge data_edge = {lvm.get_candidate_set(orders.indexing_order_[0])[0], lvm.get_candidate_set(orders.indexing_order_[1])[0]};
    embedding_[orders.indexing_order_[0]] = data_edge.first;
    embedding_[orders.indexing_order_[1]] = data_edge.second;
    visited_[data_edge.first] = true;
    visited_[data_edge.second] = true;
    uint32_t *seeds;
    local_idx_[0].first = 0;
    if (query_graph->getVertexDegree(start_vertex) == 1)
    {
        local_idx_[0].second = compute_local_candidates_for_reduced_query(query_graph, 0, order, bn_offset, bn, view_mapping, lvm, gvm);
        seeds = local_candidates_store_[0];
    }
    else
    {
        seeds = lvm.get_candidate_set(start_vertex);                     // 这里seed存放的是 真实的候选集
        local_idx_[0].second = lvm.get_candidate_set_size(start_vertex); // 表示第0个匹配结点的 索引位置就是
    }
    // 这里很容易理解，如果要匹配的结点只有一个，那么这里只需要返回这个节点的候选集规模就可以了
    if (target_depth == 1)
    {
        visited_[data_edge.first] = false;
        visited_[data_edge.second] = false;
        return local_idx_[0].second;
    }
    // 循环遍历第0层的各个候选点
    for (local_idx_[0].first = 0; local_idx_[0].first < local_idx_[0].second; local_idx_[0].first++)
    {
        // 取出候选点
        uint32_t seed = seeds[local_idx_[0].first];
        if (visited_[seed]) // 如果已经被访问，则跳过
            continue;
        embedding_[start_vertex] = seed;
        encoded_embedding_[start_vertex] = local_idx_[0].first; // 等于当前编号      
        visited_[seed] = true;  
        ui u0 = order[0];
        if(u_pruning_class[u0] == 1){
            if(current_solved[0]){
                
                visited_[seed] = false;
                ui ID1 = leaf_ID[u0];
                level_ans[0] += ID_ans[ID1]; 
                continue;
            }
            else {
                current_solved[0] = true;
            }
        }
        else if(u_pruning_class[u0] == 2){
            ui ID1 = u_v_ID[u0][seed];
            if (current_ID_flag[0][ID1]){  // 开启剪枝操作
                level_ans[0] += ID_ans[ID1];
                visited_[seed] = false;
                continue;
            }
            else{
                current_ID_flag[0][ID1] = true;
            }
        }
        spacecunt++;
        //printf("(u%u--v%u)\n",u0,seed);
        uint32_t current_depth = 1; // 当前深度达到1,计算第一层（从0开始计数）的候选集大小，当前第一个点标记为0，候选集规模存到second
        local_idx_[current_depth].first = 0;
        local_idx_[current_depth].second = compute_local_candidates_for_reduced_query(query_graph, current_depth, order,
                                                                                      bn_offset, bn,
                                                                                      view_mapping, lvm, gvm);
        if (target_depth == 2)
        {
            level_ans[current_depth] += local_idx_[current_depth].second;
            visited_[seed] = false;
        }
        else
        { // 如果不只一层，则开始递归
            while (true)
            {
                while (local_idx_[current_depth].first < local_idx_[current_depth].second)
                {
                    uint32_t u = order[current_depth]; // 取出当前的匹配点
                    uint32_t encoded_v = encoded_local_candidates_store_[current_depth][local_idx_[current_depth].first];
                    uint32_t v = local_candidates_store_[current_depth][local_idx_[current_depth].first++];
                    encoded_embedding_[u] = encoded_v;
                    embedding_[u] = v;
                    visited_[v] = true;
                    ui idxx0 = u_pruning_class[u];
                    
                    if(idxx0 == 1){ // 为叶子结点 或者为 无后驱 结点 
                        if(current_solved[current_depth]){ // 快速剪枝，用乘法;
                            ui ID = leaf_ID[u];
                            level_ans[current_depth] += ID_ans[ID];
                            visited_[v] = false;
                            continue; 
                        }
                        else{  // 要用第一个元素 构造本层的基线
                            current_solved[current_depth] = true;
                        }
                    }
                    else if(idxx0 == 2){ // 正常的可以（部分）简直的结点层
                        ui ID = u_v_ID[u][v];
                        if(current_ID_flag[current_depth][ID]){ // 开启剪枝
                            level_ans[current_depth] += ID_ans[ID]; // 把下一层的结果往上带
                            visited_[v] = false;
                            //printf("v%u--ID %u剪枝掉\n", v, ID);
                            continue;
                        }
                        else{
                            current_ID_flag[current_depth][ID] = true;
                        }
                    }
                    //printf("(u%u--v%u)\n",u,v);
                    spacecunt++;
                    uint32_t next_depth = current_depth + 1;
                    local_idx_[next_depth].first = 0;
                    local_idx_[next_depth].second = compute_local_candidates_for_reduced_query(query_graph, next_depth,
                                                                                               order,
                                                                                               bn_offset, bn,
                                                                                               view_mapping, lvm,
                                                                                               gvm); // 计算下一层候选点数量
                    if (local_idx_[next_depth].second == 0){
                        visited_[v] = false;
                        if(idxx0 == 1){
                            ui idxx1 = leaf_ID[u];
                            ID_ans[idxx1] = 0;
                        }
                        else if(idxx0 == 2){
                            ui idxx1 = u_v_ID[u][v];
                            ID_ans[idxx1] = 0;
                        }
                    }
                    else if (next_depth == target_depth - 1)
                    {
                        ui idxx1 = u_pruning_class[u];
                        if(idxx1 == 1){
                            ui ID = leaf_ID[u];
                            ID_ans[ID] = local_idx_[next_depth].second;
                            level_ans[current_depth] += ID_ans[ID];
                        }
                        else if(idxx1 == 2){
                            ui ID = u_v_ID[u][v];
                            ID_ans[ID] = local_idx_[next_depth].second;
                            level_ans[current_depth] += ID_ans[ID];
                        }
                        else{
                            // 不需要获取ID，也不需要构造基线
                            level_ans[current_depth] += local_idx_[next_depth].second;
                        }
                        visited_[v] = false;
                    }
                    else
                    { 
                        current_depth += 1;
                    }
                }
                if (g_exit)
                {
                    return 0;
                }
                // 回溯
                current_depth -= 1;
                ui u0 = order[current_depth];
                ui v0 = embedding_[u0];
                ui idxx = u_pruning_class[u0];
                if(idxx == 1){  // 回溯到叶子/无后驱节点层
                    
                    ui ID = leaf_ID[u0];
                    // 构造基线
                    ID_ans[ID] = level_ans[current_depth + 1];
                    level_ans[current_depth] += ID_ans[ID];
                    // printf("回溯到(u%u--v%u): 根结果为%zu, 层结果%zu\n",u0, v0, ID_ans[ID],level_ans[current_depth]); // 这里出了大问题

                    // 判断下一层属于什么类型，针对性进行清空操作
                    ui u1 = order[current_depth+1];
                    ui idxxx = u_pruning_class[u1];
                    if(idxxx == 1){
                        current_solved[current_depth+1] = false;
                        ui id = leaf_ID[u1];
                        ID_ans[id] = 0;
                    }
                    else if (idxxx == 2){
                        current_ID_flag[current_depth+1].clear();
                    }
                    // 清空操作
                    level_ans[current_depth+1] = 0;
                    visited_[v0] = false;
                }
                else if(idxx == 2){
                    ui ID = u_v_ID[u0][v0];
                    // 构造基线
                    ID_ans[ID] = level_ans[current_depth+1]; //这里的level_ans[current_depth+1]数值不对，可能是连续两层的问题
                    level_ans[current_depth] += ID_ans[ID];
                    // printf("回溯到(u%u--v%u): 根结果为%zu, 层结果%zu\n",u0, v0, ID_ans[ID],level_ans[current_depth]);///*****//
                    // 清空操作
                    ui u1 = order[current_depth + 1];
                    ui idxxx = u_pruning_class[u1];
                    if(idxxx == 1){
                        current_solved[current_depth+1] = false;
                        ui id = leaf_ID[u1];
                        ID_ans[id] = 0;
                    }
                    else if(idxxx == 2){
                        current_ID_flag[current_depth+1].clear();
                    }
                    level_ans[current_depth+1] = 0;
                    visited_[v0] = false;
                }
                else {
                    // 不需要获取ID，也不需要构造基线
                    level_ans[current_depth] += level_ans[current_depth+1];
                    // printf("*回溯到(u%u--v%u): 根结果为%zu, 层结果%zu\n",u0, v0, level_ans[current_depth],level_ans[current_depth]);
                    ui u1 = order[current_depth + 1];
                    ui idxxx = u_pruning_class[u1];
                    if(idxxx == 1){
                        current_solved[current_depth+1] = false;
                        ui id = leaf_ID[u1];
                        ID_ans[id] = 0;
                    }
                    else if(idxxx == 2){
                        current_ID_flag[current_depth+1].clear();
                    }      
                    // 清空操作
                    level_ans[current_depth + 1] = 0;
                    visited_[v0] = false;
                }
                
                if (current_depth < 1)
                {
                    break;
                }
            }
        }
    }
    visited_[data_edge.first] = false;
    visited_[data_edge.second] = false;
    return level_ans[0];
}

uint32_t SearchEngine::compute_local_candidates_for_reduced_query(const Graph *query_graph, uint32_t depth,
                                                                  std::vector<uint32_t> &order,
                                                                  std::vector<uint32_t> &bn_offset,
                                                                  std::vector<uint32_t> &bn,
                                                                  std::vector<uint32_t> &view_mapping,
                                                                  LocalViewManager &lvm, GlobalViewManager &gvm)
{
    uint32_t bn_begin = bn_offset[depth];
    uint32_t bn_end = bn_offset[depth + 1];
    bool one_bn = bn_end - bn_begin == 1;
    uint32_t *lc1 = nullptr;
    uint32_t lc_count1 = 0;
    uint32_t *lc2 = nullptr;
    uint32_t lc_count2 = 0;

    if (bn_end - bn_begin == 0)
    {
        // Has no backward neighbors.
        lc1 = lvm.get_candidate_set(order[depth]);
        lc_count1 = lvm.get_candidate_set_size(order[depth]);
        lc2 = local_candidates_buffer1_[depth].data();

        lc_count2 = 0;

        for (uint32_t i = 0; i < lc_count1; ++i)
        {
            uint32_t v = lc1[i];

            if (!visited_[v])
            {
                local_candidates_buffer2_[depth][lc_count2] = i;
                lc2[lc_count2++] = v;
            }
        }
        encoded_local_candidates_store_[depth] = local_candidates_buffer2_[depth].data();
        local_candidates_store_[depth] = lc2;
        return lc_count2;
    }
    else
    {
        uint32_t u = bn[bn_begin];
        uint32_t v = embedding_[u];
        uint32_t encoded_v = encoded_embedding_[u];

        uint32_t view_id = view_mapping[bn_begin];

        // If the vertex degree is 1 (i.e., a leaf), retrieve the neighbor from global view.
        lc1 = query_graph->getVertexDegree(order[depth]) == 1 ? gvm.get_view(view_id)->get_neighbor(v, lc_count1)
                                                              : lvm.get_view(view_id)->get_neighbors(encoded_v, lc_count1);
        //        lc1 = lvm.get_view(view_id)->get_neighbors(v, lc_count1);

        lc2 = local_candidates_buffer1_[depth].data();

        if (lc_count1 == 0)
        {
            goto EXIT;
        }

        // More than one backward neighbor.
        if (bn_begin + 1 < bn_end)
        {
            bn_begin += 1;
            u = bn[bn_begin];
            encoded_v = encoded_embedding_[u];
            view_id = view_mapping[bn_begin];

            lc2 = lvm.get_view(view_id)->get_neighbors(encoded_v, lc_count2);

            uint32_t temp_count;
            uint32_t *temp_buffer = local_candidates_buffer1_[depth].data();
            comcunt++;
            ComputeSetIntersection::ComputeCandidates(lc1, lc_count1, lc2, lc_count2, temp_buffer, temp_count);

            if (temp_count == 0)
            {
                lc_count1 = 0;
                goto EXIT;
            }

            lc1 = temp_buffer;
            lc_count1 = temp_count;
            temp_buffer = local_candidates_buffer2_[depth].data();

            for (bn_begin += 1; bn_begin < bn_end; ++bn_begin)
            {
                u = bn[bn_begin];
                encoded_v = encoded_embedding_[u];
                view_id = view_mapping[bn_begin];

                lc2 = lvm.get_view(view_id)->get_neighbors(encoded_v, lc_count2);
                ComputeSetIntersection::ComputeCandidates(lc1, lc_count1, lc2, lc_count2, temp_buffer, temp_count);
                comcunt++;
                if (temp_count == 0)
                    return 0;
                std::swap(temp_buffer, lc1);
                std::swap(temp_count, lc_count1);
            }

            lc2 = temp_buffer;
        }
    }

EXIT:
    if (lc_count1 == 0)
    {
        return 0;
    }

    if (query_graph->getVertexDegree(order[depth]) == 1)
    {
        lc_count2 = 0;

        for (uint32_t i = 0; i < lc_count1; ++i)
        {
            uint32_t v = lc1[i];

            if (!visited_[v])
            {
                lc2[lc_count2++] = v;
            }
        }
        encoded_local_candidates_store_[depth] = local_candidates_buffer2_[depth].data();
        local_candidates_store_[depth] = lc2;
        return lc_count2;
    }
    else
    {
        lc_count2 = 0;
        uint32_t *candidate_set = lvm.get_candidate_set(order[depth]);
        uint32_t *encoded_buffer = one_bn ? local_candidates_buffer2_[depth].data() : lc1;
        for (uint32_t i = 0; i < lc_count1; ++i)
        {
            uint32_t encoded_v = lc1[i];
            uint32_t v = candidate_set[encoded_v];

            if (!visited_[v])
            {
                encoded_buffer[lc_count2] = encoded_v;
                lc2[lc_count2++] = v;
            }
        }

        encoded_local_candidates_store_[depth] = encoded_buffer;
        local_candidates_store_[depth] = lc2;
        return lc_count2;
    }
}