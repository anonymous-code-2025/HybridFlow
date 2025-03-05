#ifndef RAPIDMATCH_SEARCH_ENGINE_H
#define RAPIDMATCH_SEARCH_ENGINE_H
#include "local_view_manager.h"
#include <unordered_set>
using spp::sparse_hash_set;

class SearchEngine {
public:
    uint64_t target_number = 1000;
public:
    // Performance counters
    uint64_t comcunt=0;
    uint64_t spacecunt=0;
    uint64_t invalid_partial_result_count_; // 无效的部分结果
    uint64_t partial_result_count_;         // 部分匹配结果
    uint64_t iso_conflict_count_;           // 同构冲突
    uint64_t si_empty_count_;               // 集合交集为空的计数
    uint64_t lc_empty_count_;
    
    void reset_performance_counters() {
        invalid_partial_result_count_ = 0;
        partial_result_count_ = 0;
        iso_conflict_count_ = 0;
        si_empty_count_ = 0;
        lc_empty_count_ = 0;
    }
private:
    std::vector<uint32_t*> local_candidates_store_;  // 存放局部候选节点的指针
    std::vector<uint32_t*> encoded_local_candidates_store_; // 存放编码后的候选节点

    std::vector<std::vector<uint32_t>> local_candidates_buffer1_; // 存放局部候选节点的缓冲区
    std::vector<std::vector<uint32_t>> local_candidates_buffer2_;
    std::vector<std::pair<uint32_t, uint32_t>> local_idx_; // 存放局部候选节点的开始索引和结束索引
    ui target_depth;
    
    bool* visited_; // 标记已经访问过的结点
    
    

    std::vector<uint64_t> level_ans; // 记录层结果
    // std::vector<spp::sparse_hash_map<ui,int64_t>> current_v_ans;
    // spp::sparse_hash_map<ui,spp::sparse_hash_map<ui,int64_t>> u_v_ans;
    // std::vector<std::vector<int64_t>> current_encodev_ans;
    std::vector<bool> current_solved;   // 替换u_solved
    // std::vector<int64_t> ID_ans;  // 替换leaf_ID_ans 和 ID_ans;
    std::vector<spp::sparse_hash_map<ui,bool>> current_ID_flag;
 


    // Map u to v
    std::vector<uint32_t> embedding_; // 存储当前匹配的结果

    // Map u to local encoded id
    std::vector<uint32_t> encoded_embedding_;//存放 编码后的 匹配结果

private:
    uint32_t compute_local_candidates_for_reduced_query(const Graph *query_graph, uint32_t depth,
                                                        std::vector<uint32_t> &order,
                                                        std::vector<uint32_t> &bn_offset,
                                                        std::vector<uint32_t> &bn,
                                                        std::vector<uint32_t> &view_mapping,
                                                        LocalViewManager &lvm, GlobalViewManager &gvm);
public:
    SearchEngine() {}
    ~SearchEngine() {}

    void initialize(const Graph *query_graph, const Graph *data_graph);
    void release();

    uint64_t search_on_reduced_query(const Graph *query_graph, OrdersPerEdge &orders, LocalViewManager &lvm,
                                     GlobalViewManager &gvm);
    uint64_t optimized_search_on_reduced_query(const Graph *query_graph, OrdersPerEdge &orders, LocalViewManager &lvm,
                                     GlobalViewManager &gvm);
    
};


#endif //RAPIDMATCH_SEARCH_ENGINE_H
