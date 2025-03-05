#ifndef RAPIDMATCH_LOCAL_VIEW_MANAGER_H
#define RAPIDMATCH_LOCAL_VIEW_MANAGER_H

#include "streaming_type.h"
#include "global_view_manager.h"
#include "relation/local_edge_view.h"

struct VectorHasher2 {
    size_t operator()(const std::vector<ui>& vec) const {
        size_t hash = 0;
        for (ui val : vec) {
            hash ^= std::hash<ui>()(val) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};

// 自定义相等性比较函数
struct VectorEqual2 {
    bool operator()(const std::vector<ui>& a, const std::vector<ui>& b) const {
        if (a.size() != b.size()) return false;
        return std::equal(a.begin(), a.end(), b.begin());
    }
};

// 使用自定义哈希函数和相等性比较函数的稀疏哈希表
typedef spp::sparse_hash_map<std::vector<ui>, std::vector<ui>, VectorHasher2, VectorEqual2> VectorMap;

class LocalViewManager {
public:
    // 下面三个像是指标
    uint64_t build_visited_neighbor_count_; //构建局部视图时，访问邻居的次数
    uint64_t generate_visited_neighbor_count_;//生成候选集时，访问邻居的次数
    uint64_t first_vertex_neighbor_;
     
    // std::vector<std::vector<uint32_t>> encoded_candidates_store2_;
    // std::vector<std::vector<uint32_t>> candidates_store2_; // 存放每个点的候选集
private:
    // Store the vertex id.
    std::vector<std::vector<uint32_t>> candidates_store_; // 存放每个点的候选集
    std::vector<std::vector<uint32_t>> encoded_candidates_store_; // 存放每个点的候选集
    std::vector<LocalEdgeView> views_; // 存储局部视图，核心数据结构
    std::vector<uint32_t> buffer_pool_;// 用于存储局部视图构建过程中的中间数据
    spp::sparse_hash_map<Edge, uint32_t> edge_view_mapping_; // 边到局部视图ID的映射

    std::vector<uint32_t> flag_array_; // 用于标记顶点在某个操作中的状态
    std::vector<uint32_t> updated_flag_; // 记录 在更新操作中涉及的顶点
    std::vector<uint32_t> si_buffer_; // 用于 存储集合交集操作的中间结果
    uint32_t updated_count_;        //  记录在更新操作中 处理的顶点数量

    std::vector<uint32_t*> candidate_set_pointer_; // 指向每个查询顶点候选集的指针
    std::vector<uint32_t> candidate_set_size_; //存储每个查询顶点u的候选集的大小
    
    // spp::sparse_hash_map<ui,std::vector<ui>> ID_vector; // ID 对应的向量组
    // std::vector<std::vector<std::vector<ui>>> u_v_ID_map;
    // spp::sparse_hash_map<ui,spp::sparse_hash_map<ui, ui>> u_v_ID_map; // u对应的候选点v的ID
    // spp::sparse_hash_map<ui,ui> leaf_ID;
    
    std::vector<bool> u_back;
public:
    std::vector<bool> u_pruning;
    std::vector<ui> u_pruning_class;
    std::vector<ui> leaf_ID;
    spp::sparse_hash_map<ui,spp::sparse_hash_map<ui, ui>> u_v_ID_map; // u对应的候选点v的ID
private:
    // 优化生成 局部候选集的过程
    bool optimized_generate_local_candidates(const Graph *query_graph, OrdersPerEdge &orders, 
                                             GlobalViewManager &gvm, Edge exclude_data_edge);

    bool optimized_generate_local_candidates_v2(const Graph *query_graph, OrdersPerEdge &orders,
                                                GlobalViewManager &gvm, Edge exclude_data_edge);
    
    // 优化生成局部视图的过程
    void optimized_build_local_view(const Graph *query_graph, OrdersPerEdge &orders, GlobalViewManager &gvm,Edge data_edge);
    
    void optimized_build_local_view_v2(const Graph *query_graph, OrdersPerEdge &orders, GlobalViewManager &gvm);


    // 剪枝操作，移除不符合条件的候选顶点
    bool prune_local_candidates(GlobalViewManager &gvm, uint32_t u, std::vector<uint32_t> &bn);

    // 函数： 用于更新与特定顶点相邻的候选点的状态
    void set_adjacent_update_candidates_flag(GlobalViewManager &gvm, uint32_t u, std::vector<uint32_t> &bn,
                                             uint32_t encoded_v0, uint32_t encoded_v1);

    // 函数： 设置更新候选集的标志
    void set_update_candidates_flag(GlobalViewManager &gvm, uint32_t u, std::vector<uint32_t> &bn,
                                    uint32_t encoded_v0, uint32_t encoded_v1);
    // 函数：从候选顶点列表中选择具有最小度数和的顶点
    uint32_t select_bn_with_minimum_degree_sum(GlobalViewManager &gvm, uint32_t u, std::vector<uint32_t> &bn);

    void deal_is_prunning2(const Graph *query_graph, OrdersPerEdge &orders, Edge data_edge);

    void create_equivalent_vertices_based_local_view(const Graph *query_graph, OrdersPerEdge &orders, GlobalViewManager &gvm,Edge data_edge);
    
public:
    LocalViewManager() {}
    ~LocalViewManager() { release(); }

    void initialize(const Graph *query_graph, const Graph *data_graph); // 初始化， 分配空间
    void release();

    bool create_view(const Graph *query_graph, OrdersPerEdge &orders, GlobalViewManager &gvm, Edge data_edge); // 根据Q、G 以及全局视图创建局部视图

    void destroy_view();
    void print_local_information(const Graph *query_graph);
    LocalEdgeView* get_view(uint32_t view_id); // 根据视图ID获得 局部视图
    LocalEdgeView* get_view(Edge query_edge); // 根据查询边获得 局部视图
    uint32_t get_view_id(Edge query_edge);    // 根据查询边获得 视图ID
    uint32_t* get_candidate_set(uint32_t u);  // 获取u的候选集
    uint32_t get_candidate_set_size(uint32_t u); //获取u的候选集规模


    
    uint32_t get_ID(ui u, ui v);
    std::vector<ui> get_equivalent_vertices(ui ID);
    inline u_int16_t get_leaf_ID(ui u){
        return leaf_ID[u];
    }
    void deal_is_prunning(const Graph *query_graph, OrdersPerEdge &orders);
    void local_reset(const Graph *query_graph);
};


#endif //RAPIDMATCH_LOCAL_VIEW_MANAGER_H
