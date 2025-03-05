#ifndef RAPIDMATCH_ORDER_MANAGER_H
#define RAPIDMATCH_ORDER_MANAGER_H

#include "graph/graph.h"
#include "streaming_type.h"

class OrderManager {
private:
    // 存储查询图中每个顶点的所有自同构映射。自同构是指图到其自身的同构，即图中顶点的排列，使得边和标签保持一致。
    std::vector<std::vector<uint32_t>> automorphisms_;
    std::vector<OrdersPerEdge> orders_; // 存储与每条查询边相关的顺序信息
    spp::sparse_hash_map<Edge, uint32_t> edge_orders_mapping_; // 存储 每条查询边 映射到 orders中的索引，便于查找特定边的顺序信息
    spp::sparse_hash_map<LabelTriple, std::vector<Edge>> label_edge_mapping_; // 标签组 到 边集合 的映射
    spp::sparse_hash_map<LabelTriple, std::vector<uint32_t>> label_automorphism_mapping_; // 标签组 到 自同构集合的索引的映射，便于查找标签快速检索相关的自同构信息
    std::vector<std::vector<Edge>> automorphism_edges_; // 存储与查询图的自同构相关的边集合。每个内部向量可能代表一个自同构类中的边集合。

private:
    // New framework.
    // 检测查询图中的自同构边
    void detect_automorphism_edges(const Graph *query_graph);
    // 为简化后的查询图创建索引顺序，通常是基于某条边生成的
    void create_indexing_order_for_reduced_graph(const Graph *query_graph, Edge edge, OrdersPerEdge &order);
    // 为简化后的查询图创建匹配顺序
    void create_matching_order_for_reduced_graph(const Graph *query_graph, Edge edge, OrdersPerEdge &order);
    // 使用回溯法生成匹配顺序
    void generate_matching_order_with_RI(const Graph* graph, std::vector<uint32_t>& matching_order);

public:
    OrderManager() {}
    ~OrderManager() { release(); }

    void initialize(const Graph *query_graph);
    void release();
    // 获取与特定边相关的顺序信息
    OrdersPerEdge* get_orders(Edge edge);
    // 根据标签组获得对应的边的集合
    std::vector<Edge>* get_mapped_edges(LabelTriple label_triple);

    // 获取与特定标签组相关联的自同构索引集合
    std::vector<uint32_t>* get_mapped_automorphism(LabelTriple labelTriple);
    // 获取与特定自同构ID相关的元数据， 包含 对应的边，以及边的顺序信息指针， 和自同构的数量
    std::tuple<Edge, OrdersPerEdge*, uint32_t> get_automorphism_meta(uint32_t id);

    // 根据查询图 创建匹配顺序
    void create_orders(const Graph *query_graph);
    void print_info();
};


#endif //RAPIDMATCH_ORDER_MANAGER_H
