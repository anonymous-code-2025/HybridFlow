#ifndef RAPIDMATCH_STREAMING_TYPE_H
#define RAPIDMATCH_STREAMING_TYPE_H

#include <cstdint>
#include <vector>

#include "sparsepp/spp.h"


/**
 * begin label, edge label, end label;
 */
// 定义一个LabelTriple的结构体，它包含三个标签
// 起始顶点标签 + 边标签 + 结束顶点标签
typedef struct LabelTriple {
    uint32_t src_label_;
    uint32_t edge_label_;
    uint32_t dst_label_;
    // 运算符重载==
    bool operator==(const LabelTriple& l) const {
        return l.src_label_ == src_label_ && l.edge_label_ == edge_label_ && l.dst_label_ == dst_label_;
    }
} LabelTriple;

/**
 * begin vertex, end vertex
 */
typedef std::pair<uint32_t, uint32_t> Edge;

// Edges and the corresponding view.
// 包含边的向量组----对应的视图索引
typedef std::pair<std::vector<Edge>, uint32_t> MappedViews;

/**
 * Orders for each edge
 */
// 枚举型，包含两种关系的边类型，普通型和全连接型
enum RelationEdgeType {
    REGULAR = 0,
    FULL_CONNECTION = 1
};

// 结构体：每个边的排序信息（存储索引顺序、匹配顺序、边类型等信息），，存储与图的边相关的多种信息
typedef struct OrdersPerEdge {
    std::vector<uint32_t> indexing_order_; // 存储与边相关的顶点索引顺序，确保顶点
    uint32_t triangle_end_index_;// 与三角形结构相关的 边的结束索引
    uint32_t adjacent_end_index_;// 与相邻边相关的 借宿索引

    std::vector<uint32_t> indexing_order_bn_offset_; // 存储与indexing_order相关的偏移量
    std::vector<uint32_t> indexing_order_bn_;// 存储索引

    std::vector<uint32_t> matching_order_; //存储匹配顺序

    std::vector<uint32_t> matching_order_bn_offset_; // 存储与matching_order相关的偏移量
    std::vector<uint32_t> matching_order_bn_; // 存储与matching_order相关的特定类型的边的匹配顺序
    
    std::vector<uint32_t> matching_order_view_mappings_; // 存储匹配顺序与视图映射之间的关系
    std::vector<RelationEdgeType> matching_order_edge_type_; // 存储与matching_order相关的边的类型信息
} OrdersPerEdge;

namespace std
{
    template<>
    struct hash<LabelTriple>{
        std::size_t operator()(LabelTriple const &l) const{
            std::size_t seed = 0;
            spp::hash_combine(seed, l.src_label_);
            spp::hash_combine(seed, l.dst_label_);
            spp::hash_combine(seed, l.edge_label_);
            return seed;
        }
    };

    template<>
    struct hash<Edge>{
        std::size_t operator()(Edge const &l) const{
            std::size_t seed = 0;
            spp::hash_combine(seed, l.first);
            spp::hash_combine(seed, l.second);
            return seed;
        }
    };
}

typedef struct Update {
    uint64_t id_;// 第几个更新
    char op_;   // 添加 or 删除
    Edge edge_; // v?-v?
    LabelTriple labels_; //标签组
} Update;


typedef struct VertexOrderingPriority {
    uint32_t bn_count_;
    uint32_t core_value_;
    uint32_t degree_;
    uint32_t vertex_id_;
    bool operator <(const VertexOrderingPriority& rhs) const {
        if (bn_count_ != rhs.bn_count_)
            return bn_count_ < rhs.bn_count_;

        if (core_value_ != rhs.core_value_)
            return core_value_ < rhs.core_value_;

        if (degree_ != rhs.degree_)
            return degree_ < rhs.degree_;

        return vertex_id_ < rhs.vertex_id_;
    }
} VertexOrderingPriority;

#endif //RAPIDMATCH_STREAMING_TYPE_H
