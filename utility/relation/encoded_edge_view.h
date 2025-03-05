#ifndef RAPIDMATCH_ENCODED_EDGE_VIEW_H
#define RAPIDMATCH_ENCODED_EDGE_VIEW_H

#include <vector>

// 这个类用于存储和操作 编码后的边 的信息
class EncodedEdgeView {
public:
    // Key: a vertex; Value: neighbors of the vertex in key.
    // 二维向量存储每个顶点的邻居列表
    std::vector<std::vector<uint32_t>> adj_list_; //存储邻接表
    uint32_t cardinality_;// 视图的基数，表示边的总数

public:
    EncodedEdgeView() { cardinality_ = 0; } // 构造和析构
    ~EncodedEdgeView() {}

    // 插边操作
    bool insert(uint32_t u, uint32_t v) {
        auto& neighbors = adj_list_[u];
        auto lb = std::lower_bound(neighbors.begin(), neighbors.end(), v);
        if (lb == neighbors.end() || *lb != v) {
            neighbors.insert(lb, v);
            cardinality_ += 1;
            return true;
        }
        return false;
    }

    bool remove(uint32_t u, uint32_t v) {
        auto& neighbors = adj_list_[u];
        auto lb = std::lower_bound(neighbors.begin(), neighbors.end(), v);
        if (lb != neighbors.end() && *lb == v) {
            neighbors.erase(lb);
            cardinality_ -= 1;
            return true;
        }
        return false;
    }

    uint32_t get_neighbor_num(uint32_t key) {
        return adj_list_[key].size();
    }

    uint32_t* get_neighbor(uint32_t key, uint32_t& count) {
        count = adj_list_[key].size();
        return adj_list_[key].data();
    }

    bool contains(uint32_t key) {
        return key < adj_list_.size();
    }

    uint32_t get_key_num() {
        return static_cast<uint32_t>(adj_list_.size());
    }

    uint32_t get_edge_num() {
        return cardinality_;
    }

    uint64_t memory_cost() {
        if (get_key_num() == 0)
            return 0;

        uint64_t memory_cost = 0;
        uint64_t per_element_size = sizeof(uint32_t);
        memory_cost += get_key_num() * (per_element_size + sizeof(std::vector<uint32_t>)) + get_edge_num() * per_element_size;
        return memory_cost;
    }
};

#endif //RAPIDMATCH_ENCODED_EDGE_VIEW_H
