#ifndef INTERNALRAPIDMATCH_VIEW_H
#define INTERNALRAPIDMATCH_VIEW_H

#include <vector>
#include <algorithm>
#include "../sparsepp/spp.h"
using spp::sparse_hash_map;

class GlobalEdgeView {
public:
    // Key: 顶点ID; Value: 该顶点的所有邻居
    // 存放一个点----对应的所有邻居
    sparse_hash_map<uint32_t, std::vector<uint32_t>> trie_;
    uint32_t cardinality_; // 表示全局视图中边的数量

public:
    GlobalEdgeView() { cardinality_ = 0; }

    ~GlobalEdgeView() {}

    // *函数：用于在全局视图中插入一条从顶点u-v的边
    bool insert(uint32_t u, uint32_t v) {
        auto it = trie_.find(u);  // 先找到 u 这个点
        if (it != trie_.end()) {  // 如果u已经存在
            // 找到合适的位置
            auto lb = std::lower_bound(it->second.begin(), it->second.end(), v);
            if (lb == it->second.end() || *lb != v) {
                it->second.insert(lb, v);
                cardinality_ += 1;
                return true;
            }
            return false;
        }
        else { // 如果顶点u不存在，使用emplace在哈希表中创建新的条目，添加v
            auto temp_it = trie_.emplace(u, std::vector<uint32_t>());
            temp_it.first->second.push_back(v);
            cardinality_ += 1;
            return true;
        }
    }
    //*函数，从全局视图中移除一条从u->v的边
    bool remove(uint32_t u, uint32_t v) {
        auto iter = trie_.find(u);
        if (iter != trie_.end()) {
            auto lb = std::lower_bound(iter->second.begin(), iter->second.end(), v);
            if (lb != iter->second.end() && *lb == v) {
                iter->second.erase(lb);
                cardinality_ -= 1;
                return true;
            }
        }

        return false; //如果顶点 u 存在，尝试找到邻居 v 并移除。
    }
    //*函数：获取顶点k的邻居数量
    uint32_t get_neighbor_num(uint32_t key) {
        auto iter = trie_.find(key);
        if (iter != trie_.end())
            return iter->second.size();
        return 0;
    }
    //*函数：获取顶点k的邻居列表
    uint32_t* get_neighbor(uint32_t key, uint32_t& count) {
        count = 0;
        auto iter = trie_.find(key);
        if (iter != trie_.end()) {
            count = iter->second.size();
            return iter->second.data();
        }
        return nullptr;
    }
    //*函数：是否包含顶点key
    bool contains(uint32_t key) {
        return trie_.contains(key);
    }
    //*函数：获取顶点key的数量
    uint32_t get_key_num() {
        return static_cast<uint32_t>(trie_.size());
    }
    //*函数：获取边的数量
    uint32_t get_edge_num() {
        return cardinality_;
    }
    //*函数：计算全局视图的内存使用量
    uint64_t memory_cost() {
        if (get_key_num() == 0)
            return 0;

        uint64_t memory_cost = 0;
        uint64_t per_element_size = sizeof(uint32_t);
        memory_cost += get_key_num() * (per_element_size + sizeof(std::vector<uint32_t>)) * trie_.load_factor() + get_edge_num() * per_element_size;
        return memory_cost;
    }
};


#endif //INTERNALRAPIDMATCH_VIEW_H
