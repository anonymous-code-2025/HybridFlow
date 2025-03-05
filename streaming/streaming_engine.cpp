#include <future>
#include <fstream>
#include <sstream>
#include "streaming_engine.h"
#include "graphoperations.h"
#include "streaming_config.h"
#include "utility/simple_command_parser.h"
#include "utility/mem_usage.h"

#ifdef MEASURE_UPDATE_COST
    std::vector<uint64_t> g_measure_time_cost_per_update; // In nanoseconds
#endif

// #define MEASURE_INDEXING_COST

#ifdef MEASURE_INDEXING_COST
    uint64_t g_update_count;
#endif

// *********函数： 初始化
void StreamingEngine::initialize(const Graph *query_graph, const Graph *data_graph, uint64_t target_embedding_num) {
    gvm_.initialize(query_graph);//gvm是StreamingEngine类中一个成员变量，为GlobalViemManager类型的实例
    lvm_.initialize(query_graph, data_graph);
    om_.initialize(query_graph);
    sm_.initialize(query_graph, data_graph);
    sm_.target_number = target_embedding_num;
}

void StreamingEngine::release() {
    gvm_.release();
    lvm_.release();
    om_.release();
    sm_.release();
}

// *********函数： 预处理模块，计算时间内存。调用gvm.create_views(query_graph, data_graph)创建全局视图，调用om.create_orders(query_graph)生成匹配顺序
void StreamingEngine::preprocess(const Graph *query_graph, const Graph *data_graph) {
    // Create global views.
    std::cout << "Preprocess...\n";
    std::cout << "Create global views...\n";
    double initial_vm_mem_usage_KB = 0;
    double initial_rss_mem_usage_KB = 0;
    double vm_mem_usage_KB = 0;
    double rss_mem_usage_KB = 0;

    mem_usage::process_mem_usage(initial_vm_mem_usage_KB, initial_rss_mem_usage_KB);

    auto start = std::chrono::high_resolution_clock::now();

    gvm_.create_views(query_graph, data_graph);

    auto end = std::chrono::high_resolution_clock::now();
    global_view_initialize_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    mem_usage::process_mem_usage(vm_mem_usage_KB, rss_mem_usage_KB);
    //打印全局视图信息
    // gvm_.print_view_info();

    // printf("Global view create time (seconds): %.6f\n", NANOSECTOSEC(global_view_initialize_time_));

    // Generate orders for each edge.
    std::cout << "Generate processing orders...";
    //start = std::chrono::high_resolution_clock::now();

    om_.create_orders(query_graph);

    //end = std::chrono::high_resolution_clock::now();
    order_generation_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    //打印顺序信息
    om_.print_info();

    printf("Global view create time (seconds): %.6f\n", NANOSECTOSEC(global_view_initialize_time_));
    printf("Preprocessing time (seconds): %.6f\n", NANOSECTOSEC(global_view_initialize_time_ + order_generation_time_));
    printf("Global view memory cost (VM in MB): %.6f\n", (vm_mem_usage_KB - initial_vm_mem_usage_KB) / 1024.0);
    printf("Global view memory cost (RSS in MB): %.6f\n", (rss_mem_usage_KB - initial_rss_mem_usage_KB) / 1024.0);
}

// *********函数：负责执行的核心函数 针对一次更新操作  第一个标志位标志是都启用局部视图，第二个标志位标志是否执行搜索，返回匹配数量
uint64_t StreamingEngine::execute(const Graph *query_graph, const Graph *data_graph, const Update &update, bool enable_local_view = true,
                                  bool enable_search = true) {
    LabelTriple label_triple = update.labels_;//从update对象中提取三元组标签和数据边
    Edge data_edge = update.edge_;

    uint64_t result_count = 0;  // 匹配结果
    is_relevant_ = false;       // 更新边是否与查询相关，是否能在全局视图中找到相关映射
    is_searched_ = false;       // 是否执行了搜索

    if (update.op_ == '+') {
        // Update global view
        auto it = gvm_.get_mapped_views(label_triple);//获取与 给定标签三元组 相关联的全局映射，如果存在则继续
        if (it != nullptr) {
            is_relevant_ = true;
            relevant_update_count_ += 1;
            for (uint32_t i = 0; i < 2; ++i) { // 循环两次以处理可能的双向边（如果查询图和数据图是无向的）
                it = gvm_.get_mapped_views(label_triple);
                uint32_t view_id = it->second;
                gvm_.update_view(update.op_, view_id, data_edge, label_triple); // 在更新视图时，一定要更新两次
                std::swap(data_edge.first, data_edge.second);
                std::swap(label_triple.src_label_, label_triple.dst_label_);
            }
            std::swap(data_edge.first, data_edge.second);               // 交换两个点和 两个点的标签，以准备下一次循环迭代
            std::swap(label_triple.src_label_, label_triple.dst_label_);
        }
        // Update nlf view.更新NLF视图
        if (is_relevant_) {  
            gvm_.update_nlf_view(update.op_, data_edge, label_triple);
        }
        // Generate local view 
        auto mapped_automorphism = om_.get_mapped_automorphism(label_triple);  // 获取与标签三元组相关联的自同构映射，如果存在则继续
        
        // 如果找到自同构结构，标记此更新相关
        if (mapped_automorphism != nullptr) {
            is_relevant_ = true;
            if (enable_local_view) { // 如果开启了局部视图，则继续

                for (auto automorphism_id: *mapped_automorphism) { // 循环遍历每个自同构映射ID 
                    // 获取自同构的元数据，包括匹配顺序和自同构大小
                    auto meta = om_.get_automorphism_meta(automorphism_id); // 获取对应的自同构边集的第一条边，对应的order[id]，和对应边集的大小
                    OrdersPerEdge *orders = std::get<1>(meta);
                    uint32_t automorphism_size = std::get<2>(meta);
                    // ****创建局部视图，根据查询图 + 对应的自同构边集的order + 全局视图 + 更新边
                    //auto start2 = std::chrono::high_resolution_clock::now();
                    bool is_valid = lvm_.create_view(query_graph, *orders, gvm_, data_edge); 
                    //auto end2 = std::chrono::high_resolution_clock::now();
                    //equivelnt_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end2 - start2).count();
                    // printf("Local view create time (seconds): %.6f\n",NANOSECTOSEC(time2));
                    if (is_valid) { // 如果局部视图创建有效，更新与搜索相关的计数器
                        search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_;
                        search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_;
                    } else { // 如果局部视图创建失败，更新与 非搜索相关的计数器
                        non_search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_;
                        non_search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_;

                        if (lvm_.generate_visited_neighbor_count_ == 0) // 如果未生成任何邻居，增加直接拒绝计数
                            direct_rejection_count_ += 1;
                    }
                    first_indexing_vertex_ += lvm_.first_vertex_neighbor_;//？？//
                    if (enable_search && is_valid) { // 如果启用本地搜索并且创建视图有效，则开启搜索
                        //auto start3 = std::chrono::high_resolution_clock::now();
                        //sm_.comcunt = 0;
                        //sm_.spacecunt = 0;
                        uint64_t local_result_count = sm_.optimized_search_on_reduced_query(query_graph, *orders, lvm_, gvm_);
                        result_count += local_result_count * automorphism_size;
                        is_searched_ = true;
                        //auto end3 = std::chrono::high_resolution_clock::now();
                        //uint64_t recu_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(end3 - start3).count();
                        //recu_time2 += recu_time_;
                        //comcount += sm_.comcunt;
                        //spacecount+= sm_.spacecunt;
                        //printf("(v%u-v%u)    result--%zu  ",data_edge.second,data_edge.first,local_result_count*automorphism_size);
                        //printf("recu time (seconds): %.6f\n",NANOSECTOSEC(recu_time_));
                    }
                    //lvm_

                    lvm_.destroy_view();

                    if (result_count >= sm_.target_number) break;
                    if (g_exit) break;
                }
                
            }
        }
    }
    else {
        // Generate local view
        auto mapped_automorphism = om_.get_mapped_automorphism(label_triple);
        if (mapped_automorphism != nullptr) {
            is_relevant_ = true;
            if (enable_local_view) {

                for (auto automorphism_id: *mapped_automorphism) {
                    auto meta = om_.get_automorphism_meta(automorphism_id);
                    OrdersPerEdge *orders = std::get<1>(meta);
                    uint32_t automorphism_size = std::get<2>(meta);

                    bool is_valid = lvm_.create_view(query_graph, *orders, gvm_, data_edge);

                    if (is_valid) {
                        search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_;
                        search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_;
                    } else {
                        non_search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_;
                        non_search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_;

                        if (lvm_.generate_visited_neighbor_count_ == 0)
                            direct_rejection_count_ += 1;
                    }
                    first_indexing_vertex_ += lvm_.first_vertex_neighbor_;
                    if (enable_search && is_valid) {
                        uint64_t local_result_count = sm_.optimized_search_on_reduced_query(query_graph, *orders, lvm_, gvm_);
                        result_count += local_result_count * automorphism_size;
                        is_searched_ = true;
                    }

                    lvm_.destroy_view();

                    if (result_count >= sm_.target_number) break;
                    if (g_exit) break;
                }
            }
        }
        // Update global view
        auto it = gvm_.get_mapped_views(label_triple);
        if (it != nullptr) {
            is_relevant_ = true;
            relevant_update_count_ += 1;
            for (uint32_t i = 0; i < 2; ++i) {
                it = gvm_.get_mapped_views(label_triple);

                uint32_t view_id = it->second;
                gvm_.update_view(update.op_, view_id, data_edge, label_triple);

                std::swap(data_edge.first, data_edge.second);
                std::swap(label_triple.src_label_, label_triple.dst_label_);
            }

            std::swap(data_edge.first, data_edge.second);
            std::swap(label_triple.src_label_, label_triple.dst_label_);
        }

        // Update nlf view.
        if (is_relevant_) {
            gvm_.update_nlf_view(update.op_, data_edge, label_triple);
        }
    }
    // Update performance counters
    edge_process_count_ += 1;
    if (is_relevant_) {

        if (result_count > 0) {
            positive_count_ += 1;
        }

        if (is_searched_) {
            search_count_ += 1;
        }

        result_count_ += result_count;
        invalid_partial_result_count_ += sm_.invalid_partial_result_count_;
        partial_result_count_ += sm_.partial_result_count_;
        iso_conflict_count_ += sm_.iso_conflict_count_;
        si_empty_count_ += sm_.si_empty_count_;
        lc_empty_count_ += sm_.lc_empty_count_;

        sm_.reset_performance_counters();
    }
    return result_count;
}

// *********函数： 评估 视图更新的时间，开启预编译时，才会调用
void StreamingEngine::evaluate_view_update(const Graph *query_graph, const Graph *data_graph,
                                           const std::vector<Update> &stream) {
    const uint32_t repetition_time = 1;
    {
        global_view_update_time_ = 0;   //初始化和重置性能计数器
        gvm_.release();
        // Get global view update time.
        for (int i = 0; i < repetition_time; ++i) {  //重复执行更新操作以测量平均时间
            gvm_.initialize(query_graph); //初始化全局视图管理器
            gvm_.create_views(query_graph, data_graph);

            auto start = std::chrono::high_resolution_clock::now();

#ifdef MEASURE_INDEXING_COST
            uint32_t local_update_count = 0;
#endif

            for (auto& update : stream) {
                execute(query_graph, data_graph, update, false, false);

#ifdef MEASURE_INDEXING_COST
                local_update_count += 1;
                if (local_update_count >= g_update_count) break;
#endif
            }

            auto end = std::chrono::high_resolution_clock::now();
            global_view_update_time_ += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

            gvm_.release();
        }

        global_view_update_time_ /= repetition_time;
        printf("Global view update time (seconds): %.6f\n", NANOSECTOSEC( global_view_update_time_));
    }
    {
        local_view_update_time_ = 0;

        search_build_neighbor_count_ = 0;
        non_search_build_neighbor_count_ = 0;
        search_generate_neighbor_count_ = 0;
        non_search_generate_neighbor_count_ = 0;
        direct_rejection_count_ = 0;
        first_indexing_vertex_ = 0;
        for (int i = 0; i < repetition_time; ++i) {
            gvm_.initialize(query_graph);
            gvm_.create_views(query_graph, data_graph);

#ifdef MEASURE_INDEXING_COST
            uint32_t local_update_count = 0;
#endif

            auto start = std::chrono::high_resolution_clock::now();
            for (auto& update : stream) {
                execute(query_graph, data_graph, update, true, false);

#ifdef MEASURE_INDEXING_COST
                local_update_count += 1;
                if (local_update_count >= g_update_count) break;
#endif
            }

            auto end = std::chrono::high_resolution_clock::now();
            local_view_update_time_ += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

            gvm_.release();
        }

        local_view_update_time_ /= repetition_time;
        local_view_update_time_ = local_view_update_time_ > global_view_update_time_ ?
                local_view_update_time_ - global_view_update_time_ : 0;
        printf("Local view update time (seconds): %.6f\n", NANOSECTOSEC(local_view_update_time_));
    }

    printf("View update time (seconds): %.6f\n", NANOSECTOSEC(global_view_update_time_ + local_view_update_time_));
}

// *********函数：评估搜索的性能，对更新流中的每个更新 执行搜索操作； 调用execute(query_graph, update, true,true)
void StreamingEngine::evaluate_search(const Graph *query_graph, const Graph *data_graph, const std::vector<Update> &stream) {
    gvm_.release();
    gvm_.initialize(query_graph);               // 全局视图初始化
    gvm_.create_views(query_graph, data_graph); // 创建全局视图

    // 性能计数器重置
    reset_performance_counters();
    auto start = std::chrono::high_resolution_clock::now(); //开始搜索时间
    // 执行更新并启用本地视图和搜索
    for (auto& update : stream) {
        execute(query_graph, data_graph, update, true, true);
        if (g_exit) { break; }
    }
    auto end = std::chrono::high_resolution_clock::now();
    query_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
}

// **********函数： 打印结果信息，输出各种时间
void StreamingEngine::print_metrics() {
    uint64_t view_update_time = global_view_update_time_ + local_view_update_time_;
    uint64_t search_time = query_time_ > view_update_time ? query_time_ - view_update_time : 0;
    printf("Execution Summary:\n");
    printf("Preprocessing time (seconds): %.6f\n", NANOSECTOSEC(global_view_initialize_time_ + order_generation_time_));
    printf("Equivelnt time (seconds): %.6f\n", NANOSECTOSEC(equivelnt_time));
    printf("Recu time (seconds): %.6f\n", NANOSECTOSEC(recu_time2));
    printf("Query time (seconds): %.6f\n", NANOSECTOSEC(query_time_));
    printf("Global view update time (seconds): %.6f\n", NANOSECTOSEC(global_view_update_time_));
    printf("Local view update time (seconds): %.6f\n", NANOSECTOSEC(local_view_update_time_));
    printf("Search time (seconds): %.6f\n", NANOSECTOSEC(search_time));
    printf("compute setintersection count: %zu\n",comcount);
    printf("space count: %zu\n",spacecount);
    printf("Edge process count: %zu\n", edge_process_count_);
    printf("Relevant update count: %zu\n", relevant_update_count_);
    printf("Search count: %zu\n", search_count_);
    printf("Positive update count: %zu\n", positive_count_);
    printf("Result count: %zu\n", result_count_);
    printf("Partial result count (exclude results): %zu\n", partial_result_count_);
    printf("Invalid partial count: %zu\n", invalid_partial_result_count_);
    printf("ISO conflict count: %zu\n", iso_conflict_count_);
    printf("SI empty count: %zu\n", si_empty_count_);
    printf("LC empty count: %zu\n", lc_empty_count_);

    
}

// **********函数：从文件中加载更新流
void load_stream(const std::string& file_path, std::vector<Update>& stream, const Graph* data_graph) {
    auto start = std::chrono::high_resolution_clock::now();

    uint32_t vertex_num = data_graph->getVerticesCount();
    spp::sparse_hash_map<uint32_t, uint32_t> new_vertex_label;
    std::vector<uint32_t> vertex_label;
    vertex_label.reserve(vertex_num);
    for (uint32_t u = 0; u < vertex_num; ++u) {
        vertex_label.push_back(data_graph->getVertexLabel(u));
    }

    std::ifstream ifs(file_path);

    if (!ifs.is_open()) {
        std::cout << "Can not open the stream file " << file_path << " ." << std::endl;
        exit(-1);
    }

    Update update;
    while (ifs.good()) {
        std::string tmp_str;
        std::stringstream ss;
        std::string op_str;
        std::getline(ifs, tmp_str);

        if (tmp_str[0] != '#') {
            ss.clear();
            ss << tmp_str;
            ss >> op_str;

            if (op_str == "v") {
                uint32_t id;
                uint32_t label;
                ss >> id >> label;
                if (id < vertex_num) {
                    vertex_label[id] = label;
                }
                else {
                    new_vertex_label[id] = label;
                }
            }
            else if (op_str == "e" || op_str == "-e") {
                update.op_ = op_str == "e" ? '+' : '-';
                ss >> update.edge_.first >> update.edge_.second >> update.labels_.edge_label_;

                if (update.edge_.first < vertex_num) {
                    update.labels_.src_label_ = vertex_label[update.edge_.first];
                }
                else {
                    update.labels_.src_label_ = new_vertex_label[update.edge_.first];
                }

                if (update.edge_.second < vertex_num) {
                    update.labels_.dst_label_ = vertex_label[update.edge_.second];
                }
                else {
                    update.labels_.dst_label_ = new_vertex_label[update.edge_.second];
                }

                stream.emplace_back(update);
            }
        }
    }

    ifs.close();

    auto end = std::chrono::high_resolution_clock::now();
    auto load_stream_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    printf("Load stream file time (seconds): %.6f\n", NANOSECTOSEC(load_stream_time));
}

bool g_exit = false; // 全局退出标志，标志搜索的过程是否要提前退出，标志完成还是超时

// **********函数：sm是执行StreamingEngine对象的指针，该对象负责协调整个搜索过程，调用evaluate_search函数
// 该函数通过使用 C++11的异步特性 来实现对搜索过程的时间限制，用g_exit来标记
void execute_within_time_limit(const Graph *query_graph, const Graph *data_graph, std::vector<Update> &stream,
                               StreamingEngine *sm, uint64_t time_limit) {
    g_exit = false;
    auto future = std::async(std::launch::async, [query_graph, data_graph, &stream, sm](){//异步执行 搜索评估
        sm->evaluate_search(query_graph, data_graph, stream);
    });

    std::future_status status;
    do {
        status = future.wait_for(std::chrono::seconds(time_limit));
        if (status == std::future_status::deferred) {
            std::cout << "Deferred\n";
            exit(-1);
        } else if (status == std::future_status::timeout) {
            g_exit = true;
        }
    } while (status != std::future_status::ready);
}

std::string getFileNameWithoutExtension(const std::string& s) {
    char sep = '/';
#ifdef _WIN32
    sep = '\\';
#endif
    size_t i = s.rfind(sep, s.length());
    if (i != std::string::npos)
    {
        std::string filename = s.substr(i+1, s.length() - i);
        size_t lastindex = filename.find_last_of(".");
        std::string rawname = filename.substr(0, lastindex);
        return rawname;
    }
    return "";
}

int main(int argc, char** argv) {
    InputParser cmd_parser(argc, argv); // 创建命令行参数解析器
    // 从命令行参数中获取数据图文件路径、更新流文件路径、查询图文件路径、目标数量、步长、时间限制
    std::string input_data_graph_file = cmd_parser.get_cmd_option("-d");
    std::string input_data_graph_update_file = cmd_parser.get_cmd_option("-u");
    std::string input_query_graph_file = cmd_parser.get_cmd_option("-q");
    std::string input_target_embedding_number = cmd_parser.get_cmd_option("-num");
    std::string input_step_length = cmd_parser.get_cmd_option("-s");
    std::string input_time_limit = cmd_parser.get_cmd_option("-time_limit");

    // 将字符串参数转换为数值，若转换失败则使用默认值
    uint64_t target_embedding_number = std::numeric_limits<uint64_t>::max();
    if (!input_target_embedding_number.empty()) {
        if (input_target_embedding_number != "MAX") {
            target_embedding_number = std::stoll(input_target_embedding_number);
        }
    }

    uint32_t step_length = 10000;
    if (!input_step_length.empty()) {
        step_length = std::stoul(input_step_length);
    }

    uint64_t time_limit = 3600;
    if (!input_time_limit.empty()) {
        time_limit = std::stoll(input_time_limit);
    }
    /**
    * Output the command line information.
    */
    //打印命令行信息
    std::cout << "Command Line:" << '\n';
    std::cout << "\tData Graph: " << input_data_graph_file << '\n';
    std::cout << "\tData Graph Update: " << input_data_graph_update_file << '\n';
    std::cout << "\tQuery Graph: " << input_query_graph_file << '\n';
    std::cout << "\tTarget Number of Embeddings: " << target_embedding_number << '\n';
    std::cout << "\tTime limit (seconds): " << time_limit << '\n';
    std::cout << "--------------------------------------------------------------------" << std::endl;

#ifdef MEASURE_UPDATE_COST
    g_measure_time_cost_per_update.reserve(1000000);
#endif
    // 创建图对象并加载数据
    Graph* data_graph = new Graph(false);
    data_graph->is_edge_labeled = true;
    data_graph->loadGraphFromFileWithoutMeta(input_data_graph_file);
    Graph* query_graph = new Graph(false);
    query_graph->is_edge_labeled = true;
    query_graph->loadGraphFromFileWithoutMeta(input_query_graph_file);

    // 从文件中加载更新流
    std::vector<Update> stream;
    load_stream(input_data_graph_update_file, stream, data_graph);
    //打印图的元数据信息
    std::cout << "Data Graph:\n";
    data_graph->printGraphMetaData();
    std::cout << "Query Graph:\n";
    query_graph->printGraphMetaData();
    std::cout << "Update:\nUpdates num: " << stream.size() << '\n';
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "Preprocess....\n";

    //初始化StreamingEngine类对象、并进行初始化和预处理，负责协调整个搜索过程
    StreamingEngine streaming_engine;
    streaming_engine.initialize(query_graph, data_graph, target_embedding_number);
    streaming_engine.preprocess(query_graph, data_graph);

    //执行搜索性能评估、打印评估结果
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "Evaluate the performance of search...\n";
    execute_within_time_limit(query_graph, data_graph, stream, &streaming_engine, time_limit);

#ifdef MEASURE_INDEXING_COST
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "Evaluate the performance of view update...\n";
    streaming_engine.evaluate_view_update(query_graph, data_graph, stream);
#endif

    std::cout << "--------------------------------------------------------------------" << std::endl;
    streaming_engine.print_metrics();
    std::cout << "--------------------------------------------------------------------" << std::endl;

#ifdef MEASURE_UPDATE_COST
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "Dump measure update results..." << std::endl;

    {
        auto file_name = getFileNameWithoutExtension(input_query_graph_file);
        std::string file_path = file_name + "_time_cost_update.bin";
        std::ofstream ofs(file_path, std::ios::binary);

        // Format: The first 8 byte records the number of elements, the second 4 byte records the batch size, then append the vector content.
        uint64_t num_element = g_measure_time_cost_per_update.size();
        uint32_t batch_size = MEASURE_BATCH_SIZE;
        int64_t vec_size = sizeof(uint64_t) * num_element;

        ofs.write(reinterpret_cast<const char *>(&num_element), 8);
        ofs.write(reinterpret_cast<const char *>(&batch_size), 4);
        ofs.write(reinterpret_cast<const char *>(&g_measure_time_cost_per_update.front()), vec_size);
    }

    std::cout << "--------------------------------------------------------------------" << std::endl;
#endif

    if (g_exit)
        std::cout << "Time out..." << std::endl;
    std::cout << "Done." << std::endl;

    delete query_graph;
    delete data_graph;

    return 0;
}
