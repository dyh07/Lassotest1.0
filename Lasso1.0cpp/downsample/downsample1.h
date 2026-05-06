#ifndef DOWNSAMPLE_H
#define DOWNSAMPLE_H

#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>

// 模拟AnnData数据结构（不变）
struct AnnData
{
    // 1. 细胞 obs 信息
    std::vector<std::string> cell_ids;
    std::vector<int> leiden_labels;

    // 2. 基因 var 信息 
    std::vector<std::string> gene_ids;

    // 3. 表达矩阵 X [cell][gene] 
    std::vector<std::vector<float>> X;

    // 4. 降维坐标
    std::vector<std::vector<float>> umap_coords;
    std::vector<std::vector<float>> pca_coords;

    // 基础信息 
    int n_obs = 0;  // 细胞数
    int n_vars = 0; // 基因数

    // 一次性加载整套5个CSV
    bool load_full_from_csv(const std::string& folder_path);
};

class downsample
{
public:
    downsample() = default;

    // 核心下采样函数
    void runSampling(AnnData& adata,
        AnnData& adata_sampled,
        std::vector<int>& nearest_ids,
        double sample_rate = 0.01,
        double leiden_r = 1.0,
        double uniform_rate = 0.5);
    std::vector<int> get_nearest_ids(AnnData& adata, AnnData& sampled_adata);
    std::vector<int> recover_full_list(std::vector<int>& selected_list, std::vector<int>& nearest_ids);

private:
    float euclidean_distance(const std::vector<float>& a, const std::vector<float>& b);
    std::vector<int> random_choice(std::vector<int> source, int size);
};

#endif
