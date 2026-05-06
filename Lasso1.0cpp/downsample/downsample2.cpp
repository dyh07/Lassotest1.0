#include "downsample1.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <random>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace std;

// 工具：读取单行CSV
static std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> res;
    std::stringstream ss(line);
    std::string val;
    while (std::getline(ss, val, ',')) res.push_back(val);
    return res;
}

// 核心：加载全套 5 个 CSV
bool AnnData::load_full_from_csv(const std::string& folder_path)
{
    try {
        //1. 读取 obs.csv (细胞信息)
        std::ifstream obs(folder_path + "/obs.csv");
        if (!obs.is_open()) return false;
        std::string line;
        std::getline(obs, line);
        while (std::getline(obs, line)) {
            auto parts = split_csv(line);
            cell_ids.push_back(parts[0]);
            leiden_labels.push_back(std::stoi(parts.back()));
        }
        n_obs = cell_ids.size();

        //2. 读取 var.csv (基因信息)
        std::ifstream var(folder_path + "/var.csv");
        if (var.is_open()) {
            std::getline(var, line);
            while (std::getline(var, line)) {
                auto parts = split_csv(line);
                gene_ids.push_back(parts[0]);
            }
        }
        n_vars = gene_ids.size();

        //3. 读取 X.csv (表达矩阵)
        std::ifstream X_mat(folder_path + "/X.csv");
        if (X_mat.is_open()) {
            std::getline(X_mat, line);
            while (std::getline(X_mat, line)) {
                auto parts = split_csv(line);
                std::vector<float> row;
                for (auto& p : parts) row.push_back(std::stof(p));
                X.push_back(row);
            }
        }

        //4. 读取 X_umap.csv
        std::ifstream umap(folder_path + "/X_umap.csv");
        if (umap.is_open()) {
            std::getline(umap, line);
            while (std::getline(umap, line)) {
                auto parts = split_csv(line);
                umap_coords.push_back({ std::stof(parts[0]), std::stof(parts[1]) });
            }
        }

        //5. 读取 X_pca.csv
        std::ifstream pca(folder_path + "/X_pca.csv");
        if (pca.is_open()) {
            std::getline(pca, line);
            while (std::getline(pca, line)) {
                auto parts = split_csv(line);
                std::vector<float> vec;
                for (auto& p : parts) vec.push_back(std::stof(p));
                pca_coords.push_back(vec);
            }
        }

        return n_obs > 0;
    }
    catch (...) {
        return false;
    }
}

// 欧氏距离
float downsample::euclidean_distance(const std::vector<float>& a, const std::vector<float>& b)
{
    float sum = 0;
    for (int i = 0; i < a.size(); i++) sum += pow(a[i] - b[i], 2);
    return sqrt(sum);
}

// 随机不重复采样
std::vector<int> downsample::random_choice(std::vector<int> source, int size)
{
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(source.begin(), source.end(), g);
    std::vector<int> res;
    for (int i = 0; i < size && i < source.size(); i++) res.push_back(source[i]);
    return res;
}

// 完整下采样
void downsample::runSampling(
    AnnData& adata,
    AnnData& adata_sampled,
    std::vector<int>& nearest_ids,
    double sample_rate,
    double leiden_r,
    double uniform_rate
)
{
    std::unordered_set<int> clusters(adata.leiden_labels.begin(), adata.leiden_labels.end());
    int num_clusters = clusters.size();

    int sample_num = std::min((int)(adata.n_obs * sample_rate) + num_clusters * 3, adata.n_obs);
    uniform_rate = std::clamp(uniform_rate, 0.0, 1.0);

    int n_uniform = sample_num * uniform_rate;
    int n_balanced = sample_num - n_uniform;

    std::vector<int> all_idx(adata.n_obs);
    for (int i = 0; i < adata.n_obs; i++) all_idx[i] = i;

    // 均匀采样
    std::vector<int> uniform_idx = random_choice(all_idx, n_uniform);

    // 均衡采样
    std::vector<int> balanced_idx;
    if (n_balanced > 0) {
        int n_per = std::max(1, n_balanced / num_clusters);
        for (int c : clusters) {
            std::vector<int> idx;
            for (int i = 0; i < adata.n_obs; i++)
                if (adata.leiden_labels[i] == c) idx.push_back(i);
            auto s = random_choice(idx, n_per);
            balanced_idx.insert(balanced_idx.end(), s.begin(), s.end());
        }
    }

    // 合并去重
    std::unordered_set<int> final_set(uniform_idx.begin(), uniform_idx.end());
    final_set.insert(balanced_idx.begin(), balanced_idx.end());
    std::vector<int> final_idx(final_set.begin(), final_set.end());

    // 构建采样后 AnnData
    adata_sampled.n_obs = final_idx.size();
    adata_sampled.n_vars = adata.n_vars;
    adata_sampled.gene_ids = adata.gene_ids;

    for (int i : final_idx) {
        adata_sampled.cell_ids.push_back(adata.cell_ids[i]);
        adata_sampled.leiden_labels.push_back(adata.leiden_labels[i]);
        adata_sampled.umap_coords.push_back(adata.umap_coords[i]);
        adata_sampled.pca_coords.push_back(adata.pca_coords[i]);
        adata_sampled.X.push_back(adata.X[i]);
    }

    nearest_ids = get_nearest_ids(adata, adata_sampled);
}

// 近邻查找
std::vector<int> downsample::get_nearest_ids(AnnData& adata, AnnData& sampled_adata)
{
    std::vector<int> res;
    for (auto& p : adata.umap_coords) {
        float min_d = 1e9;
        int best = 0;
        for (int i = 0; i < sampled_adata.n_obs; i++) {
            float d = euclidean_distance(p, sampled_adata.umap_coords[i]);
            if (d < min_d) { min_d = d; best = i; }
        }
        res.push_back(best);
    }
    return res;
}

// 恢复全列表
std::vector<int> downsample::recover_full_list(std::vector<int>& selected_list, std::vector<int>& nearest_ids)
{
    std::unordered_set<int> sel(selected_list.begin(), selected_list.end());
    std::vector<int> res;
    for (int i = 0; i < nearest_ids.size(); i++)
        if (sel.count(nearest_ids[i])) res.push_back(i);
    return res;
}