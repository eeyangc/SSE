#ifndef _SNP_DATA_H
#define _SNP_DATA_H

#include "common.h"

namespace SNP
{

#define DATA_LOCK std::lock_guard<std::mutex> lock(_mtx)

struct data_block_t
{
    // beta: p * K predictor
    std::vector<std::vector<real_t> > beta;

    // standard error: p * K for predictor
    std::vector<std::vector<real_t> > s;

    // correlation matrix: p * p for variable x
    std::vector<std::vector<real_t> > R;

    // sample size matrix: p * K
    std::vector<std::vector<real_t> > sample_size;

    // mask
    std::vector<bool> mask;

    // range information
    uint32_t beg;
    uint32_t end;

    data_block_t()
    {
        beta.clear();
        s.clear();
        R.clear();
        mask.clear();
        beg = 0;
        end = 0;
    }

    int init_size(uint32_t p)
    {
        if (p == 0)
        {
            return 0;
        }

        beta.clear();
        beta.resize(2);
        beta[0].resize(p, 0.0);
        beta[1].resize(p, 0.0);

        s.clear();
        s.resize(2);
        s[0].resize(p, 0.0);
        s[1].resize(p, 0.0);

        R.clear();
        R.resize(p);
        for (uint32_t i = 0; i < p; i++) {
            R[i].resize(p, 0.0);
        }

        mask.clear();
        mask.resize(p, true);

        return 0;
    }

    uint32_t size()
    {
        return beta.size();
    }

    bool operator < (const data_block_t& block) const
    {
        if (beta.size() == block.beta.size()) {
            return false;
        } else if (beta.size() < block.beta.size()) {
            return true;
        }

        return false;
    }

    bool operator > (const data_block_t& block) const
    {
        if (beta.size() == block.beta.size()) {
            return false;
        } else if (beta.size() > block.beta.size()) {
            return true;
        }

        return false;
    }
};

typedef std::vector<data_block_t>::iterator SNPDataIterator;

class SNPData
{
public:
    static SNPData* instance()
    {
        if (_snp_data == nullptr) {
            _snp_data = new (std::nothrow) SNPData();
        }

        return _snp_data;
    }

public:
    SNPData(void);
    ~SNPData(void);

    int load();

    const uint32_t dim() const
    {
        return _dim;
    }

    const uint32_t effective_dim() const
    {
        return _effective_dim;
    }

    const uint32_t col_num() const
    {
        return _col_num;
    }

    data_block_t* next()
    {
        DATA_LOCK;

        if (_current_idx >= _data_blocks.size()) {
            return nullptr;
        }

        _current_idx++;
        return &_data_blocks[_current_idx - 1];
    }

    int init_access()
    {
        DATA_LOCK;

        _current_idx = 0;

        return 0;
    }
private:
    int sort()
    {
        // TODO:确定是按照block的大小从小到大排序的
        std::sort(_data_blocks.begin(), _data_blocks.end(), std::less<data_block_t>());
        return 0;
    }

    int load_R(const SNPConfs& confs);
    int load_beta(const SNPConfs& confs);
    int load_s(const SNPConfs& confs);
    int load_sample_size(const SNPConfs& confs);

    void debug_print();

private:
    std::mutex _mtx;
    std::vector<data_block_t> _data_blocks;
    uint32_t _current_idx;

    // 记录了每个block的end index，index从0开始
    // 比如：10,21,50 (其中数据是从0,1,2,...,49，维度是50维)
    // 代表有3个range：[0, 10), [10, 21), [21, 50)
    std::vector<uint32_t> _block_range;

    // number of instance  = _data_blocks.size()
    uint32_t _dim;
    // effective size
    uint32_t _effective_dim;

    // number of cols
    uint32_t _col_num;

    // shrinkage param
    real_t _shrinkage;

    static SNPData* _snp_data;
};

}
#endif
