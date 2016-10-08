#ifndef _SNP_WORKER_H
#define _SNP_WORKER_H

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "common.h"
#include "SNPData.h"

namespace SNP
{

#define WORKER_LOCK std::lock_guard<std::mutex> lock(_worker_mtx)

// register parameters of worker for different tasks here
struct worker_param_t
{
    real_t heritability_beta;
    std::vector<std::vector<real_t> > co_heritability_beta;
    gsl_rng* random_gen;

    worker_param_t() {
        random_gen = nullptr;
        heritability_beta = 0.0;
        co_heritability_beta.resize(2);
        co_heritability_beta[0].resize(2, 0.0);
        co_heritability_beta[1].resize(2, 0.0);
    }
};

// register parameters of merger for different tasks here
struct merger_param_t
{
    std::vector<real_t> heritability_heritability;
    std::vector<real_t> heritability_heritability1;
    std::vector<real_t> heritability_heritability2;
};

class SNPWorker
{
public:
    static SNPWorker* intance()
    {
        if (_snp_worker == nullptr) {
            _snp_worker = new (std::nothrow) SNPWorker();
        }
        
        return _snp_worker;
    }

    int init();

public:
    SNPWorker(void);
    ~SNPWorker(void);

    int run();
private:
    void set_error_no(uint32_t error_no)
    {
        WORKER_LOCK;

        _error_no = error_no;
    }

    const uint32_t get_error_no() const
    {
        return _error_no;
    }

    /*
    ** worker routines
    */
    void worker(uint32_t iter, worker_param_t& param);
    // routine for gibbs sampling
    void gibbs_summary_stat(data_block_t* block, uint32_t iter, worker_param_t& param);
    void gibbs_coheritability(data_block_t* block, uint32_t iter, worker_param_t& param);
    void update_prior(uint32_t iter, std::vector<worker_param_t>& params);

    /*
    ** merger routines
    */
    void merger(merger_param_t& param);
    int dump_result(std::vector<merger_param_t>& params);

    int get_param_iter_index(uint32_t iter)
    {
        if (iter < _num_burnin) {
            return -1;
        }

        int ret_res = 0;
        uint32_t index = iter - _num_burnin;
        if (index % NUM_ITER_PER_GROUP == 0) {
            ret_res = static_cast<int>(index / NUM_ITER_PER_GROUP);
        } else {
            ret_res = -1;
        }

#ifdef _DEBUG
        fprintf(stdout, "iter-[%u] to param index-[%d]\n", iter, ret_res);
#endif

        return ret_res;
    }

private:
    uint32_t _num_threads;
    uint32_t _num_iter;
    uint32_t _num_burnin;

    // result for N iterations: heritability
    std::vector<std::vector<real_t> > _beta;

    // result for co_heritability
    std::vector<std::vector<real_t> > _beta1;
    std::vector<std::vector<real_t> > _beta2;

    static SNPWorker* _snp_worker;

    gsl_rng* _r_global;
    std::mutex _worker_mtx;
    uint32_t _error_no;

private:
    // for heritability
    bool _flag_single_trait;
    std::vector<real_t> _compute_buffer;
    real_t _sigma_2_beta;

    // for heritability
    bool _flag_heritability;
    // for co_heritability:
    bool _flag_two_traits; 

    // for cache
    bool _flag_cache;
    //2 by dim
    std::vector<std::vector<real_t> > _co_compute_buffer;
    // for co_heritability: 2 by 2
    std::vector<std::vector<real_t> > _co_lambda;
    // co-heritability: _num_iter * 1
    std::vector<real_t> _rho;
    std::vector<real_t> _sigma1;
    std::vector<real_t> _sigma2;

    // for sampling
    real_t _nu;

    std::vector<worker_param_t> _worker_params;
    std::vector<merger_param_t> _merger_params;
};

}

#endif
