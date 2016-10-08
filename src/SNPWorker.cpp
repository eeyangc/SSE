#include <math.h>
#include <chrono>

#include "Cacher.h"
#include "SNPWorker.h"
#include "Conf.h"

namespace SNP
{
    
SNPWorker* SNPWorker::_snp_worker = nullptr;

#define eps 1e-5
#define sign(value) (value) > 0: 1: -1

int invserse2by2matrix(std::vector<std::vector<real_t> >& in, std::vector<std::vector<real_t> >& out)
{
    CHECK(in.size() == 2);
    CHECK(in[0].size() == 2);
    CHECK(in[1].size() == 2);
    CHECK(out.size() == 2);
    CHECK(out[0].size() == 2);
    CHECK(out[1].size() == 2);

    real_t det = in[0][0] * in[1][1] - in[0][1] * in[1][0];
//    if (fabs(det) < eps) {
//       det = sign(det) * eps;
//    }

    out[0][0] = in[1][1] / det;
    out[0][1] = -in[0][1] / det;
    out[1][0] = -in[1][0] / det;
    out[1][1] = in[0][0] / det;

    return 0;    
}

int chol2by2matrix(std::vector<std::vector<real_t> >& in, std::vector<std::vector<real_t> >& out)
{
    CHECK(in.size() == 2);
    CHECK(in[0].size() == 2);
    CHECK(in[1].size() == 2);
    CHECK(out.size() == 2);
    CHECK(out[0].size() == 2);
    CHECK(out[1].size() == 2);

    real_t det = in[0][0] * in[1][1] - in[0][1] * in[1][0];
    CHECK(det > 0);

    out[0][0] = sqrt(in[0][0]);
    out[0][1] = 0.0;
    out[1][0] = in[1][0] / out[0][0];
    out[1][1] = sqrt(in[1][1] - out[1][0] * out[1][0]);

    return 0;
}

int qr2by2matrix(std::vector<std::vector<real_t> >& in, std::vector<std::vector<real_t> >& R)
{
    //QR decomposition of 2-by-2 matrix, output R
    //Gram-Schmidt procedure
    //used in inverse-wishart distribution
    CHECK(in.size() == 2);
    CHECK(in[0].size() == 2);
    CHECK(in[1].size() == 2);
    CHECK(R.size() == 2);
    CHECK(R[0].size() == 2);
    CHECK(R[1].size() == 2);

    //real_t det = in[0][0] * in[1][1] - in[0][1] * in[1][0];
    //CHECK(det > 0);

    real_t norm1 = sqrt(in[0][0] * in[0][0] + in[1][0] * in[1][0]);
    real_t q00 = in[0][0] / norm1;
    real_t q10 = in[1][0] / norm1;

    R[0][0] = in[0][0] * q00 + in[1][0] * q10;
    R[0][1] = in[0][1] * q00 + in[1][1] * q10;

    real_t q01 = in[0][1] - R[0][1] * q00;
    real_t q11 = in[1][1] - R[0][1] * q10;

    real_t norm2 = sqrt(q01 * q01 + q11 * q11);
    q01 = q01 / norm2;
    q11 = q11 / norm2;

    R[1][1] = in[0][1] * q01 + in[1][1] * q11;
    R[1][0] = 0.0;
    return 0;
}

int times2by2matrix(std::vector<std::vector<real_t> >& in1, std::vector<std::vector<real_t> >& in2, std::vector<std::vector<real_t> >& out)
{
    //QR decomposition of 2-by-2 matrix, output R
    //Gram-Schmidt procedure
    //used in inverse-wishart distribution
    CHECK(in1.size() == 2);
    CHECK(in1[0].size() == 2);
    CHECK(in1[1].size() == 2);
    CHECK(in2.size() == 2);
    CHECK(in2[0].size() == 2);
    CHECK(in2[1].size() == 2);
    CHECK(out.size() == 2);
    CHECK(out[0].size() == 2);
    CHECK(out[1].size() == 2);

    out[0][0] = in1[0][0] * in2[0][0] + in1[0][1] * in2[1][0];
    out[0][1] = in1[0][0] * in2[0][1] + in1[0][1] * in2[1][1];
    out[1][0] = in1[1][0] * in2[0][0] + in1[1][1] * in2[1][0];
    out[1][1] = in1[1][0] * in2[0][1] + in1[1][1] * in2[1][1];

    return 0;
}

int gsl_ran_iwishart(
    const gsl_rng* r,
    std::vector<std::vector<real_t> >&in,
    real_t param,
    std::vector<std::vector<real_t> >&out)
{
    CHECK(in.size() == 2);
    CHECK(in[0].size() == 2);
    CHECK(in[1].size() == 2);
    CHECK(out.size() == 2);
    CHECK(out[0].size() == 2);
    CHECK(out[1].size() == 2);

    /*
    std::vector<std::vector<real_t> > d;
    d.resize(2);
    d[0].resize(2, 0);
    d[1].resize(2, 0);

    chol2by2matrix(in, d);

    std::vector<std::vector<real_t> > x;
    x.resize(2);
    x[0].resize(2, 0);
    x[1].resize(2, 0);

    x[0][0] = sqrt(gsl_ran_chisq(r, param));
    x[1][1] = sqrt(gsl_ran_chisq(r, param - 1));
    x[0][1] = gsl_ran_gaussian(r, 1);
    x[1][0] = 0.0;

    std::vector<std::vector<real_t> > R;
    R.resize(2);
    R[0].resize(2, 0);
    R[1].resize(2, 0);
    std::vector<std::vector<real_t> > invR;
    invR.resize(2);
    invR[0].resize(2, 0);
    invR[1].resize(2, 0);
    std::vector<std::vector<real_t> > T;
    T.resize(2);
    T[0].resize(2, 0);
    T[1].resize(2, 0);

    qr2by2matrix(x, R);
    invserse2by2matrix(R, invR);
    times2by2matrix(d, invR, T);

    out[0][0] = T[0][0] * T[0][0] + T[0][1] * T[0][1];
    out[0][1] = T[0][0] * T[1][0] + T[0][1] * T[1][1];
    out[1][0] = T[1][0] * T[0][0] + T[1][1] * T[0][1];
    out[1][1] = T[1][0] * T[1][0] + T[1][1] * T[1][1];
    */

    /*
    std::vector<std::vector<real_t> > Sigma_beta;
    Sigma_beta.resize(2);
    Sigma_beta[0].resize(2, 0);
    Sigma_beta[1].resize(2, 0);

    Sigma_beta[0][0] = T[0][0] * T[0][0] + T[0][1] * T[0][1];
    Sigma_beta[0][1] = T[0][0] * T[1][0] + T[0][1] * T[1][1];
    Sigma_beta[1][0] = T[1][0] * T[0][0] + T[1][1] * T[0][1];
    Sigma_beta[1][1] = T[1][0] * T[1][0] + T[1][1] * T[1][1];    

    invserse2by2matrix(Sigma_beta, out);
    */
    
    
    out[0][0] = in[0][0] / param;
    out[1][0] = in[1][0] / param;
    out[0][1] = in[0][1] / param;
    out[1][1] = in[1][1] / param;
    

    //invserse2by2matrix(in, out);
    

    return 0;
}
SNPWorker::SNPWorker(void)
{
    gsl_rng_env_setup();
    _r_global = gsl_rng_alloc(gsl_rng_default);
    //fprintf(stdout, "gsl seed %d\n", _r_global);
    gsl_rng_set(_r_global, time(NULL));
    _flag_single_trait = false;
    _flag_two_traits = false;
    _flag_cache = false;
}

SNPWorker::~SNPWorker(void)
{
    gsl_rng_free(_r_global);
    _beta.clear();

    if (_snp_worker != nullptr) {
        delete _snp_worker;
        _snp_worker = nullptr;
    }
}

int SNPWorker::init()
{
    const SNPConfs& confs = Conf::instance()->get_confs();
    CHECK_CONFS(confs, "num_thread");

    _num_threads = boost::lexical_cast<uint32_t>(confs.at("num_thread"));

    CHECK_CONFS(confs, "num_iteration");
    _num_iter = boost::lexical_cast<uint32_t>(confs.at("num_iteration"));

    CHECK_CONFS(confs, "num_burnin");
    _num_burnin = boost::lexical_cast<uint32_t>(confs.at("num_burnin"));

    fprintf(stdout, "begin to init worker\n");

    uint32_t dim = SNPData::instance()->dim();
    CHECK(_num_iter > 0);
    CHECK(dim > 0);

    CHECK_CONFS(confs, "enable_single_trait");
    uint32_t single_trait_flag = boost::lexical_cast<uint32_t>(confs.at("enable_single_trait"));
    _flag_single_trait = single_trait_flag == 1? true: false;

    CHECK_CONFS(confs, "enable_two_traits");
    uint32_t two_traits_flag = boost::lexical_cast<uint32_t>(confs.at("enable_two_traits"));
    _flag_two_traits = two_traits_flag == 1 ? true : false;
    if (SNPData::instance()->col_num() < 2) {
        fprintf(stdout, "co heritability is closed because col_num is smaller than 2.\n");
        _flag_two_traits = false;
    }

    CHECK_CONFS(confs, "enable_cache");
    uint32_t cache_flag = boost::lexical_cast<uint32_t>(confs.at("enable_cache"));
    _flag_cache = cache_flag == 1? true: false;
    if (_flag_cache == true) {
        _num_burnin = 0;
    }

    if (_flag_two_traits){
        fprintf(stdout, "single trait analysis is turned off because two-trait joint analysis is turned on.\n");
        _flag_single_trait = false;
    }
    
    if (_flag_single_trait) {
        _beta.resize(_num_iter);
        for (uint32_t i = 0; i < _num_iter; i++) {
            _beta[i].resize(dim);
        }
        _compute_buffer.resize(dim, 0);

        SNPCache::instance()->add_to_cache(_compute_buffer.data(), dim);
        SNPCache::instance()->add_to_cache(&_sigma_2_beta, 1);
    }

    if (_flag_two_traits) {
        _beta1.resize(_num_iter);
        _beta2.resize(_num_iter);
        for (uint32_t i = 0; i < _num_iter; i++) {
            _beta1[i].resize(dim);
            _beta2[i].resize(dim);
        }
        _co_compute_buffer.resize(2);
        _co_compute_buffer[0].resize(dim, 0);
        _co_compute_buffer[1].resize(dim, 0);

        SNPCache::instance()->add_to_cache(_co_compute_buffer[0].data(), dim);
        SNPCache::instance()->add_to_cache(_co_compute_buffer[1].data(), dim);
    
        // init co_heritability
        _co_lambda.clear();
        _co_lambda.resize(2);
        _co_lambda[0].resize(2, 0.0);
        _co_lambda[1].resize(2, 0.0);
        // init co lambda as an identiy matrix
        _co_lambda[0][0] = 1.0;
        _co_lambda[1][1] = 1.0;
        SNPCache::instance()->add_to_cache(_co_lambda[0].data(), 2);
        SNPCache::instance()->add_to_cache(_co_lambda[1].data(), 2);
        
       //CHECK_CONFS(confs, "wishart_nu");
       //_nu = boost::lexical_cast<real_t>(confs.at("wishart_nu"));
        _nu = 0;
    }

    return 0;
}

void SNPWorker::worker(uint32_t iter, worker_param_t& param)
{
    // init parameters
    param.heritability_beta = 0.0;
    param.co_heritability_beta[0][0] = 0.0;
    param.co_heritability_beta[0][1] = 0.0;
    param.co_heritability_beta[1][0] = 0.0;
    param.co_heritability_beta[1][1] = 0.0;

    SNPData* data = SNPData::instance();
    if (data == nullptr) {
        return;
    }

    while (true) {
        data_block_t* block = data->next();
#ifdef _DEBUG
        std::thread::id thread_id = std::this_thread::get_id();
        fprintf(stdout, "thread %x accesses block at %x\n", &thread_id, block);
#endif
        if (block == nullptr) {
            break;
        }

        // TODO: process data with one datablock
        
        // gibbs sampling with summary stat (single trait)
        gibbs_summary_stat(block, iter, param);

        // gibbs sampling for co-heritability (two traits)
        gibbs_coheritability(block, iter, param);
    }
}

void SNPWorker::gibbs_summary_stat(data_block_t* block, uint32_t iter, worker_param_t& param)
{
    /*
    * (1) perform calculation
    */
    if (block == nullptr) {
        return;
    }

    // const SNPConfs& confs = Conf::instance()->get_confs();
    if (_flag_single_trait == false) {
        return;
    }

    // obtain parameters to update
    int param_index = get_param_iter_index(iter);
    uint32_t beg_index = block->beg;
    uint32_t end_index = block->end;
    uint32_t p = end_index - beg_index;

    real_t sigma_beta2 = _sigma_2_beta;
    for (uint32_t k = 0; k < p; k++) {
        if (block->mask[k] == false) {
            continue;
        }

        //real_t sk = sqrt(block->s[k][0]);
        real_t sigma_k2 = 1.0 / (1.0 / (block->s[k][0]) + sigma_beta2);
        real_t corr_sum = 0.0;
        for (uint32_t j = 0; j < p; j++) {
            if (j == k) {
                continue;
            }

            if (block->mask[j] == false) {
                continue;
            }

            corr_sum += block->R[k][j] * _compute_buffer[beg_index + j] / sqrt(block->s[j][0]);
        }

        real_t mu_k = sigma_k2 * (block->beta[k][0] / (block->s[k][0]) - corr_sum / sqrt(block->s[k][0]));

        real_t beta_k = gsl_ran_gaussian(param.random_gen, sqrt(sigma_k2)) + mu_k;
        //real_t beta_k = mu_k;
      
        _compute_buffer[beg_index + k] = beta_k;
        param.heritability_beta += beta_k * beta_k;
        // push back the sampling data
        if (param_index != -1) {
            _beta[param_index][beg_index + k] = beta_k;
        }
    }
}

void SNPWorker::gibbs_coheritability(data_block_t* block, uint32_t iter, worker_param_t& param)
{
    if (block == nullptr) {
        return;
    }

    if (_flag_two_traits == false) {
        return;
    }

    // obtain parameters to update
    int param_index = get_param_iter_index(iter);
    uint32_t beg_index = block->beg;
    uint32_t end_index = block->end;
    uint32_t p = end_index - beg_index;

    std::vector<std::vector<real_t> > lambda_k;
    lambda_k.clear();
    lambda_k.resize(2);
    lambda_k[0].resize(2);
    lambda_k[1].resize(2);
    std::vector<std::vector<real_t> > inverse_lambda_k;
    inverse_lambda_k.clear();
    inverse_lambda_k.resize(2);
    inverse_lambda_k[0].resize(2);
    inverse_lambda_k[1].resize(2);
    
    for (uint32_t k = 0; k < p; k++) {
        if (block->mask[k] == false) {
            continue;
        }

        // compute lambda_k
        lambda_k[0][0] = 1.0 / block->s[k][0] + _co_lambda[0][0];
        lambda_k[0][1] = _co_lambda[0][1];
        lambda_k[1][0] = _co_lambda[1][0];
        lambda_k[1][1] = 1.0 / block->s[k][1] + _co_lambda[1][1];

        // get invserse
        invserse2by2matrix(lambda_k, inverse_lambda_k);
        real_t sigma_x = sqrt(inverse_lambda_k[0][0]);
        real_t sigma_y = sqrt(inverse_lambda_k[1][1]);
        real_t rho = inverse_lambda_k[0][1] / sigma_x / sigma_y;

        // get expected value
        real_t v1 = block->beta[k][0] / block->s[k][0];
        real_t v2 = block->beta[k][1] / block->s[k][1];

        real_t prod_sum1 = 0.0;
        real_t prod_sum2 = 0.0;
        for (uint32_t j = 0; j < p; j++) {
            if (j == k || block->mask[j] == false) {
                continue;
            }

            prod_sum1 += block->R[k][j] * _co_compute_buffer[0][beg_index + j] / sqrt(block->s[j][0]);
            prod_sum2 += block->R[k][j] * _co_compute_buffer[1][beg_index + j] / sqrt(block->s[j][1]);
        }

        v1 -= prod_sum1 / sqrt(block->s[k][0]);
        v2 -= prod_sum2 / sqrt(block->s[k][1]);

        real_t mu_k1 = inverse_lambda_k[0][0] * v1 + inverse_lambda_k[0][1] * v2;
        real_t mu_k2 = inverse_lambda_k[1][0] * v1 + inverse_lambda_k[1][1] * v2;

        // sampling bivariant gaussian
        double x, y;
        gsl_ran_bivariate_gaussian(param.random_gen, sigma_x, sigma_y, rho, &x, &y);
        x += mu_k1;
        y += mu_k2;

        _co_compute_buffer[0][beg_index + k] = x;
        _co_compute_buffer[1][beg_index + k] = y;

        // update param
        param.co_heritability_beta[0][0] += x * x;
        param.co_heritability_beta[0][1] += x * y;
        param.co_heritability_beta[1][0] += x * y;
        param.co_heritability_beta[1][1] += y * y;

        // push back the sampling data
        if (param_index != -1) {
            _beta1[param_index][beg_index + k] = x;
            _beta2[param_index][beg_index + k] = y;
        }  
    }
}

void SNPWorker::update_prior(uint32_t iter, std::vector<worker_param_t>& params)
{
    // update a and b
    SNPData* data = SNPData::instance();

    if (_flag_single_trait) {
        real_t a = 0;
        real_t b = 0;
        real_t a_hat = a + static_cast<real_t>(data->effective_dim()) ;
        real_t b_hat = 0;
    
        for (uint32_t i = 0; i < params.size(); i++) {
            b_hat += params[i].heritability_beta;
        }

        b_hat = b_hat + b;

        // sampling gamma
        //_sigma_2_beta = gsl_ran_gamma(_r_global, a_hat/2.0, 2.0/b_hat);
        _sigma_2_beta = a_hat / b_hat; //MLE
         
//#ifdef _DEBUG
        int param_index = get_param_iter_index(iter);
        if (param_index != -1) {
            fprintf(stdout, "a_hat=%e,b_hat=%e,sigma_2_beta=%e\n", a_hat, b_hat, 1.0 / _sigma_2_beta);
        }
//#endif
    }

    if (_flag_two_traits) {
        std::vector<std::vector<real_t> > lambda;
        lambda.resize(2);
        lambda[0].resize(2, 0);
        lambda[1].resize(2, 0);

        for (uint32_t i = 0; i < params.size(); i++) {
            lambda[0][0] += params[i].co_heritability_beta[0][0];
            lambda[0][1] += params[i].co_heritability_beta[0][1];
            lambda[1][0] += params[i].co_heritability_beta[1][0];
            lambda[1][1] += params[i].co_heritability_beta[1][1];
        }

        std::vector<std::vector<real_t> > co_Sigma;
        co_Sigma.resize(2);
        co_Sigma[0].resize(2, 0);
        co_Sigma[1].resize(2, 0);
        // sampling wishart
        gsl_ran_iwishart(_r_global, lambda, data->effective_dim() + _nu, co_Sigma);

        invserse2by2matrix(co_Sigma, _co_lambda);

        int param_index = get_param_iter_index(iter);
        if (param_index != -1) {
#ifdef _DEBUG
            fprintf(stdout, "Lambda[0][0]=%f,Lambda[0][1]=%f\n", _co_lambda[0][0], _co_lambda[0][1]);
            fprintf(stdout, "Lambda[1][0]=%f,Lambda[1][1]=%f\n", _co_lambda[1][0], _co_lambda[1][1]);
#endif
            real_t sigma1 = sqrt(co_Sigma[0][0]);
            real_t sigma2 = sqrt(co_Sigma[1][1]);
            real_t rho = co_Sigma[0][1] / sqrt(co_Sigma[0][0] * co_Sigma[1][1]);
#ifdef _DEBUG
            fprintf(stdout, "sigma_1 = %e, sigma_2 = %e, rho = %e\n", 
                sigma1, sigma2, rho);
#endif
            //fprintf(stdout, "Iteration: %u\n", iter);
            _rho[param_index] = rho;
            _sigma1[param_index] = sigma1;
            _sigma2[param_index] = sigma2;
        }
    }
}

int SNPWorker::run()
{
    SNPData* data = SNPData::instance();

    CHECK(data != nullptr);

    if (_num_threads <= 0) {
        return 0;
    }

    std::thread* threads = new (std::nothrow) std::thread[_num_threads];
    if (threads == nullptr) {
        fprintf(stderr, "Fail to allocate thread\n");
        return -1;
    }

    _sigma_2_beta = 1.0;

    fprintf(stdout, "burnin=[%u], Gibbs_samples = [%u],  Nthin = [%d]\n", 
        _num_burnin, _num_iter, NUM_ITER_PER_GROUP);
    fprintf(stdout, "total number of iters = burnin + (Gibbs_samples-1) * Nthin = [%u]\n", _num_burnin + (_num_iter-1) * NUM_ITER_PER_GROUP);

    _worker_params.clear();
    _worker_params.resize(_num_threads);

    for (uint32_t i = 0; i < _num_threads; i++) {
        // sleep for a while to let different thread uses different seeds
        // sleep(xxxx); sleep for 1 seconds
        fprintf(stdout, "Wait 0.01s for %u-th seed\n", i);
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        // use time as a generator
        _worker_params[i].random_gen = gsl_rng_alloc(gsl_rng_default);
        gsl_rng_set(_worker_params[i].random_gen, time(NULL));
    }

    _merger_params.clear();
    _merger_params.resize(_num_threads);

    _rho.clear();
    _rho.resize(_num_iter);
    _sigma1.clear();
    _sigma1.resize(_num_iter);
    _sigma2.clear();
    _sigma2.resize(_num_iter);

    set_error_no(0);

    if (_flag_cache == true) {
        SNPCache::instance()->get_cache_info();
    }

    for (uint32_t i = 0; i <= (_num_iter - 1) * NUM_ITER_PER_GROUP + _num_burnin; i++) {
        fprintf(stdout, "Iteration: %u\n", i);
        
        // iteration begin
        data->init_access();
        
        for (uint32_t j = 0; j < _num_threads; j++) {
            threads[j] = std::thread(&SNPWorker::worker, this, i, std::ref(_worker_params[j]));
        }

        // wait for thread to finish
        for (uint32_t j = 0; j < _num_threads; j++) {
            threads[j].join();
        }

        //CHECK_CMD(get_error_no() == 0, {delete[] threads; threads = nullptr;});

        // update prior
        update_prior(i, _worker_params);
    }

    // release random generator
    for (uint32_t i = 0; i < _num_threads; i++) {
        gsl_rng_free(_worker_params[i].random_gen);
    }

    fprintf(stdout, "gibbs sampling done (%d iters) \n", (_num_iter-1) * NUM_ITER_PER_GROUP + _num_burnin);
    // merge and dump the result
    data->init_access();
    for (uint32_t i = 0; i < _num_threads; i++) {
        threads[i] = std::thread(&SNPWorker::merger, this, std::ref(_merger_params[i]));
    }

    // wait for all thread to finish
    for (uint32_t i = 0; i < _num_threads; i++) {
        threads[i].join();
    }
    //CHECK_CMD(get_error_no() == 0, {delete[] threads; threads = nullptr;});

    fprintf(stdout, "start dumping result ... \n");
    CHECK_CMD(dump_result(_merger_params) == 0, {delete[] threads; threads = nullptr;});
    fprintf(stdout, "done! \n");
    // release threads
    delete[] threads;
    threads = nullptr;

    return 0;
}

void SNPWorker::merger(merger_param_t& param)
{
    // begin to perform data analysis
    std::vector<real_t>& heritability = param.heritability_heritability;
    std::vector<real_t>& heritability1 = param.heritability_heritability1;
    std::vector<real_t>& heritability2 = param.heritability_heritability2;
    // init params
    if (_flag_single_trait) {
        heritability.resize(_num_iter, 0);
    }

    if (_flag_two_traits) {
        // compute heritability based on two traits
        heritability1.clear();
        heritability1.resize(_num_iter);
        heritability2.clear();
        heritability2.resize(_num_iter);
    }

    SNPData* data = SNPData::instance();
    if (data == nullptr) {
        fprintf(stderr, "Cannot access data\n");
        return;
    }

    while (true) {
        data_block_t* block = data->next();

        if (block == nullptr) {
            break;
        }

        uint32_t beg = block->beg;
        uint32_t end = block->end;
        uint32_t p = end - beg;
        for (uint32_t k = 0; k < p; k++) {
            if (block->mask[k] == false) {
                continue;
            }

            if (_flag_single_trait) {
                real_t n_sk_2 = block->sample_size[k][0] * block->s[k][0];
                real_t beta = block->beta[k][0];
                real_t var_sum = n_sk_2 + beta * beta;

                for (uint32_t k_hat = 0; k_hat < p; k_hat++) {
                    if (block->mask[k_hat] == false) {
                        continue;
                    }

                    real_t n_sk_2_hat = block->sample_size[k_hat][0] * block->s[k_hat][0];
                    real_t beta_hat = block->beta[k_hat][0];
                    real_t var_sum_hat = n_sk_2_hat + beta_hat * beta_hat;

                    for (uint32_t iter = 0; iter < _num_iter; iter++) {
                        heritability[iter] += 
                            block->R[k][k_hat] * _beta[iter][beg + k] * _beta[iter][beg+k_hat] / sqrt(var_sum * var_sum_hat);
                    
                    }  // for iter
                }  // for k_hat
            }
            
            if (_flag_two_traits) {
                // heritability-based on two traits
                // trait1
                real_t n_sk_2_0 = block->sample_size[k][0] * block->s[k][0];
                real_t beta_0 = block->beta[k][0];
                real_t var_sum_0 = n_sk_2_0 + beta_0 * beta_0;
                // trait2
                real_t n_sk_2_1 = block->sample_size[k][1] * block->s[k][1];
                real_t beta_1 = block->beta[k][1];
                real_t var_sum_1 = n_sk_2_1 + beta_1 * beta_1;

                for (uint32_t k_hat = 0; k_hat < p; k_hat++) {
                    if (block->mask[k_hat] == false) {
                        continue;
                    }

                    real_t n_sk_2_hat_0 = block->sample_size[k_hat][0] * block->s[k_hat][0];
                    real_t beta_hat_0 = block->beta[k_hat][0];
                    real_t var_sum_hat_0 = n_sk_2_hat_0 + beta_hat_0 * beta_hat_0;

                    real_t n_sk_2_hat_1 = block->sample_size[k_hat][1] * block->s[k_hat][1];
                    real_t beta_hat_1 = block->beta[k_hat][1];
                    real_t var_sum_hat_1 = n_sk_2_hat_1 + beta_hat_1 * beta_hat_1;

                    for (uint32_t iter = 0; iter < _num_iter; iter++) {
                        heritability1[iter] +=
                            block->R[k][k_hat] * _beta1[iter][beg + k] * _beta1[iter][beg + k_hat] / sqrt(var_sum_0 * var_sum_hat_0);
                        heritability2[iter] +=
                            block->R[k][k_hat] * _beta2[iter][beg + k] * _beta2[iter][beg + k_hat] / sqrt(var_sum_1 * var_sum_hat_1);

                    }  // for iter
                }  // for k_hat
            }  // _flag_two_traits
        }  // for k
        
    }
}

int SNPWorker::dump_result(std::vector<merger_param_t>& params)
{
    // begin to perform data analysis
    const SNPConfs& confs = Conf::instance()->get_confs();
    CHECK_CONFS(confs, "output_beta_path");
    CHECK_CONFS(confs, "output_heritability_path");
    CHECK_CONFS(confs, "output_co_heritability_path");
    FILE* f_writer_beta = fopen(confs.at("output_beta_path").c_str(), "w");
    FILE* f_writer_her = fopen(confs.at("output_heritability_path").c_str(), "w");
    FILE* f_writer_co_her = fopen(confs.at("output_co_heritability_path").c_str(), "w");

    CHECK(f_writer_beta != nullptr);
    CHECK(f_writer_her != nullptr);
    CHECK(f_writer_co_her != nullptr);

    if (_flag_single_trait) {
        std::vector<real_t> heritability;
        heritability.resize(_num_iter, 0);

        for (uint32_t i = 0; i < params.size(); i++) {
            for (uint32_t j = 0; j < _num_iter; j++) {
                heritability[j] += params[i].heritability_heritability[j];
            }
        }

        CHECK_CONFS(confs, "output_beta");
        if (boost::lexical_cast<uint32_t>(confs.at("output_beta")) == 1) {
            for (uint32_t i = 0; i < _num_iter; i++) {
                //fprintf(f_writer_beta, "beta_%u = \n", i);
                for (uint32_t j = 0; j < SNPData::instance()->dim(); j++) {
                        fprintf(f_writer_beta, "%e\n", _beta[i][j]);                    
                }
                fprintf(f_writer_beta, "\n");
            }
        }

        CHECK_CONFS(confs, "output_heritability");
        if (boost::lexical_cast<uint32_t>(confs.at("output_heritability")) == 1) {
            for (uint32_t i = 0; i < heritability.size(); i++) {
                //fprintf(stdout, "heritability[%u]=%f\n", i, heritability[i]);
                fprintf(f_writer_her, "%e\n", heritability[i]);
            }
        }
    }

    if (_flag_two_traits) {
        // heritability based on two traits
        std::vector<real_t> heritability1;
        std::vector<real_t> heritability2;
        heritability1.resize(_num_iter, 0);
        heritability2.resize(_num_iter, 0);

        for (uint32_t i = 0; i < params.size(); i++) {
            for (uint32_t j = 0; j < _num_iter; j++) {
                heritability1[j] += params[i].heritability_heritability1[j];
                heritability2[j] += params[i].heritability_heritability2[j];
            }
        }

        CHECK_CONFS(confs, "output_beta");
        uint32_t output_co_beta = boost::lexical_cast<uint32_t>(confs.at("output_beta"));
        if (output_co_beta) {
            for (uint32_t i = 0; i < _num_iter; i++) {
                //fprintf(f_writer_beta, "co_beta=\n");
                for (uint32_t j = 0; j < SNPData::instance()->dim(); j++) {
                    fprintf(f_writer_beta, "%e,%e\n", _beta1[i][j], _beta2[i][j]);
                }
            }
        }

        CHECK_CONFS(confs, "output_co_heritability");
        uint32_t output_co_heritability = boost::lexical_cast<uint32_t>(confs.at("output_co_heritability"));
        //fprintf(f_writer_co_her, "--------------co-heritability--------------\n");
        if (output_co_heritability) {
            for (uint32_t i = 0; i < _num_iter; i++) {
                //fprintf(f_writer_co_her, "sigma1=%f,sigma2=%f,rho=%f\n", _sigma1[i], _sigma2[i], _rho[i]);
                fprintf(f_writer_co_her, "%e\n", _rho[i]);
            }
        }

        //fprintf(f_writer_co_her, "--------------heritability of two traits--------------\n");
        CHECK_CONFS(confs, "output_heritability");
        if (boost::lexical_cast<uint32_t>(confs.at("output_heritability")) == 1) {
            for (uint32_t i = 0; i < heritability1.size(); i++) {
                fprintf(f_writer_her, "%e,%e\n", heritability1[i], heritability2[i]);
            }
        }
    }

    fclose(f_writer_beta);
    fclose(f_writer_her);
    fclose(f_writer_co_her);

    return 0;
}

}
