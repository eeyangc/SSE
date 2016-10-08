#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "Conf.h"
#include "common.h"
#include "SNPData.h"
#include "SNPWorker.h"
#include "Cacher.h"

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        fprintf(stderr, "%s conf_file\n", argv[0]);
        return -1;
    }

	time_t start_t, end_t;

	
    // load configuration
    CHECK(SNP::Conf::instance()->load(argv[1]) == 0);

    // load data
    CHECK(SNP::SNPData::instance()->load() == 0);

    // begin to perform the task
    CHECK(SNP::SNPWorker::intance()->init() == 0);

    // init cache
    CHECK(SNP::SNPCache::instance()->init() == 0);

    start_t = time(NULL);
    CHECK(SNP::SNPWorker::intance()->run() == 0);
    end_t = time(NULL);

    fprintf(stdout, "Elapsed time:%ld s\n", (end_t - start_t));
    
    CHECK(SNP::SNPCache::instance()->dump_cache() == 0);

    return 0;
}
