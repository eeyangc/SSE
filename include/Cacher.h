#ifndef _SNP_CACHE_H
#define _SNP_CACHE_H

#include <fstream>

#include "common.h"

namespace SNP {

struct cache_type_t
{
    real_t* data;
    uint32_t size;
};

class SNPCache
{
public:
    SNPCache() {}
    ~SNPCache();

    static SNPCache* instance()
    {
        if (_cacher == nullptr) {
            _cacher = new (std::nothrow) SNPCache();
        }
        
        return _cacher;
    }

    int init()
    {
        return 0;
    }

    int add_to_cache(real_t* data, uint32_t size)
    {
        cache_type_t item;
        item.data = data;
        item.size = size;

        _data.push_back(item);

        return 0;
    }

    int get_cache_info()
    {
        return load_info();
    }

    int dump_cache() 
    {
        return dump_info();
    }

private:
    int load_info();
    int dump_info();
private:
    std::vector<cache_type_t> _data;

    static SNPCache* _cacher;
};

}

#endif
