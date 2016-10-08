#include "Cacher.h"
#include "Conf.h"

namespace SNP
{
SNPCache* SNPCache::_cacher = nullptr;

SNPCache::~SNPCache()
{
    if (_cacher != nullptr) {
        delete _cacher;
        _cacher = nullptr;
    }
}

int SNPCache::load_info()
{
    std::vector<real_t> tmp_buffer;
    const SNPConfs& confs = Conf::instance()->get_confs();
    CHECK_CONFS(confs, "cache_file");
    const char* f_name = confs.at("cache_file").c_str();
    std::fstream f_reader(f_name, std::ios::in | std::ios::binary);
    if (f_reader.good() == false) {
        return -1;
     }

    for (uint32_t k = 0; k < _data.size(); k++) {
        for (uint32_t i = 0; i < _data[k].size; i++) {
            real_t value;
            f_reader.read((char*)&value, sizeof(real_t));
            CHECK(f_reader.good() == true);
            tmp_buffer.push_back(value);
        }
    }

    f_reader.close();
    
    // use buffer content
    uint32_t offset = 0;
    for (uint32_t k = 0; k < _data.size(); k++) {
        for (uint32_t i = 0; i < _data[k].size; i++) {
            _data[k].data[i] = tmp_buffer[offset];
            offset += 1;
        }
    }

    return 0;
}

int SNPCache::dump_info()
{
    const SNPConfs& confs = Conf::instance()->get_confs();
    CHECK_CONFS(confs, "cache_file");
    const char* f_name = confs.at("cache_file").c_str();
    std::fstream f_writer(f_name, std::ios::out | std::ios::binary);
    CHECK(f_writer.good() == true);

    for (uint32_t k = 0; k < _data.size(); k++) {
        cache_type_t& data = _data[k];
        for (uint32_t i = 0; i < data.size; i++) {
            f_writer.write((const char*)&data.data[i], sizeof(real_t));
        }
    } 

    f_writer.close();

    return 0;
}

}
