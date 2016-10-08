#include "Conf.h"

namespace SNP
{

Conf* Conf::_p_conf = nullptr;

Conf::Conf(void)
{
}

Conf::~Conf(void)
{
    if (_p_conf != nullptr) {
        delete _p_conf;
        _p_conf = nullptr;
    }
}

int Conf::load(const char* f_name)
{
    if (f_name == nullptr || f_name[0] == '\0') {
        fprintf(stderr, "Invalid conf file name\n");
        return -1;
    }

    FILE* f_reader = fopen(f_name, "r");
    if (f_reader == nullptr) {
        fprintf(stderr, "Fail to open the conf file: %s\n", f_name);
        return -1;
    }

    char* read_buffer = new (std::nothrow) char[READ_BUFFER_LEN];
    uint32_t line_count = 0;
    while (fgets(read_buffer, READ_BUFFER_LEN, f_reader) != NULL) {
        line_count++;

        if (read_buffer == NULL || read_buffer[0] == '\0' || read_buffer[0] == '#'
            || read_buffer[0] == '\n' || read_buffer[0] == '\r') {
            continue;
        }

        std::vector<std::string> res_str;
        boost::split(res_str, read_buffer, boost::is_any_of("="));

        CHECK(res_str.size() >= 2);

        std::string key = boost::trim_copy(res_str[0]);
        std::string value = boost::trim_copy(res_str[1]);

#ifdef _DEBUG
        fprintf(stdout, "load key=[%s] with value=[%s]\n", key.c_str(), value.c_str());
#endif
        _confs[key] = value;
    }

    fclose(f_reader);
    return 0;
}

}
