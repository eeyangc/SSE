#ifndef _SNP_CONF_H
#define _SNP_CONF_H

#include "common.h"

namespace SNP
{

class Conf
{
public:
    static Conf* instance()
    {
        if (_p_conf == nullptr) {
            _p_conf = new (std::nothrow) Conf();
        }

        return _p_conf;
    }

public:
    Conf(void);
    ~Conf(void);

    int load(const char* f_name);

    const SNPConfs& get_confs() const {
        return _confs;
    }

private:
    SNPConfs _confs;

    static Conf* _p_conf;
};

}

#endif
