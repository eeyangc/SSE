#include "SNPData.h"
#include "Conf.h"

namespace SNP
{

SNPData* SNPData::_snp_data = nullptr;

SNPData::SNPData(void)
{
    _current_idx = 0;
    _col_num = 0;
    _effective_dim = 0;
}

SNPData::~SNPData(void)
{
    if (_snp_data != nullptr)
    {
        delete _snp_data;
        _snp_data = nullptr;
    }
}

int SNPData::load()
{
    const SNPConfs& confs = Conf::instance()->get_confs();

    // get data dimension
    CHECK_CONFS(confs, "dim");
    _dim = boost::lexical_cast<uint32_t>(confs.at("dim"));

    CHECK_CONFS(confs, "num_block");
    uint32_t num_blocks = boost::lexical_cast<uint32_t>(confs.at("num_block"));

    // init the size of data blocks
    _data_blocks.clear();
    _data_blocks.resize(num_blocks);

    // load Correlation information: from which to get block infomration
    CHECK(load_R(confs) == 0);
    CHECK(load_beta(confs) == 0);
    CHECK(load_s(confs) == 0);
    CHECK(load_sample_size(confs) == 0);

    // sort data according to block size
    sort();

    // get effective dim
    _effective_dim = 0;
    for (uint32_t i = 0; i < _data_blocks.size(); i++) {
        for (uint32_t j = 0; j < _data_blocks[i].mask.size(); j++) {
            if (_data_blocks[i].mask[j] == true) {
                _effective_dim += 1;
            }
        }
    }

    debug_print();

    return 0;
}

int SNPData::load_s(const SNPConfs& confs)
{
    CHECK_CONFS(confs, "f_s_path");     // s
    const char* f_path = confs.at("f_s_path").c_str();

    FILE* f_reader = fopen(f_path, "r");
    CHECK(f_reader != nullptr);

    char read_buffer[READ_BUFFER_LEN];
    uint32_t line_id = 0;
    uint32_t range_id = 0;
    uint32_t range_var = _block_range[range_id];
    std::vector<std::vector<real_t> >* sel_block = &_data_blocks[range_id].s;
    std::vector<bool>* sel_mask = &_data_blocks[range_id].mask;
    sel_block->clear();
    uint32_t mask_id = 0;

    fprintf(stdout, "begin to load the s matrix\n");
    while (fgets(read_buffer, READ_BUFFER_LEN, f_reader) != nullptr) {
        if (read_buffer == NULL || read_buffer[0] == '\0' || read_buffer[0] == '\n' ||
                read_buffer[0] == '\r') {
                    continue;
        }

        std::vector<std::string> res_list;
        std::string str_read = read_buffer;
        std::string tmp_str = boost::trim_copy(str_read);
        boost::split(res_list, tmp_str, boost::is_any_of(" "));

        if (line_id >= range_var) {
            range_id++;
            range_var = _block_range[range_id];
            sel_block = &_data_blocks[range_id].s;
            sel_mask = &_data_blocks[range_id].mask;
            sel_block->clear();
            mask_id = 0;
        }

        CHECK_CMD(res_list.size() <= 2, fclose(f_reader));

        bool mask_value = true;
        std::vector<real_t> row_var;
        row_var.clear();
        for (uint32_t idx = 0; idx < res_list.size(); idx++) {
            //row_var.push_back(boost::lexical_cast<real_t>(res_list[idx]));
            real_t value = 0.0;
            if (strcmp(res_list[idx].c_str(), "NA") == 0) {
                mask_value = false;
            }
            else {
                value = boost::lexical_cast<real_t>(res_list[idx]);
            }

            row_var.push_back(value);
        }
        sel_block->push_back(row_var);
        CHECK(mask_id < sel_mask->size());
        mask_value &=  sel_mask->operator[](mask_id);
        sel_mask->operator[](mask_id) = mask_value;
        mask_id++;
        line_id++;
    }

    fclose(f_reader);

    CHECK(_dim == line_id);

    return 0;
}

int SNPData::load_R(const SNPConfs& confs)
{
    CHECK_CONFS(confs, "f_corr_path");  // R
    uint32_t num_blocks = _data_blocks.size();

    char file_name[STR_BUFFER_LEN];
    char read_buffer[READ_BUFFER_LEN];

    uint32_t p = 0;
    _block_range.clear();

    for (uint32_t i = 0; i < num_blocks; i++) {
        int ret = snprintf(file_name, STR_BUFFER_LEN, "%s/block_%u",
            confs.at("f_corr_path").c_str(), i);

        if (ret < 0 || ret >= STR_BUFFER_LEN) {
            fprintf(stderr, "Fail to get file full path\n");
            return -1;
        }

        FILE* f_reader = fopen(file_name, "r");
        fprintf(stdout, "Begin to load correlation matrix: %u\n", i);
        CHECK(f_reader != nullptr);
        std::vector<std::vector<real_t> >& sel_R = _data_blocks[i].R;
        std::vector<bool>& sel_mask = _data_blocks[i].mask;

        uint32_t read_line = 0;
        uint32_t line_dim = 0;
        while(fgets(read_buffer, READ_BUFFER_LEN, f_reader) != nullptr) {
            if (read_buffer == NULL || read_buffer[0] == '\0' || read_buffer[0] == '\n' ||
                read_buffer[0] == '\r') {
                    continue;
            }

            if (read_line == 0) {
                std::vector<std::string> res_str;
                std::string read_str = read_buffer;
                std::string tmp_str = boost::trim_copy(read_str);
                boost::split(res_str, tmp_str, boost::is_any_of(" "));

                line_dim = res_str.size();
                _data_blocks[i].beg = p;
                p += line_dim;
                _data_blocks[i].end = p;
                _block_range.push_back(p);

                sel_R.clear();
                sel_R.resize(line_dim);
                sel_mask.clear();
                sel_mask.resize(line_dim, true);

                for (uint32_t idx = 0; idx < line_dim; idx++) {
                    sel_R[idx].resize(line_dim);
                }
            }

            CHECK_CONFS(confs, "shrinkage");
            _shrinkage = boost::lexical_cast<real_t>(confs.at("shrinkage"));
            CHECK(_shrinkage >= 0);
            CHECK(_shrinkage <= 1);

            // begin to process the file
            CHECK(line_dim > 0);

            char* pch = nullptr;
            char* s_ptr = nullptr;
            pch = strtok_r(read_buffer, " \n\r", &s_ptr);
            uint32_t col_idx = 0;
            while (pch != nullptr) {
                if (read_line == col_idx){
                    sel_R[read_line][col_idx] = atof(pch);
                }
                else
                {
                    sel_R[read_line][col_idx] = atof(pch) * _shrinkage;
                }
                col_idx++;

                pch = strtok_r(nullptr, " \n\r", &s_ptr);
            }
            
            read_line++;
        }
        
        CHECK_CMD(line_dim == read_line, fclose(f_reader));

        fclose(f_reader);
    }

    CHECK(p == _dim);

    return 0;
}

int SNPData::load_beta(const SNPConfs& confs)
{
    CHECK_CONFS(confs, "f_beta_path");  // hat(beta)

    const char* f_path = confs.at("f_beta_path").c_str();

    FILE* f_reader = fopen(f_path, "r");
    CHECK(f_reader != nullptr);

    char read_buffer[READ_BUFFER_LEN];
    uint32_t line_id = 0;
    uint32_t range_id = 0;
    uint32_t range_var = _block_range[range_id];
    std::vector<std::vector<real_t> >* sel_block = &_data_blocks[range_id].beta;
    std::vector<bool>* sel_mask = &_data_blocks[range_id].mask;
    sel_block->clear();
    uint32_t mask_id = 0;

    fprintf(stdout, "begin to load the beta matrix\n");
    while (fgets(read_buffer, READ_BUFFER_LEN, f_reader) != nullptr) {
        if (read_buffer == NULL || read_buffer[0] == '\0' || read_buffer[0] == '\n' ||
                read_buffer[0] == '\r') {
                    continue;
        }

        std::string tmp_str = read_buffer;
        std::string read_data = boost::trim_copy(tmp_str);

        if (line_id >= range_var) {
            range_id++;
            range_var = _block_range[range_id];
            sel_block = &_data_blocks[range_id].beta;
            sel_mask = &_data_blocks[range_id].mask;
            sel_block->clear();
            mask_id = 0;
        }

        std::vector<std::string> res_list;

        boost::split(res_list, read_data, boost::is_any_of(" "));
        uint32_t col_num = res_list.size();
        CHECK_CMD(col_num <= 2, fclose(f_reader));

        if (_col_num == 0 || col_num < _col_num) {
            _col_num = col_num;
        }

        bool mask_value = true;
        std::vector<real_t> row_var;
        row_var.clear();
        for (uint32_t idx = 0; idx < col_num; idx++) {
            real_t value = 0.0;
            if (strcmp(res_list[idx].c_str(), "NA") == 0) {
                mask_value = false;
            } else {
                value = boost::lexical_cast<real_t>(res_list[idx]);
            }

            row_var.push_back(value);
        }

        sel_block->push_back(row_var);
        CHECK(mask_id < sel_mask->size());
        mask_value &= sel_mask->operator[](mask_id);
        sel_mask->operator[](mask_id) = mask_value;
        mask_id ++;
        line_id++;
    }

    fclose(f_reader);

    CHECK(line_id == _dim);

    return 0;
}

int SNPData::load_sample_size(const SNPConfs& confs)
{
    CHECK_CONFS(confs, "sample_size_path");     // s
    const char* f_path = confs.at("sample_size_path").c_str();

    FILE* f_reader = fopen(f_path, "r");
    CHECK(f_reader != nullptr);

    char read_buffer[READ_BUFFER_LEN];
    uint32_t line_id = 0;
    uint32_t range_id = 0;
    uint32_t range_var = _block_range[range_id];
    std::vector<std::vector<real_t> >* sel_block = &_data_blocks[range_id].sample_size;
    std::vector<bool>* sel_mask = &_data_blocks[range_id].mask;
    sel_block->clear();
    uint32_t mask_id = 0;

    fprintf(stdout, "begin to load the sample size matrix\n");
    while (fgets(read_buffer, READ_BUFFER_LEN, f_reader) != nullptr) {
        if (read_buffer == NULL || read_buffer[0] == '\0' || read_buffer[0] == '\n' ||
                read_buffer[0] == '\r') {
                    continue;
        }

        std::vector<std::string> res_list;
        std::string str_read = read_buffer;
        std::string tmp_str = boost::trim_copy(str_read);
        boost::split(res_list, tmp_str, boost::is_any_of(" "));

        if (line_id >= range_var) {
            range_id++;
            range_var = _block_range[range_id];
            sel_block = &_data_blocks[range_id].sample_size;
            sel_mask = &_data_blocks[range_id].mask;
            sel_block->clear();
            mask_id = 0;
        }

        CHECK_CMD(res_list.size() <= 2, fclose(f_reader));
        bool mask_value = true;
        std::vector<real_t> row_var;
        row_var.clear();
        for (uint32_t idx = 0; idx < res_list.size(); idx++) {
            //row_var.push_back(boost::lexical_cast<real_t>(res_list[idx]));
            real_t value = 0.0;
            if (strcmp(res_list[idx].c_str(), "NA") == 0) {
                mask_value = false;
            }
            else {
                value = boost::lexical_cast<real_t>(res_list[idx]);
            }

            row_var.push_back(value);
        }
        sel_block->push_back(row_var);
        CHECK(mask_id < sel_mask->size());
        mask_value &= sel_mask->operator[](mask_id);
        sel_mask->operator[](mask_id) = mask_value;
        mask_id++;
        line_id++;

    }

    fclose(f_reader);

    CHECK(_dim == line_id);

    return 0;
}

void SNPData::debug_print()
{
#ifdef _DEBUG
    uint32_t block_num = _data_blocks.size();
    fprintf(stdout, "the number of block is: %u\n", block_num);
    for (uint32_t i = 0; i < block_num; i++) {
        fprintf(stdout, "-------------------\n");
        fprintf(stdout, "beg=%u, end=%d\n", _data_blocks[i].beg, _data_blocks[i].end);
        // print s
        fprintf(stdout, "S_matrix = \n");
        std::vector<std::vector<real_t> >& sel_s = _data_blocks[i].s;
        for (uint32_t j = 0; j < sel_s.size(); j++) {
            for (uint32_t k = 0; k < sel_s[j].size(); k++) {
                if (k == 0) {
                    fprintf(stdout, "%f", sel_s[j][k]);
                } else {
                    fprintf(stdout, " %f", sel_s[j][k]);
                }
            }
            fprintf(stdout, "\n");
        }

        // print sample size
        fprintf(stdout, "Sample_size_matrix = \n");
        std::vector<std::vector<real_t> >& sel_sample_size = _data_blocks[i].sample_size;
        for (uint32_t j = 0; j < sel_sample_size.size(); j++) {
            for (uint32_t k = 0; k < sel_sample_size[j].size(); k++) {
                if (k == 0) {
                    fprintf(stdout, "%f", sel_sample_size[j][k]);
                } else {
                    fprintf(stdout, " %f", sel_sample_size[j][k]);
                }
            }
            fprintf(stdout, "\n");
        }

        // print betahat
        fprintf(stdout, "betahat_matrix = \n");
        std::vector<std::vector<real_t> >& sel_beta = _data_blocks[i].beta;
        for (uint32_t j = 0; j < sel_beta.size(); j++) {
            for (uint32_t k = 0; k < sel_beta[j].size(); k++) {
                if (k == 0) {
                    fprintf(stdout, "%f", sel_beta[j][k]);
                } else {
                    fprintf(stdout, " %f", sel_beta[j][k]);
                }
            }
            fprintf(stdout, "\n");
        }

        // print R
        /*
        fprintf(stdout, "R = \n");
        std::vector<std::vector<real_t> >& sel_R = _data_blocks[i].R;
        for (uint32_t j = 0; j < sel_R.size(); j++) {
            for (uint32_t k = 0; k < sel_R[j].size(); k++) {
                if (k == 0) {
                    fprintf(stdout, "%f", sel_R[j][k]);
                } else {
                    fprintf(stdout, " %f", sel_R[j][k]);
                }
            }
            fprintf(stdout, "\n");
        }
        */
        // print mask
        fprintf(stdout, "mask = \n");
        std::vector<bool>& mask = _data_blocks[i].mask;
        for (uint32_t j = 0; j < mask.size(); j++) {
            fprintf(stdout, "%d\n", mask[j] == true?1:0);
        }
    }
#endif
}

}
