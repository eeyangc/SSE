 /* audit -- definition of audit_context structure and supporting types
  *
  * Copyright 2003-2004 Red Hat, Inc.
  * Copyright 2005 Hewlett-Packard Development Company, L.P.
  * Copyright 2005 IBM Corporation
  *
  * This program is free software; you can redistribute it and/or modify
  * it under the terms of the GNU General Public License as published by
  * the Free Software Foundation; either version 2 of the License, or
  * (at your option) any later version.
  *
  * This program is distributed in the hope that it will be useful,
  * but WITHOUT ANY WARRANTY; without even the implied warranty of
  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  * GNU General Public License for more details.
  *
  * You should have received a copy of the GNU General Public License
  * along with this program; if not, write to the Free Software
  * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
  */

#ifndef _SNP_COMMON_H
#define _SNP_COMMON_H

#include <stdio.h>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <algorithm>
#include <functional> 
#include <stdint.h>
#include <thread>
#include <mutex>
#include <boost/format.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

namespace SNP
{
#define READ_BUFFER_LEN 102400
#define STR_BUFFER_LEN 512
#define NUM_ITER_PER_GROUP 20
#define real_t double


typedef std::map<std::string, std::string> SNPConfs;

#define CHECK_CONFS(confs, item) do\
{ \
    SNPConfs::const_iterator it = confs.find(item); \
    if (it == confs.end()) { \
        fprintf(stderr, "Cannot find [%s] in your configure at %d of %s\n", \
            item, __LINE__, __FILE__); \
        return -1; \
    } \
} while (0)

#define CHECK_CONFS_VOID(confs, item, cmd) do\
{ \
    SNPConfs::const_iterator it = confs.find(item); \
    if (it == confs.end()) { \
        fprintf(stderr, "Cannot find [%s] in your configure at %d of %s\n", \
            item, __LINE__, __FILE__); \
        cmd; \
        return; \
    } \
} while (0)

#define CHECK(expression) do {\
    if (!(expression)) {fprintf(stderr, \
    "Fail to check expression at line %d of %s\n", __LINE__, __FILE__); return -1;}\
} while (0)

#define CHECK_CMD(expression, cmd) do {\
    if (!(expression)) {fprintf(stderr, \
    "Fail to check expression at line %d of %s\n", __LINE__, __FILE__); cmd; return -1;}\
} while (0)

}

#ifdef _WIN32
#define snprintf _snprintf
#define strtok_r strtok_s
#endif

#endif