
#include <Rcpp.h>

namespace sc {
    std::string random_id (size_t len);
}

Rcpp::CharacterVector rcpp_gen_hash (const int n, const size_t hash_len);
