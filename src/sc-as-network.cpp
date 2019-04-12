#include "sc-as-network.h"

// Function to generate IDs for the edges in each way
std::string sc::random_id (size_t len) {
    auto randchar = []() -> char
    {
        const char charset[] = \
            "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
        const size_t max_index = (sizeof(charset) - 1);
        //return charset [ rand() % max_index ];
        size_t i = static_cast <size_t> (floor (Rcpp::runif (1) [0] * max_index));
        return charset [i];
    };
    std::string str (len, 0);
    std::generate_n (str.begin(), len, randchar);
    return str;
}

//' rcpp_gen_hash
//'
//' Efficient generation of long sequences of hash keys
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::CharacterVector rcpp_gen_hash (const int n, const int hash_len)
{
    Rcpp::CharacterVector res (n);
    for (int i = 0; i < n; i++)
        res [i] = sc::random_id (hash_len);
    return res;
}
