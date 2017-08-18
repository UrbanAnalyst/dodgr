#pragma once

#include <Rcpp.h>
#include <algorithm> // std::find
#include <vector>
#include <map>
#include <limits>
#include <random>

const float INFINITE_FLOAT =  std::numeric_limits<float>::max ();
const int INFINITE_INT =  std::numeric_limits<int>::max ();

typedef std::string vertex_id_t;
typedef int edge_id_t;

struct vertex_t
{
    private:
        std::unordered_set <vertex_id_t> in, out;

    public:
        void add_neighbour_in (vertex_id_t vert_id) { in.insert (vert_id); }
        void add_neighbour_out (vertex_id_t vert_id) { out.insert (vert_id); }
        int get_degree_in () { return in.size (); }
        int get_degree_out () { return out.size (); }

        std::unordered_set <vertex_id_t> get_all_neighbours ()
        {
            std::unordered_set <vertex_id_t> all_neighbours = in;
            all_neighbours.insert (out.begin (), out.end ());
            return all_neighbours;
        }

        std::unordered_set <vertex_id_t> get_in_neighbours ()
        {
            return in;
        }

        std::unordered_set <vertex_id_t> get_out_neighbours ()
        {
            return out;
        }

        void replace_neighbour (vertex_id_t n_old, vertex_id_t n_new)
        {
            if (in.find (n_old) != in.end ())
            {
                in.erase (n_old);
                in.insert (n_new);
            }
            if (out.find (n_old) != out.end ())
            {
                out.erase (n_old);
                out.insert (n_new);
            }
        }

        // in can equal out, so get_all_neighbours is vital here:
        bool is_intermediate_single ()
        {
            return (in.size () == 1 && out.size () == 1 &&
                    get_all_neighbours ().size () == 2);
        }

        bool is_intermediate_double ()
        {
            return (in.size () == 2 && out.size () == 2 &&
                    get_all_neighbours ().size () == 2);
        }
};

struct edge_t
{
    private:
        vertex_id_t from, to;
        edge_id_t id;
        std::set <unsigned int> old_edges;
        bool in_original_graph;

    public:
        float dist;
        float weight;
        bool replaced_by_compact = false;

        vertex_id_t get_from_vertex () { return from; }
        vertex_id_t get_to_vertex () { return to; }
        edge_id_t getID () { return id; }
        std::set <unsigned int> get_old_edges () { return old_edges; }
        bool in_original () { return in_original_graph; }

        edge_t (vertex_id_t from_id, vertex_id_t to_id, float dist, float weight,
                   unsigned int id, std::set <unsigned int> replacement_edges)
        {
            this -> to = to_id;
            this -> from = from_id;
            this -> dist = dist;
            this -> weight = weight;
            this -> id = id;
            this -> old_edges.insert (replacement_edges.begin (),
                    replacement_edges.end ());
        }
};


typedef std::unordered_map <vertex_id_t, vertex_t> vertex_map_t;
typedef std::unordered_map <unsigned int, edge_t> edge_map_t;
typedef std::unordered_map <vertex_id_t, std::set <unsigned int>> vert2edge_map_t;

//----------------------------
//----- functions in graph.cpp
//----------------------------
void add_to_edge_map (vert2edge_map_t &vert2edge_map, vertex_id_t vid,
        unsigned int eid);

void erase_from_edge_map (vert2edge_map_t &vert2edge_map, vertex_id_t vid,
        unsigned int eid);

void graph_from_df (Rcpp::DataFrame gr, vertex_map_t &vm,
        edge_map_t &edge_map, vert2edge_map_t &vert2edge_map,
        bool is_spatial);

int identify_graph_components (vertex_map_t &v,
        std::unordered_map <vertex_id_t, int> &com);

Rcpp::NumericVector rcpp_get_component_vector (Rcpp::DataFrame graph,
        bool is_spatial);

//----------------------------
//----- functions in graph-sample.cpp
//----------------------------
std::vector <unsigned int>  sample_one_edge_no_comps (vertex_map_t &vertices,
        edge_map_t &edge_map);

unsigned int sample_one_edge_with_comps (Rcpp::DataFrame graph);

bool graph_has_components (Rcpp::DataFrame graph);

vertex_id_t sample_one_vertex (Rcpp::DataFrame graph, vertex_map_t &vertices,
        edge_map_t &edge_map);

Rcpp::NumericVector rcpp_sample_graph (Rcpp::DataFrame graph,
        unsigned int nverts_to_sample, bool is_spatial);

//----------------------------
//----- functions in graph-contract.cpp
//----------------------------
void contract_graph (vertex_map_t &vertex_map, edge_map_t &edge_map,
        vert2edge_map_t &vert2edge_map);

Rcpp::List rcpp_contract_graph (Rcpp::DataFrame graph, bool is_spatial,
        bool quiet);

Rcpp::List rcpp_insert_vertices (Rcpp::DataFrame fullgraph,
        Rcpp::DataFrame compactgraph, std::vector <int> pts_to_insert,
        bool is_spatial);
