#pragma once

#include <Rcpp.h>
#include <algorithm> // std::find
#include <vector>
#include <map>
#include <limits>
#include <string> // stoi
#include <cmath> // round
#include <math.h> // isnan

const float INFINITE_FLOAT =  std::numeric_limits <float>::max ();
const double INFINITE_DOUBLE =  std::numeric_limits <double>::max ();
const long int INFINITE_INT =  std::numeric_limits <long int>::max ();

typedef std::string vertex_id_t, edge_id_t;
typedef std::unordered_map <size_t,
    std::unordered_set <size_t> > int2ints_map_t;

struct edge_component
{
    // used only for edge sampling on graphs without component numbers
    edge_id_t edge;
    size_t component;
};

struct vertex_t
{
    private:
        std::unordered_set <vertex_id_t> in, out;

    public:
        void add_neighbour_in (vertex_id_t vert_id) { in.insert (vert_id); }
        void add_neighbour_out (vertex_id_t vert_id) { out.insert (vert_id); }
        size_t get_degree_in () { return in.size (); }
        size_t get_degree_out () { return out.size (); }

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
        std::set <edge_id_t> old_edges;

    public:
        double dist;
        double weight;
        double time;
        double timew;
        bool replaced_by_compact;

        vertex_id_t get_from_vertex () { return from; }
        vertex_id_t get_to_vertex () { return to; }
        edge_id_t getID () { return id; }
        std::set <edge_id_t> get_old_edges () { return old_edges; }

        edge_t (vertex_id_t from_id, vertex_id_t to_id,
                std::vector <double> weights_in, edge_id_t id_in,
                std::set <edge_id_t> replacement_edges)
        {
            replaced_by_compact = false;
            this -> to = to_id;
            this -> from = from_id;
            this -> dist = this -> weight = this -> time = weights_in [0];
            if (weights_in.size () > 1)
                this -> weight = weights_in [1];
            if (weights_in.size () > 2)
                this -> time = weights_in [2];
            if (weights_in.size () > 3)
                this -> timew = weights_in [3];
            this -> id = id_in;
            this -> old_edges.insert (replacement_edges.begin (),
                    replacement_edges.end ());
        }
};


typedef std::unordered_map <vertex_id_t, vertex_t> vertex_map_t;
typedef std::unordered_map <edge_id_t, edge_t> edge_map_t;
typedef std::unordered_map <vertex_id_t,
        std::unordered_set <edge_id_t> > vert2edge_map_t;

//----------------------------
//----- functions in graph.cpp
//----------------------------

namespace graph {

void add_to_v2e_map (vert2edge_map_t &vert2edge_map, const vertex_id_t vid,
        const edge_id_t eid);

void erase_from_v2e_map (vert2edge_map_t &vert2edge_map, const vertex_id_t vid,
        const edge_id_t eid);

bool graph_has_components (const Rcpp::DataFrame &graph);

bool graph_from_df (const Rcpp::DataFrame &gr, vertex_map_t &vm,
        edge_map_t &edge_map, vert2edge_map_t &vert2edge_map);

size_t identify_graph_components (vertex_map_t &v,
        std::unordered_map <vertex_id_t, size_t> &com);

} // end namespace graph

Rcpp::List rcpp_get_component_vector (const Rcpp::DataFrame &graph);

Rcpp::DataFrame rcpp_unique_rownames (Rcpp::DataFrame xyfrom, Rcpp::DataFrame xyto,
        const int precision);

//----------------------------
//----- functions in graph-sample.cpp
//----------------------------

namespace graph_sample {

edge_component sample_one_edge_no_comps (vertex_map_t &vertices,
        edge_map_t &edge_map);

edge_id_t sample_one_edge_with_comps (Rcpp::DataFrame graph,
        edge_map_t &edge_map);

vertex_id_t sample_one_vertex (Rcpp::DataFrame graph, vertex_map_t &vertices,
        edge_map_t &edge_map);

vertex_id_t select_random_vert (Rcpp::DataFrame graph,
        edge_map_t &edge_map, vertex_map_t &vertices);

} // end namespace graph_sample

Rcpp::StringVector rcpp_sample_graph (Rcpp::DataFrame graph,
        size_t nverts_to_sample);

//----------------------------
//----- functions in graph-contract.cpp
//----------------------------

namespace graph_contract {

void get_to_from (const edge_map_t &edge_map,
        const std::unordered_set <edge_id_t> &edges,
        const std::vector <vertex_id_t> &two_nbs,
        vertex_id_t &vt_from, vertex_id_t &vt_to,
        edge_id_t &edge_from_id, edge_id_t &edge_to_id);

void contract_one_edge (vert2edge_map_t &vert2edge_map,
        vertex_map_t &vertex_map, edge_map_t &edge_map,
        const std::unordered_set <edge_id_t> &edgelist,
        const vertex_id_t vtx_id, const vertex_id_t vt_from,
        const vertex_id_t vt_to,
        const edge_id_t edge_from_id, const edge_id_t edge_to_id,
        const edge_id_t new_edge_id,
        bool has_times);

bool same_hwy_type (const edge_map_t &edge_map, const edge_id_t &e1,
        const edge_id_t &e2);

void contract_graph (vertex_map_t &vertex_map, edge_map_t &edge_map,
        vert2edge_map_t &vert2edge_map,
        std::unordered_set <vertex_id_t> verts_to_keep,
        bool has_times);

} // end namespace graph_contract

Rcpp::List rcpp_contract_graph (const Rcpp::DataFrame &graph,
        Rcpp::Nullable <Rcpp::StringVector> &vertlist_in);

Rcpp::NumericVector rcpp_merge_flows (Rcpp::DataFrame graph);
