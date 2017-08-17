#pragma once

#include <Rcpp.h>
#include <algorithm> // std::find
#include <vector>
#include <map>
#include <limits>
#include <random>

const float INFINITE_FLOAT =  std::numeric_limits<float>::max ();
const int INFINITE_INT =  std::numeric_limits<int>::max ();

typedef std::string osm_id_t;
typedef int osm_edge_id_t;

struct osm_vertex_t
{
    private:
        std::unordered_set <osm_id_t> in, out;
        unsigned int component;
        double lat, lon;

    public:
        void add_neighbour_in (osm_id_t osm_id) { in.insert (osm_id); }
        void add_neighbour_out (osm_id_t osm_id) { out.insert (osm_id); }
        int get_degree_in () { return in.size (); }
        int get_degree_out () { return out.size (); }

        void set_component (unsigned int) { this -> component = component;  }

        void set_lat (double lat) { this -> lat = lat; }
        void set_lon (double lon) { this -> lon = lon; }
        double getLat () { return lat; }
        double getLon () { return lon; }

        std::unordered_set <osm_id_t> get_all_neighbours ()
        {
            std::unordered_set <osm_id_t> all_neighbours = in;
            all_neighbours.insert (out.begin (), out.end ());
            return all_neighbours;
        }

        std::unordered_set <osm_id_t> get_in_neighbours ()
        {
            return in;
        }

        std::unordered_set <osm_id_t> get_out_neighbours ()
        {
            return out;
        }

        void replace_neighbour (osm_id_t n_old, osm_id_t n_new)
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

struct osm_edge_t
{
    private:
        osm_id_t from, to;
        osm_edge_id_t id;
        std::set <unsigned int> old_edges;
        bool in_original_graph;

    public:
        float dist;
        float weight;
        unsigned int component;
        bool replaced_by_compact = false;
        std::string highway;

        osm_id_t get_from_vertex () { return from; }
        osm_id_t get_to_vertex () { return to; }
        osm_edge_id_t getID () { return id; }
        std::set <unsigned int> get_old_edges () { return old_edges; }
        bool in_original () { return in_original_graph; }

        osm_edge_t (osm_id_t from_id, osm_id_t to_id, float dist, float weight,
                   std::string highway, unsigned int id,
                   std::set <unsigned int> replacement_edges)
        {
            this -> to = to_id;
            this -> from = from_id;
            this -> dist = dist;
            this -> weight = weight;
            this -> highway = highway;
            this -> id = id;
            this -> old_edges.insert (replacement_edges.begin (),
                    replacement_edges.end ());
        }
};


typedef std::unordered_map <osm_id_t, osm_vertex_t> vertex_map_t;
typedef std::unordered_map <unsigned int, osm_edge_t> edge_map_t;
typedef std::unordered_map <osm_id_t, std::set <unsigned int>> vert2edge_map_t;

//----------------------------
//----- functions in graph.cpp
//----------------------------
void add_to_edge_map (vert2edge_map_t &vert2edge_map, osm_id_t vid,
        unsigned int eid);

void erase_from_edge_map (vert2edge_map_t &vert2edge_map, osm_id_t vid,
        unsigned int eid);

void graph_from_df (Rcpp::DataFrame gr, vertex_map_t &vm,
        edge_map_t &edge_map, vert2edge_map_t &vert2edge_map,
        bool is_spatial);

int identify_graph_components (vertex_map_t &v,
        std::unordered_map <osm_id_t, int> &com);

Rcpp::NumericVector rcpp_get_component_vector (Rcpp::DataFrame graph);

void add_components_to_graph (Rcpp::NumericVector &components, vertex_map_t &vm,
        edge_map_t &edge_map, vert2edge_map_t &vert2edge_map);

void contract_graph (vertex_map_t &vertex_map, edge_map_t &edge_map,
        vert2edge_map_t &vert2edge_map);

Rcpp::NumericVector rcpp_sample_graph (Rcpp::DataFrame graph,
        unsigned int nverts_to_sample, unsigned int e0, bool is_spatial);

Rcpp::List rcpp_contract_graph (Rcpp::DataFrame graph, bool is_spatial,
        bool quiet);

Rcpp::List rcpp_insert_vertices (Rcpp::DataFrame fullgraph,
        Rcpp::DataFrame compactgraph, std::vector <int> pts_to_insert,
        bool is_spatial);
