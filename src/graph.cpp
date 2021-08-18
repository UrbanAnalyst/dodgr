#include "graph.h"

void graph::add_to_v2e_map (vert2edge_map_t &vert2edge_map, const vertex_id_t vid,
        const edge_id_t eid)
{
    std::unordered_set <edge_id_t> edge_ids;
    if (vert2edge_map.find (vid) == vert2edge_map.end ())
    {
        edge_ids.emplace (eid);
        vert2edge_map.emplace (vid, edge_ids);
    } else
    {
        edge_ids = vert2edge_map [vid];
        edge_ids.insert (eid);
        vert2edge_map [vid] = edge_ids;
    }
}

void graph::erase_from_v2e_map (vert2edge_map_t &vert2edge_map,
        const vertex_id_t vid, const edge_id_t eid)
{
    std::unordered_set <edge_id_t> edge_ids = vert2edge_map [vid];
    if (edge_ids.find (eid) != edge_ids.end ())
    {
        edge_ids.erase (eid);
        vert2edge_map [vid] = edge_ids;
    }
}

//' graph_has_components
//'
//' Does a graph have a vector of connected component IDs? Only used in
//' \code{sample_one_vertex}
//' @noRd
bool graph::graph_has_components (const Rcpp::DataFrame &graph)
{
    Rcpp::CharacterVector graph_names = graph.attr ("names");
    for (auto n: graph_names)
        if (n == "component")
            return true;

    return false;
}


//' @name graph_from_df
//'
//' Convert a standard graph data.frame into an object of class graph. Graphs 
//' are standardised with the function \code{dodgr_convert_graph()$graph}, and
//' contain only the four columns [from, to, d, w]
//'
//' @noRd
bool graph::graph_from_df (const Rcpp::DataFrame &gr, vertex_map_t &vm,
        edge_map_t &edge_map, vert2edge_map_t &vert2edge_map)
{
    Rcpp::StringVector edge_id = gr ["edge_id"];
    Rcpp::StringVector from = gr ["from"];
    Rcpp::StringVector to = gr ["to"];
    Rcpp::NumericVector dist = gr ["d"];
    Rcpp::NumericVector weight = gr ["d_weighted"];
    Rcpp::StringVector colnames = gr.attr ("names");
    bool has_times = false;
    Rcpp::NumericVector time, timew;
    if (gr.ncol () == 7)
    {
        has_times = true;
        time = gr ["time"];
        timew = gr ["time_weighted"];
    }

    std::set <edge_id_t> replacement_edges; // all empty here

    for (R_xlen_t i = 0; i < to.length (); i ++)
    {
        Rcpp::checkUserInterrupt ();
        vertex_id_t from_id = std::string (from [i]);
        vertex_id_t to_id = std::string (to [i]);

        if (vm.find (from_id) == vm.end ())
        {
            vertex_t fromV = vertex_t ();
            vm.emplace (from_id, fromV);
        }
        vertex_t from_vtx = vm.at (from_id);
        from_vtx.add_neighbour_out (to_id);
        vm [from_id] = from_vtx;

        if (vm.find (to_id) == vm.end ())
        {
            vertex_t toV = vertex_t ();
            vm.emplace (to_id, toV);
        }
        vertex_t to_vtx = vm.at (to_id);
        to_vtx.add_neighbour_in (from_id);
        vm [to_id] = to_vtx;

        edge_id_t edge_id_str = Rcpp::as <edge_id_t> (edge_id [i]);

        double wt = weight [i];
        if (weight [i] < 0.0)
            wt = INFINITE_DOUBLE; // # nocov

        std::vector <double> weights;
        if (!has_times)
            weights = {dist [i], wt};
        else
            weights = {dist [i], wt, time [i], timew [i]};

        edge_t edge = edge_t (from_id, to_id, weights,
                edge_id_str, replacement_edges);

        edge_map.emplace (edge_id_str, edge);
        graph::add_to_v2e_map (vert2edge_map, from_id, edge_id_str);
        graph::add_to_v2e_map (vert2edge_map, to_id, edge_id_str);
    }

    return has_times;
}

//' identify_graph_components
//'
//' Identify initial graph components for each **vertex**
//' Identification for edges is subsequently perrformed with 
//' \code{rcpp_get_component_vector}.
//'
//' @param v unordered_map <vertex_id_t, vertex_t>
//' @param com component map from each vertex to component numbers
//' @noRd
size_t graph::identify_graph_components (vertex_map_t &v,
        std::unordered_map <vertex_id_t, size_t> &com)
{
    // initialize components map
    std::unordered_set <vertex_id_t> all_verts;
    for (auto it: v)
        all_verts.emplace (it.first);
    com.clear ();

    std::unordered_set <vertex_id_t> nbs_todo, nbs_done;
    nbs_todo.insert (*all_verts.begin ());
    size_t compnum = 0;
    while (all_verts.size () > 0)
    {
        Rcpp::checkUserInterrupt ();
        vertex_id_t vt = (*nbs_todo.begin ());
        all_verts.erase (vt);

        vertex_t vtx = v.find (vt)->second;
        std::unordered_set <vertex_id_t> nbs = vtx.get_all_neighbours ();
        for (auto nvtx: nbs)
        {
            com.emplace (nvtx, compnum);
            if (nbs_done.find (nvtx) == nbs_done.end ())
                nbs_todo.emplace (nvtx);
        }
        nbs_done.emplace (vt);
        com.emplace (vt, compnum);
        nbs_todo.erase (vt);

        if (nbs_todo.size () == 0 && all_verts.size () > 0)
        {
            nbs_todo.insert (*all_verts.begin ());
            compnum++;
        }
    }

    long int largest_id = 0;
    if (compnum > 0)
    {
        std::vector <size_t> comp_sizes (compnum + 1, 0);
        for (auto c: com)
            comp_sizes [c.second]++;
        auto maxi = std::max_element (comp_sizes.begin (), comp_sizes.end ());

        largest_id = std::distance (comp_sizes.begin (), maxi);
    }

    return static_cast <size_t> (largest_id);
}


//' rcpp_get_component_vector
//'
//' Get component numbers for each edge of graph
//'
//' @param graph graph to be processed; stripped down and standardised to five
//' columns
//'
//' @return Two vectors: one of edge IDs and one of corresponding component
//' numbers
//' @noRd
// [[Rcpp::export]]
Rcpp::List rcpp_get_component_vector (const Rcpp::DataFrame &graph)
{
    vertex_map_t vertices;
    edge_map_t edge_map;
    vert2edge_map_t vert2edge_map;

    bool has_times = graph::graph_from_df (graph, vertices, edge_map, vert2edge_map);
    has_times = false; // rm unused variable warning

    std::unordered_map <vertex_id_t, size_t> components;
    size_t largest_component =
        graph::identify_graph_components (vertices, components);
    largest_component++; // suppress unused variable warning

    // Then map component numbers of vertices onto edges
    std::unordered_map <edge_id_t, size_t> comp_nums;
    for (auto ve: vert2edge_map)
    {
        vertex_id_t vi = ve.first;
        std::unordered_set <edge_id_t> edges = ve.second;
        for (edge_id_t e: edges)
            comp_nums.emplace (e, components.find (vi)->second);
    }

    Rcpp::StringVector edge_id (comp_nums.size ());
    Rcpp::IntegerVector comp_num (comp_nums.size ());
    size_t i = 0;
    for (auto cn: comp_nums)
    {
        edge_id (i) = cn.first;
        comp_num (i) = static_cast <int> (cn.second) + 1L; // 1-indexed for R
        i++;
    }

    return Rcpp::List::create (
            Rcpp::Named ("edge_id") = edge_id,
            Rcpp::Named ("edge_component") = comp_num);
}

//' rcpp_unique_rownames
//'
//' Construct vertex (from, to) ID values from unique pairs of coordinates
//' rounded to <precision>. Used when vertices have no ID values.
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::DataFrame rcpp_unique_rownames (Rcpp::DataFrame xyfrom,
        Rcpp::DataFrame xyto,
        const int precision = 10)
{
    const size_t prec_size_t = static_cast <size_t> (precision);

    std::vector <double> xf = xyfrom ["x"],
    yf = xyfrom ["y"],
    xt = xyto ["x"],
    yt = xyto ["y"];
    const size_t n = static_cast <size_t> (xyfrom.nrow ());

    std::vector <std::string> s_f (n), s_t (n);
    std::vector <std::string> i_f (n), i_t (n); // indices as rownames
    std::unordered_map <std::string, std::string> xynames;
    size_t count = 0;
    for (size_t i = 0; i < n; i++)
    {
        std::string xfs = std::to_string (xf [i]),
            yfs = std::to_string (yf [i]),
            xts = std::to_string (xt [i]),
            yts = std::to_string (yt [i]);
        s_f [i] = xfs.substr(0, xfs.find(".") + prec_size_t + 1) +
            yfs.substr (0, yfs.find(".") + prec_size_t + 1);
        s_t [i] = xts.substr (0, xts.find(".") + prec_size_t + 1) +
            yts.substr (0, yts.find(".") + prec_size_t + 1);

        if (xynames.find (s_f [i]) != xynames.end ())
        {
            i_f [i] = xynames.at (s_f [i]);
        } else
        {
            i_f [i] = std::to_string (count);
            xynames.emplace (s_f [i], i_f [i]);
            count++;
        }

        if (xynames.find (s_t [i]) != xynames.end ())
        {
            i_t [i] = xynames.at (s_t [i]);
        } else
        {
            i_t [i] = std::to_string (count);
            xynames.emplace (s_t [i], i_t [i]);
            count++;
        }
    }

    Rcpp::DataFrame res = Rcpp::DataFrame::create (
            Rcpp::Named ("from_id") = i_f,
            Rcpp::Named ("to_id") = i_t,
            Rcpp::_["stringsAsFactors"] = false);
    return res;
}
