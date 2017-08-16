/***************************************************************************
 *  Project:    osmdata
 *  File:       lines-as-network.cpp
 *  Language:   C++
 *
 *  osmdata is free software: you can redistribute it and/or modify it under
 *  the terms of the GNU General Public License as published by the Free
 *  Software Foundation, either version 3 of the License, or (at your option)
 *  any later version.
 *
 *  osmdata is distributed in the hope that it will be useful, but WITHOUT ANY
 *  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 *  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 *  details.
 *
 *  You should have received a copy of the GNU General Public License along with
 *  osm-router.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Author:     Mark Padgham 
 *  E-Mail:     mark.padgham@email.com 
 *
 *  Description:    Convert sf linestring collection to data.frame of network
 *                  connections
 *
 *  Limitations:
 *
 *  Dependencies:       none (rapidXML header included in osmdatar)
 *
 *  Compiler Options:   -std=c++11
 ***************************************************************************/

#include <string>
#include <cmath>

#include <Rcpp.h>

const float INFINITE_FLOAT =  std::numeric_limits<float>::max ();
const int INFINITE_INT =  std::numeric_limits<int>::max ();

// Haversine great circle distance between two points
float haversine (float x1, float y1, float x2, float y2)
{
    float xd = (x2 - x1) * M_PI / 180.0;
    float yd = (y2 - y1) * M_PI / 180.0;
    float d = sin (yd / 2.0) * sin (yd / 2.0) + cos (y2 * M_PI / 180.0) *
        cos (y1 * M_PI / 180.0) * sin (xd / 2.0) * sin (xd / 2.0);
    d = 2.0 * 3671.0 * asin (sqrt (d));
    return (d);
}

//' rcpp_sf_as_network
//'
//' Return OSM data in Simple Features format
//'
//' @param sf_lines An sf collection of LINESTRING objects
//' @param pr Rcpp::DataFrame containing the weighting profile
//'
//' @return Rcpp::List objects of OSM data
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::List rcpp_sf_as_network (const Rcpp::List &sf_lines,
        Rcpp::DataFrame pr)
{
    std::map <std::string, float> profile;
    Rcpp::StringVector hw = pr [1];
    Rcpp::NumericVector val = pr [2];
    for (int i = 0; i != hw.size (); i ++)
        profile.insert (std::make_pair (std::string (hw [i]), val [i]));

    Rcpp::CharacterVector nms = sf_lines.attr ("names");
    if (nms [nms.size () - 1] != "geometry")
        throw std::runtime_error ("sf_lines have no geometry component");
    if (nms [0] != "osm_id")
        throw std::runtime_error ("sf_lines have no osm_id component");
    int one_way_index = -1;
    int one_way_bicycle_index = -1;
    int highway_index = -1;
    for (int i = 0; i < nms.size (); i++)
    {
        if (nms [i] == "oneway")
            one_way_index = i;
        if (nms [i] == "oneway.bicycle")
            one_way_bicycle_index = i;
        if (nms [i] == "highway")
            highway_index = i;
    }
    Rcpp::CharacterVector ow = NULL;
    Rcpp::CharacterVector owb = NULL;
    Rcpp::CharacterVector highway = NULL;
    if (one_way_index >= 0)
        ow = sf_lines [one_way_index];
    if (one_way_bicycle_index >= 0)
        owb = sf_lines [one_way_bicycle_index];
    if (highway_index >= 0)
        highway = sf_lines [highway_index];
    if (ow.size () > 0)
    {
        if (ow.size () == owb.size ())
        {
            for (unsigned i = 0; i != ow.size (); ++ i)
                if (ow [i] == "NA" && owb [i] != "NA")
                    ow [i] = owb [i];
        } else if (owb.size () > ow.size ())
            ow = owb;
    }

    Rcpp::List geoms = sf_lines [nms.size () - 1];
    std::vector<bool> isOneWay (geoms.length ());
    std::fill (isOneWay.begin (), isOneWay.end (), false);
    // Get dimension of matrix
    size_t nrows = 0;
    int ngeoms = 0;
    for (auto g = geoms.begin (); g != geoms.end (); ++g)
    {
        // Rcpp uses an internal proxy iterator here, NOT a direct copy
        Rcpp::NumericMatrix gi = (*g);
        int rows = gi.nrow () - 1;
        nrows += rows;
        if (ngeoms < ow.size ())
        {
            if (!(ow [ngeoms] == "yes" || ow [ngeoms] == "-1"))
            {
                nrows += rows;
                isOneWay [ngeoms] = true;
            }
        }
        ngeoms ++;
    }

    Rcpp::NumericMatrix nmat = Rcpp::NumericMatrix (Rcpp::Dimension (nrows, 6));
    Rcpp::CharacterMatrix idmat = Rcpp::CharacterMatrix (Rcpp::Dimension (nrows,
                3));

    nrows = 0;
    ngeoms = 0;
    int fake_id = 0;
    for (auto g = geoms.begin (); g != geoms.end (); ++ g)
    {
        Rcpp::NumericMatrix gi = (*g);
        std::string hway = std::string (highway [ngeoms]);
        float hw_factor = profile [hway];
        if (hw_factor == 0.0) hw_factor = 1e-5;
        hw_factor = 1.0 / hw_factor;

        Rcpp::List ginames = gi.attr ("dimnames");
        Rcpp::CharacterVector rnms;
        if (ginames.length () > 0)
            rnms = ginames [0];
        else
        {
            rnms = Rcpp::CharacterVector (gi.nrow ());
            for (int i = 0; i < gi.nrow (); i ++)
                rnms [i] = fake_id ++;
        }
        if (rnms.size () != gi.nrow ())
            throw std::runtime_error ("geom size differs from rownames");

        for (int i = 1; i < gi.nrow (); i ++)
        {
            float d = haversine (gi (i-1, 0), gi (i-1, 1), gi (i, 0),
                    gi (i, 1));
            nmat (nrows, 0) = gi (i-1, 0);
            nmat (nrows, 1) = gi (i-1, 1);
            nmat (nrows, 2) = gi (i, 0);
            nmat (nrows, 3) = gi (i, 1);
            nmat (nrows, 4) = d;
            nmat (nrows, 5) = d * hw_factor;
            idmat (nrows, 0) = rnms (i-1);
            idmat (nrows, 1) = rnms (i);
            idmat (nrows, 2) = hway;
            nrows ++;
            if (isOneWay [ngeoms])
            {
                nmat (nrows, 0) = gi (i, 0);
                nmat (nrows, 1) = gi (i, 1);
                nmat (nrows, 2) = gi (i-1, 0);
                nmat (nrows, 3) = gi (i-1, 1);
                nmat (nrows, 4) = d;
                nmat (nrows, 5) = d * hw_factor;
                idmat (nrows, 0) = rnms (i);
                idmat (nrows, 1) = rnms (i-1);
                idmat (nrows, 2) = hway;
                nrows ++;
            }
        }
        ngeoms ++;
    }

    Rcpp::List res (2);
    res [0] = nmat;
    res [1] = idmat;

    return res;
}

//' rcpp_points_index
//'
//' Get index of nearest vertices to list of points
//'
//' @param graph Rcpp::DataFrame containing the graph
//' @param pts Rcpp::DataFrame containing the routing points
//'
//' @return Rcpp::NumericVector index into graph of nearest points
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_points_index (const Rcpp::DataFrame &xy,
        Rcpp::DataFrame &pts)
{
    Rcpp::NumericVector ptx = pts ["x"];
    Rcpp::NumericVector pty = pts ["y"];

    Rcpp::NumericVector vtx = xy ["x"];
    Rcpp::NumericVector vty = xy ["y"];

    Rcpp::NumericVector index (pts.nrow ());

    for (int i = 0; i < pts.nrow (); i++) // Rcpp::nrow is int!
    {
        float dmin = INFINITE_FLOAT;
        int jmin = INFINITE_INT;
        for (int j = 0; j < xy.nrow (); j++)
        {
            float dij = (vtx [j] - ptx [i]) * (vtx [j] - ptx [i]) +
                (vty [j] - pty [i]) * (vty [j] - pty [i]);
            if (dij < dmin)
            {
                dmin = dij;
                jmin = j;
            }
        }
        if (jmin == INFINITE_INT)
            Rcpp::Rcout << "---ERROR---" << std::endl;
        index (i) = jmin;
    }

    return index;
}
