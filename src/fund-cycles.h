#pragma once

#include <iostream>
#include <vector>
#include <list>
#include <stack>
#include <memory> // unique_ptr

#include <Rcpp.h>

namespace graph
{

class AdjacencyMatrix
{
    private:

        const size_t m_nNodes;
        size_t m_nEdges;
        std::vector <bool> m_Adjacencies;
        const long long m_nRows; // Equivalent #rows of adjacency matrix

    public:
        // ---------- Constructor and copy-constructor:
        inline AdjacencyMatrix (const size_t nNodes) :
            m_nNodes (nNodes),
            m_nEdges (0),
            m_Adjacencies ((nNodes * (nNodes - 1)) / 2),
            m_nRows (static_cast<const long long>(1 + 2 * (nNodes - 2)))
            { }

        inline AdjacencyMatrix (const AdjacencyMatrix& m) :
            m_nNodes (m.m_nNodes),
            m_nEdges (m.m_nEdges),
            m_Adjacencies (m.m_Adjacencies),
            m_nRows (m.m_nRows)
            { }

        // ---------- Primary functions
        inline void connect (size_t i, size_t j);
        inline void disconnect (size_t i, size_t j);
        inline bool is_connected (size_t i, size_t j) const;
        inline bool operator() (size_t i, size_t j) const
        {
            return is_connected (i, j);
        }
        inline size_t getNumEdges() const { return m_nEdges; }

        // ---------- Primary operators
        // ---------- (1) XOR
        inline AdjacencyMatrix operator ^ (const AdjacencyMatrix& adj_comp) const
        {
            AdjacencyMatrix result (m_nNodes);
            for (size_t i = 0; i < m_Adjacencies.size(); ++i)
            {
                if ((m_Adjacencies [i] || adj_comp.m_Adjacencies [i]) &&
                        (m_Adjacencies [i] != adj_comp.m_Adjacencies [i]))
                {
                    result.m_Adjacencies [i] = true;
                    result.m_nEdges++;
                }
            }
            return result;
        }

        // ---------- (1) XOR each component
        inline AdjacencyMatrix& operator ^= (const AdjacencyMatrix& adj_comp)
        {
            m_nEdges = 0;
            for (size_t i = 0; i < m_Adjacencies.size(); ++i)
            {
                if ((m_Adjacencies [i] || adj_comp.m_Adjacencies [i]) &&
                        (m_Adjacencies [i] != adj_comp.m_Adjacencies [i]))
                {
                    m_Adjacencies[i] = true;
                    ++m_nEdges;
                }
                else
                    m_Adjacencies [i] = false;
            }
            return *this;
        }

        inline AdjacencyMatrix& operator = (const AdjacencyMatrix& adj_comp)
        {
            m_Adjacencies = adj_comp.m_Adjacencies;
            m_nEdges = adj_comp.m_nEdges;
            return *this;
        }

    private:
        inline size_t get_adj_index (const size_t i, const size_t j) const
        {
            if (i >= m_nNodes || j >= m_nNodes || i == j)
                throw std::out_of_range("get_adj_index: (i, j) must be < nNodes AND they must not be equal!"); // # nocov

            // (j - i * (i - 1 - 2 * (N - 2)) / 2) - 1
            long long li = static_cast<long long>(i),
                 lj = static_cast<long long>(j);
            if (i < j)
                return static_cast<size_t>((lj - li * (li - m_nRows) / 2) - 1);
            else // swap indices
                return static_cast<size_t>((li - lj * (lj - m_nRows) / 2) - 1);
        }

};

inline void AdjacencyMatrix::connect (size_t i, size_t j)
{
    if (AdjacencyMatrix::m_Adjacencies [AdjacencyMatrix::get_adj_index (i, j)])
        return;

    AdjacencyMatrix::m_Adjacencies [AdjacencyMatrix::get_adj_index (i, j)] = true;
    AdjacencyMatrix::m_nEdges++;
}

inline void AdjacencyMatrix::disconnect (size_t i, size_t j)
{
    if (!AdjacencyMatrix::m_Adjacencies [AdjacencyMatrix::get_adj_index (i, j)])
        return;

    AdjacencyMatrix::m_Adjacencies [AdjacencyMatrix::get_adj_index (i, j)] = false;
    AdjacencyMatrix::m_nEdges--;
}

inline bool AdjacencyMatrix::is_connected (size_t i, size_t j) const
{
    if (i == j)
        return false;
    return AdjacencyMatrix::m_Adjacencies [AdjacencyMatrix::get_adj_index (i, j)];
}

template <class TObject>
class Graph
{
    public:

        typedef std::vector <TObject> NodeArray;
        typedef std::vector <AdjacencyMatrix> CycleArray;
        typedef std::list <const TObject*> NodePath;

    private:

        NodeArray 				m_nodeArray;
        AdjacencyMatrix 	    m_adjMat;
        struct TreeNode
        {
            size_t index;
            TreeNode* parent;
        };

    public:

        CycleArray 				m_fundamentalCycles;

        inline Graph (const std::vector <TObject> &nodeArray,
                const size_t nNodes,
                std::vector <size_t> &edgeArray,
                const size_t nEdges) :
            m_adjMat (nNodes)
        {
            m_nodeArray.reserve (nNodes);
            for (size_t i = 0; i < nNodes; ++i)
                m_nodeArray.push_back (nodeArray[i]);
            for (size_t i = 0; i < nEdges; ++i)
                m_adjMat.connect (edgeArray [2 * i], edgeArray [2 * i + 1]);
        }

        void computeFundamentalCycles();
        NodePath cycleMatrix2nodePath (const AdjacencyMatrix& m) const;
        inline const CycleArray& getFundamentalCycles() const
        {
            return m_fundamentalCycles;
        }
        inline size_t getNumNodes() const { return m_nodeArray.size(); }

    private:

        void cycleMatrix2nodePath_recursion(const AdjacencyMatrix& m,
                NodePath& path,
                const size_t i,
                const size_t prev,
                const size_t start_node) const;

        inline void path_to_tree_root(TreeNode* pNode, AdjacencyMatrix& adjMat);
};

template <class TObject>
void Graph<TObject>::computeFundamentalCycles()
{
    if (!Graph<TObject>::m_fundamentalCycles.empty())
        return;

    std::unique_ptr <TreeNode []>
        aTree (new TreeNode [Graph::m_nodeArray.size ()]);
    std::stack<size_t> nodeStack;

    nodeStack.push(0); // Start with first node

    AdjacencyMatrix adjMat = Graph<TObject>::m_adjMat; // copy of matrix

    // Initialise tree:
    for (size_t i = 0; i < m_nodeArray.size(); ++i)
    {
        aTree [i].parent = &aTree [i];
        aTree [i].index = i;
    }

    while (nodeStack.size() > 0)
    {
        size_t currentNodeIndex = nodeStack.top ();
        nodeStack.pop ();
        TreeNode& currentTreeNode = aTree [currentNodeIndex];

        // Iterate though edges connecting this node:
        for (size_t j = 0; j < Graph<TObject>::m_nodeArray.size(); ++j)
        {
            if (!adjMat.is_connected (currentNodeIndex, j))
                continue;

            if (aTree [j].parent != &aTree [j])
            {
                // Fundamental Cycle found: Get unique paths from
                // both nodes within the spanning tree!
                AdjacencyMatrix pi (m_nodeArray.size ()),
                                pj (m_nodeArray.size ());
                path_to_tree_root (&aTree [currentNodeIndex], pi);
                path_to_tree_root (&aTree [j], pj);

                // Connection between currentNodeIndex and j has to
                // be inserted to ONE of the two paths (which one
                // does not matter)
                pi.connect (currentNodeIndex, j);

                // XOR the 2 matrices to get fundamental cycle
                Graph<TObject>::m_fundamentalCycles.push_back (pi ^ pj);
            }
            else // node not in tree, so add it
            {
                aTree [j].parent = &currentTreeNode;
                nodeStack.push (j);
            }
            adjMat.disconnect (currentNodeIndex, j);
        }
    }
}

/* Transform AdjacencyMatrix containing indices of the graph to a path
 * containing pointers to the actual objects. Each node in AdjacencyMatrix must
 * contain either 0 or 2 edges, and path will always return to starting point.
 */
template <class TObject>
typename Graph<TObject>::NodePath Graph<TObject>::cycleMatrix2nodePath
                            (const AdjacencyMatrix& m) const
{
    Graph<TObject>::NodePath path;
    for (size_t i = 0; i < Graph<TObject>::m_nodeArray.size(); ++i)
    {
        for (size_t j = 0; j < Graph<TObject>::m_nodeArray.size(); ++j)
        {
            if (m.is_connected (i, j))
            {
                path.push_back (&Graph<TObject>::m_nodeArray [i]);
                path.push_back (&Graph<TObject>::m_nodeArray [j]);
                Graph<TObject>::cycleMatrix2nodePath_recursion (m, path, j, i, i);

                return path;
            }
        }
    }
    // Only !return is if matrix does not contain any edges
    throw std::runtime_error("Graph::cycleMatrix2nodePath(): Given Cycle Matrix does not contain any edges!"); // # nocov
}

template <class TObject>
void Graph<TObject>::cycleMatrix2nodePath_recursion(const AdjacencyMatrix& m,
        Graph<TObject>::NodePath& path,
        const size_t i, // desired node
        const size_t prev,
        const size_t start_node) const // start_node defines stop criterion
{
    for (size_t j = 0; j < Graph<TObject>::m_nodeArray.size(); ++j)
    {
        if (m.is_connected (i, j) && j != prev)
        {
            path.push_back (&m_nodeArray [j]);
            if (j != start_node)
                Graph<TObject>::cycleMatrix2nodePath_recursion (m, path, j, i,
                        start_node);
            return;
        }
    }
    // Only not return if in a dead end:
    throw std::runtime_error("Graph::cycleMatrix2nodePath_recursion(): Found a dead end!"); // # nocov
}

// Function recursively finds the unique path within the tree from the given node to the root of the tree
template <class TObject>
inline void Graph<TObject>::path_to_tree_root(Graph<TObject>::TreeNode* pNode, AdjacencyMatrix& adjMat)
{
    if (pNode->parent != pNode)
    {
        adjMat.connect(pNode->index, pNode->parent->index);
        Graph<TObject>::path_to_tree_root(pNode->parent, adjMat);
    }
}


} // end namespace

Rcpp::List rcpp_fundamental_cycles (Rcpp::DataFrame graph,
        Rcpp::DataFrame verts);
