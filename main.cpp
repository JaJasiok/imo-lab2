#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include <climits>
#include <algorithm>
#include <unordered_set>
#include <random>
#include <numeric>
#include <variant>
#include <chrono>

using namespace std;

vector<vector<int>> readKroaFile(const string &filename)
{
    ifstream file(filename);
    string line;
    vector<vector<int>> verticesCoords;

    if (file.is_open())
    {
        // Skip the header lines
        for (int i = 0; i < 6; i++)
        {
            getline(file, line);
        }

        // Read the coordinates and populate the cost matrix
        while (getline(file, line) && line != "EOF")
        {
            istringstream iss(line);
            int index, x, y;
            iss >> index >> x >> y;
            verticesCoords.push_back({x, y});
        }

        verticesCoords.pop_back();
        
        file.close();
    }
    else
    {
        cout << "Failed to open file: " << filename << endl;
    }

    return verticesCoords;
}

vector<vector<int>> createDistanceMatrix(const vector<vector<int>> &verticesCoords)
{
    int numVertices = verticesCoords.size();
    vector<vector<int>> distanceMatrix(numVertices, vector<int>(numVertices));

    for (int i = 0; i < numVertices; i++)
    {
        for (int j = 0; j < numVertices; j++)
        {
            double x1 = verticesCoords[i][0];
            double y1 = verticesCoords[i][1];
            double x2 = verticesCoords[j][0];
            double y2 = verticesCoords[j][1];

            double distance = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
            distanceMatrix[i][j] = int(round(distance));
        }
    }

    return distanceMatrix;
}

class Vertex
{
public:
    int id;
    vector<Vertex *> neighbours;

    Vertex(int _id) : id(_id) {}
};

class Edge
{
public:
    Vertex *src;
    Vertex *dest;
    int distance;

    Edge(Vertex *_src, Vertex *_dest, int _distance) : src(_src), dest(_dest), distance(_distance) {}
};

class Graph
{
public:
    vector<Vertex *> vertices;
    vector<Edge *> edges;
    int distance = 0;

    void addVertex(Vertex *v)
    {
        if (v)
        {
            if (find(vertices.begin(), vertices.end(), v) == vertices.end())
            {
                vertices.push_back(v);
            }
        }
    }

    void addEdge(Vertex *src, Vertex *dest, int distance = 0)
    {
        if (src && dest && src != dest)
        {
            for (Edge *e : edges)
            {
                if ((e->src == src && e->dest == dest) || (e->src == dest && e->dest == src))
                {
                    return;
                }
            }

            Edge *e = new Edge(src, dest, distance);
            src->neighbours.push_back(dest);
            dest->neighbours.push_back(src);
            this->distance += distance;
            edges.push_back(e);
        }
    }

    void removeVertex(Vertex *v)
    {
        if (v)
        {
            vertices.erase(remove(vertices.begin(), vertices.end(), v), vertices.end());
            for (Vertex *neighbour : v->neighbours)
            {
                neighbour->neighbours.erase(remove(neighbour->neighbours.begin(), neighbour->neighbours.end(), v), neighbour->neighbours.end());
            }
            for (Edge *e : edges)
            {
                if (e->src == v || e->dest == v)
                {
                    removeEdge(e->src, e->dest);
                }
            }
        }
    }

    void normalizeEdges()
    {
        for (Edge *e1 : edges)
        {
            for (Edge *e2 : edges)
            {
                if (e1 != e2)
                {
                    if (e1->src == e2->src || e1->dest == e2->dest)
                    {
                        swap(e2->dest, e2->src);
                        break;
                    }
                }
            }
        }
    }

    void normalizeVertices()
    {
        for (Vertex *v : vertices)
        {
            for (Edge *e : edges)
            {
                if (e->src == v)
                {
                    if (v->neighbours[1] != e->dest)
                    {
                        swap(v->neighbours[0], v->neighbours[1]);
                    }
                }
            }
        }
    }

    void normalizeGraph()
    {
        normalizeEdges();
        normalizeEdges();
        normalizeVertices();
    }

    // void removeEdge(Vertex *src, Vertex *dest)
    // {
    //     if (src && dest)
    //     {
    //         src->neighbours.erase(remove(src->neighbours.begin(), src->neighbours.end(), dest), src->neighbours.end());
    //         dest->neighbours.erase(remove(dest->neighbours.begin(), dest->neighbours.end(), src), dest->neighbours.end());
    //         for (Edge *e : edges)
    //         {
    //             if ((e->src == src && e->dest == dest) || (e->src == dest && e->dest == src))
    //             {
    //                 edges.erase(remove(edges.begin(), edges.end(), e), edges.end());
    //                 this->distance -= e->distance;
    //                 break;
    //             }
    //         }
    //     }
    // }

    void removeEdge(Vertex *src, Vertex *dest)
    {
        if (src && dest)
        {
            src->neighbours.erase(remove(src->neighbours.begin(), src->neighbours.end(), dest), src->neighbours.end());
            dest->neighbours.erase(remove(dest->neighbours.begin(), dest->neighbours.end(), src), dest->neighbours.end());
            for (Edge *e : edges)
            {
                if (e->src == src && e->dest == dest || e->src == dest && e->dest == src)
                {
                    edges.erase(remove(edges.begin(), edges.end(), e), edges.end());
                    this->distance -= e->distance;
                }
            }
        }
    }

    Vertex *findVertex(int id)
    {
        for (Vertex *v : vertices)
        {
            if (v->id == id)
                return v;
        }
        return nullptr;
    }

    // bool hasCycle()
    // {
    //     unordered_set<Vertex *> visited;
    //     for (Vertex *v : vertices)
    //     {
    //         if (!visited.count(v) && hasCycleUtil(v, nullptr, visited))
    //             return true;
    //     }
    //     return false;
    // }

    // bool hasCycleUtil(Vertex *v, Vertex *parent, unordered_set<Vertex *> &visited)
    // {
    //     visited.insert(v);
    //     for (Vertex *neighbour : v->neighbours)
    //     {
    //         if (!visited.count(neighbour))
    //         {
    //             if (hasCycleUtil(neighbour, v, visited))
    //                 return true;
    //         }
    //         else if (neighbour != parent)
    //         {
    //             return true;
    //         }
    //     }
    //     return false;
    // }
};

void saveGraphs(const vector<Graph> &graphs, const string &filename)
{
    ofstream file(filename);
    if (file.is_open())
    {
        for (Graph g : graphs)
        {
            for (Edge *e : g.edges)
            {
                file << e->src->id << " " << e->dest->id << endl;
            }
            file << endl;
        }
        file.close();
    }
    else
    {
        cout << "Failed to open file: " << filename << endl;
    }
}

vector<Graph> randomCycles(const vector<vector<int>> &distanceMatrix)
{
    std::vector<int> numbers(100);
    iota(numbers.begin(), numbers.end(), 0);

    random_shuffle(numbers.begin(), numbers.end());

    vector<Graph> cycles(2);

    std::vector<int> group1(numbers.begin(), numbers.begin() + 50);
    std::vector<int> group2(numbers.begin() + 50, numbers.end());

    for (int i : group1)
    {
        Vertex *v = new Vertex(i);
        cycles[0].addVertex(v);
    }

    for (int i : group2)
    {
        Vertex *v = new Vertex(i);
        cycles[1].addVertex(v);
    }

    for (int i = 0; i < 50; i++)
    {
        if (i < 49)
        {
            Vertex *v1 = cycles[0].vertices[i];
            Vertex *v2 = cycles[0].vertices[i + 1];
            cycles[0].addEdge(v1, v2, distanceMatrix[v1->id][v2->id]);
        }
        else
        {
            Vertex *v1 = cycles[0].vertices[i];
            Vertex *v2 = cycles[0].vertices[0];
            cycles[0].addEdge(v1, v2, distanceMatrix[v1->id][v2->id]);
        }
    }

    for (int i = 0; i < 50; i++)
    {
        if (i < 49)
        {
            Vertex *v1 = cycles[1].vertices[i];
            Vertex *v2 = cycles[1].vertices[i + 1];
            cycles[1].addEdge(v1, v2, distanceMatrix[v1->id][v2->id]);
        }
        else
        {
            Vertex *v1 = cycles[1].vertices[i];
            Vertex *v2 = cycles[1].vertices[0];
            cycles[1].addEdge(v1, v2, distanceMatrix[v1->id][v2->id]);
        }
    }

    return cycles;
}

vector<Graph> greedyCycles(const vector<vector<int>> &distanceMatrix, int startId)
{
    int numVertices = distanceMatrix.size();
    vector<Graph> cycles(2);
    vector<bool> visited(numVertices, false);

    // Choose the first vertex randomly
    Vertex *startVertex1 = new Vertex(startId);
    cycles[0].addVertex(startVertex1);
    visited[startVertex1->id] = true;

    // Find the furthest vertex from startVertex1
    int furthestVertex = -1;
    int maxDistance = -1;
    for (int j = 0; j < numVertices; j++)
    {
        if (!visited[j] && distanceMatrix[startVertex1->id][j] > maxDistance)
        {
            maxDistance = distanceMatrix[startVertex1->id][j];
            furthestVertex = j;
        }
    }

    Vertex *startVertex2 = new Vertex(furthestVertex);
    cycles[1].addVertex(startVertex2);
    visited[startVertex2->id] = true;

    // Find the nearest neighbour for each starting vertex
    for (Graph &cycle : cycles)
    {
        Vertex *vertex = cycle.vertices.front();
        int vertexId = vertex->id;
        int minDistance = INT_MAX;
        int nearestNeighbourId = -1;

        for (int j = 0; j < numVertices; j++)
        {
            if (!visited[j])
            {
                if (distanceMatrix[vertexId][j] < minDistance)
                {
                    minDistance = distanceMatrix[vertexId][j];
                    nearestNeighbourId = j;
                }
            }
        }

        Vertex *nearestNeighbour = new Vertex(nearestNeighbourId);

        cycle.addVertex(nearestNeighbour);
        cycle.addEdge(vertex, nearestNeighbour, minDistance);
        visited[nearestNeighbourId] = true;
    }

    // Building the rest of the cycle
    for (int i = 0; i < numVertices - 4; i++)
    {
        int minDistance = INT_MAX;
        int vertexId = -1;
        pair<Edge *, Graph *> minPair;

        vector<pair<Edge *, Graph *>> edgesInGraphs;
        if (cycles[0].vertices.size() == numVertices / 2)
        {
            for (Edge *e : cycles[1].edges)
            {
                edgesInGraphs.push_back({e, &cycles[1]});
            }
        }
        else if (cycles[1].vertices.size() == numVertices / 2)
        {
            for (Edge *e : cycles[0].edges)
            {
                edgesInGraphs.push_back({e, &cycles[0]});
            }
        }
        else
        {
            for (Graph &cycle : cycles)
            {
                for (Edge *e : cycle.edges)
                {
                    edgesInGraphs.push_back({e, &cycle});
                }
            }
        }

        for (int j = 0; j < numVertices; j++)
        {
            if (!visited[j])
            {
                for (pair<Edge *, Graph *> singlePair : edgesInGraphs)
                {
                    Edge *e = singlePair.first;
                    int distanceSum = cycles[0].distance + cycles[1].distance;
                    if (distanceSum - e->distance + distanceMatrix[e->dest->id][j] + distanceMatrix[e->src->id][j] < minDistance)
                    {
                        minDistance = distanceSum - e->distance + distanceMatrix[e->dest->id][j] + distanceMatrix[e->src->id][j];
                        vertexId = j;
                        minPair = singlePair;
                    }
                }
            }
        }

        Graph *cyclePtr = minPair.second;
        Graph &cycle = *cyclePtr;

        Edge *minEdge = minPair.first;

        Vertex *newVertex = new Vertex(vertexId);
        cycle.addVertex(newVertex);
        cycle.addEdge(minEdge->src, newVertex, distanceMatrix[minEdge->src->id][vertexId]);
        cycle.addEdge(newVertex, minEdge->dest, distanceMatrix[newVertex->id][minEdge->dest->id]);

        if (cycle.vertices.size() > 3)
        {
            cycle.removeEdge(minEdge->src, minEdge->dest);
        }

        visited[vertexId] = true;
    }

    for (Graph &cycle : cycles)
    {
        cycle.normalizeGraph();
    }
    return cycles;
}

void swapVerticesBetweenCycles(int i, int j, vector<Graph> &graphs, const vector<vector<int>> &distanceMatrix)
{
    Graph *graph1 = &graphs[0];
    Graph *graph2 = &graphs[1];

    Vertex *vertex1 = graph1->findVertex(i);
    Vertex *vertex2 = graph2->findVertex(j);

    vector<Vertex *> neighbours1 = vertex1->neighbours;
    vector<Vertex *> neighbours2 = vertex2->neighbours;

    Vertex *vertex1Neighbour1 = neighbours1[0];
    Vertex *vertex1Neighbour2 = neighbours1[1];

    Vertex *vertex2Neighbour1 = neighbours2[0];
    Vertex *vertex2Neighbour2 = neighbours2[1];

    graph1->removeVertex(vertex1);
    graph1->removeEdge(vertex1Neighbour1, vertex1);
    graph1->removeEdge(vertex1, vertex1Neighbour2);

    graph2->removeVertex(vertex2);
    graph2->removeEdge(vertex2Neighbour1, vertex2);
    graph2->removeEdge(vertex2, vertex2Neighbour2);

    graph1->addVertex(vertex2);
    graph1->addEdge(vertex1Neighbour1, vertex2, distanceMatrix[vertex1Neighbour1->id][vertex2->id]);
    graph1->addEdge(vertex2, vertex1Neighbour2, distanceMatrix[vertex2->id][vertex1Neighbour2->id]);

    graph2->addVertex(vertex1);
    graph2->addEdge(vertex2Neighbour1, vertex1, distanceMatrix[vertex2Neighbour1->id][vertex1->id]);
    graph2->addEdge(vertex1, vertex2Neighbour2, distanceMatrix[vertex1->id][vertex2Neighbour2->id]);
}

void swapVerticesInCycle(int i, int j, Graph &graph, const vector<vector<int>> &distanceMatrix)
{
    Vertex *vertex1 = graph.findVertex(i);
    Vertex *vertex2 = graph.findVertex(j);

    vector<Vertex *> neighbours1 = vertex1->neighbours;
    vector<Vertex *> neighbours2 = vertex2->neighbours;

    Vertex *vertex1Neighbour1 = neighbours1[0];
    Vertex *vertex1Neighbour2 = neighbours1[1];

    Vertex *vertex2Neighbour1 = neighbours2[0];
    Vertex *vertex2Neighbour2 = neighbours2[1];

    if (vertex1Neighbour1 == vertex2 || vertex1Neighbour2 == vertex2)
    {
        if (vertex1Neighbour1 == vertex2 && vertex2Neighbour1 == vertex1)
        {
            graph.removeEdge(vertex1Neighbour2, vertex1);
            graph.removeEdge(vertex2, vertex2Neighbour2);

            graph.addEdge(vertex1, vertex2Neighbour2, distanceMatrix[vertex1->id][vertex2Neighbour2->id]);
            graph.addEdge(vertex1Neighbour2, vertex2, distanceMatrix[vertex1Neighbour2->id][vertex2->id]);
        }
        else if (vertex1Neighbour2 == vertex2 && vertex2Neighbour1 == vertex1)
        {
            graph.removeEdge(vertex1Neighbour1, vertex1);
            graph.removeEdge(vertex2, vertex2Neighbour2);

            graph.addEdge(vertex1, vertex2Neighbour2, distanceMatrix[vertex1->id][vertex2Neighbour2->id]);
            graph.addEdge(vertex1Neighbour1, vertex2, distanceMatrix[vertex1Neighbour1->id][vertex2->id]);
        }
        else if (vertex1Neighbour1 == vertex2 && vertex2Neighbour2 == vertex1)
        {
            graph.removeEdge(vertex1Neighbour2, vertex1);
            graph.removeEdge(vertex2Neighbour1, vertex2);

            graph.addEdge(vertex1, vertex2Neighbour1, distanceMatrix[vertex1->id][vertex2Neighbour1->id]);
            graph.addEdge(vertex1Neighbour2, vertex2, distanceMatrix[vertex1Neighbour2->id][vertex2->id]);
        }
        else if (vertex1Neighbour2 == vertex2 && vertex2Neighbour2 == vertex1)
        {
            graph.removeEdge(vertex1Neighbour1, vertex1);
            graph.removeEdge(vertex2Neighbour1, vertex2);

            graph.addEdge(vertex1, vertex2Neighbour1, distanceMatrix[vertex1->id][vertex2Neighbour1->id]);
            graph.addEdge(vertex1Neighbour1, vertex2, distanceMatrix[vertex1Neighbour1->id][vertex2->id]);
        }
    }
    else
    {
        graph.removeEdge(vertex1Neighbour1, vertex1);
        graph.removeEdge(vertex1, vertex1Neighbour2);
        graph.removeEdge(vertex2Neighbour1, vertex2);
        graph.removeEdge(vertex2, vertex2Neighbour2);

        graph.addEdge(vertex1Neighbour1, vertex2, distanceMatrix[vertex1Neighbour1->id][vertex2->id]);
        graph.addEdge(vertex2, vertex1Neighbour2, distanceMatrix[vertex2->id][vertex1Neighbour2->id]);
        graph.addEdge(vertex2Neighbour1, vertex1, distanceMatrix[vertex2Neighbour1->id][vertex1->id]);
        graph.addEdge(vertex1, vertex2Neighbour2, distanceMatrix[vertex1->id][vertex2Neighbour2->id]);
    }
}

vector<Graph> steepyLocalSearchNeighbourhood1(vector<Graph> &cycles, const vector<vector<int>> &distanceMatrix)
{
    int bestDelta = 0;
    
    do
    {
        vector<pair<pair<Vertex *, Vertex *>, Graph *>> moves;

        for (Vertex *vertex1 : cycles[0].vertices)
        {
            for (Vertex *vertex2 : cycles[1].vertices)
            {
                moves.push_back(make_pair(make_pair(vertex1, vertex2), nullptr));
            }
        }
        for (Graph &cycle : cycles)
        {
            for (Vertex *vertex1 : cycle.vertices)
            {
                for (Vertex *vertex2 : cycle.vertices)
                {
                    if (vertex1 != vertex2)
                    {
                        moves.push_back(make_pair(make_pair(vertex1, vertex2), &cycle));
                    }
                }
            }
        }

        // random_shuffle(moves.begin(), moves.end());

        pair<pair<Vertex *, Vertex *>, Graph *> bestMove;

        bestDelta = 0;

        for (auto move : moves)
        {
            Vertex *vertex1 = move.first.first;
            Vertex *vertex2 = move.first.second;
            Graph *graphToSwap = move.second;

            int delta = 0;

            if (vertex1->neighbours[0] == vertex2 || vertex1->neighbours[1] == vertex2)
            {
                if (vertex1->neighbours[0] == vertex2 && vertex2->neighbours[0] == vertex1)
                {
                    delta = distanceMatrix[vertex1->neighbours[1]->id][vertex2->id] + distanceMatrix[vertex2->neighbours[1]->id][vertex1->id] - distanceMatrix[vertex1->neighbours[1]->id][vertex1->id] - distanceMatrix[vertex2->neighbours[1]->id][vertex2->id];
                }
                else if (vertex1->neighbours[1] == vertex2 && vertex2->neighbours[1] == vertex1)
                {
                    delta = distanceMatrix[vertex1->neighbours[0]->id][vertex2->id] + distanceMatrix[vertex2->neighbours[0]->id][vertex1->id] - distanceMatrix[vertex1->neighbours[0]->id][vertex1->id] - distanceMatrix[vertex2->neighbours[0]->id][vertex2->id];
                }
                else if (vertex1->neighbours[0] == vertex2 && vertex2->neighbours[1] == vertex1)
                {
                    delta = distanceMatrix[vertex1->neighbours[1]->id][vertex2->id] + distanceMatrix[vertex2->neighbours[0]->id][vertex1->id] - distanceMatrix[vertex1->neighbours[1]->id][vertex1->id] - distanceMatrix[vertex2->neighbours[0]->id][vertex2->id];
                }
                else if (vertex1->neighbours[1] == vertex2 && vertex2->neighbours[0] == vertex1)
                {
                    delta = distanceMatrix[vertex1->neighbours[0]->id][vertex2->id] + distanceMatrix[vertex2->neighbours[1]->id][vertex1->id] - distanceMatrix[vertex1->neighbours[0]->id][vertex1->id] - distanceMatrix[vertex2->neighbours[1]->id][vertex2->id];
                }
            }
            else
            {
                delta = distanceMatrix[vertex1->neighbours[0]->id][vertex2->id] + distanceMatrix[vertex2->id][vertex1->neighbours[1]->id] + distanceMatrix[vertex2->neighbours[0]->id][vertex1->id] + distanceMatrix[vertex1->id][vertex2->neighbours[1]->id] - distanceMatrix[vertex1->neighbours[0]->id][vertex1->id] - distanceMatrix[vertex1->id][vertex1->neighbours[1]->id] - distanceMatrix[vertex2->neighbours[0]->id][vertex2->id] - distanceMatrix[vertex2->id][vertex2->neighbours[1]->id];
            }

            if (delta < bestDelta)
            {
                bestMove = move;
                bestDelta = delta;
            }
        }

        if (bestDelta < 0)
        {
            Vertex *vertex1 = bestMove.first.first;
            Vertex *vertex2 = bestMove.first.second;
            Graph *graphToSwap = bestMove.second;

            if (graphToSwap == nullptr)
            {
                swapVerticesBetweenCycles(vertex1->id, vertex2->id, cycles, distanceMatrix);
            }
            else
            {
                swapVerticesInCycle(vertex1->id, vertex2->id, *graphToSwap, distanceMatrix);
            }

            for(Graph &cycle : cycles)
            {
                cycle.normalizeGraph();
            }
        }
    } while(bestDelta < 0);

    return cycles;
}

vector<Graph> greedyLocalSearchNeighbourhood1(vector<Graph> &cycles, const vector<vector<int>> &distanceMatrix)
{
    int delta;

    do
    {
        vector<pair<pair<Vertex *, Vertex *>, Graph *>> moves;

        for (Vertex *vertex1 : cycles[0].vertices)
        {
            for (Vertex *vertex2 : cycles[1].vertices)
            {
                moves.push_back(make_pair(make_pair(vertex1, vertex2), nullptr));
            }
        }

        for (Graph &cycle : cycles)
        {
            for (Vertex *vertex1 : cycle.vertices)
            {
                for (Vertex *vertex2 : cycle.vertices)
                {
                    if (vertex1 != vertex2)
                    {
                        moves.push_back(make_pair(make_pair(vertex1, vertex2), &cycle));
                    }
                }
            }
        }

        random_shuffle(moves.begin(), moves.end());

        for (auto move : moves)
        {
            Vertex *vertex1 = move.first.first;
            Vertex *vertex2 = move.first.second;
            Graph *graphToSwap = move.second;

            delta = 0;

            if (vertex1->neighbours[0] == vertex2 && vertex2->neighbours[0] == vertex1)
            {
                delta = distanceMatrix[vertex1->neighbours[1]->id][vertex2->id] + distanceMatrix[vertex2->neighbours[1]->id][vertex1->id] - distanceMatrix[vertex1->neighbours[1]->id][vertex1->id] - distanceMatrix[vertex2->neighbours[1]->id][vertex2->id];
            }
            else if (vertex1->neighbours[1] == vertex2 && vertex2->neighbours[1] == vertex1)
            {
                delta = distanceMatrix[vertex1->neighbours[0]->id][vertex2->id] + distanceMatrix[vertex2->neighbours[0]->id][vertex1->id] - distanceMatrix[vertex1->neighbours[0]->id][vertex1->id] - distanceMatrix[vertex2->neighbours[0]->id][vertex2->id];
            }
            else if (vertex1->neighbours[0] == vertex2 && vertex2->neighbours[1] == vertex1)
            {
                delta = distanceMatrix[vertex1->neighbours[1]->id][vertex2->id] + distanceMatrix[vertex2->neighbours[0]->id][vertex1->id] - distanceMatrix[vertex1->neighbours[1]->id][vertex1->id] - distanceMatrix[vertex2->neighbours[0]->id][vertex2->id];
            }
            else if (vertex1->neighbours[1] == vertex2 && vertex2->neighbours[0] == vertex1)
            {
                delta = distanceMatrix[vertex1->neighbours[0]->id][vertex2->id] + distanceMatrix[vertex2->neighbours[1]->id][vertex1->id] - distanceMatrix[vertex1->neighbours[0]->id][vertex1->id] - distanceMatrix[vertex2->neighbours[1]->id][vertex2->id];
            }
            else
            {
                delta = distanceMatrix[vertex1->neighbours[0]->id][vertex2->id] + distanceMatrix[vertex2->id][vertex1->neighbours[1]->id] + distanceMatrix[vertex2->neighbours[0]->id][vertex1->id] + distanceMatrix[vertex1->id][vertex2->neighbours[1]->id] - distanceMatrix[vertex1->neighbours[0]->id][vertex1->id] - distanceMatrix[vertex1->id][vertex1->neighbours[1]->id] - distanceMatrix[vertex2->neighbours[0]->id][vertex2->id] - distanceMatrix[vertex2->id][vertex2->neighbours[1]->id];
            }

            if (delta < 0)
            {
                Vertex *vertex1 = move.first.first;
                Vertex *vertex2 = move.first.second;
                Graph *graphToSwap = move.second;

                if (graphToSwap == nullptr)
                {
                    swapVerticesBetweenCycles(vertex1->id, vertex2->id, cycles, distanceMatrix);
                }
                else
                {
                    swapVerticesInCycle(vertex1->id, vertex2->id, *graphToSwap, distanceMatrix);
                }

                for(Graph &cycle : cycles)
                {
                    cycle.normalizeGraph();
                }

                break;
            }
        }
    } while(delta < 0);

    return cycles;
}

vector<Graph> steepyLocalSearchNeighbourhood2(vector<Graph> &cycles, const vector<vector<int>> &distanceMatrix)
{
    int bestDelta = 0;

    do
    {
        vector<pair<variant<pair<Vertex *, Vertex *>, pair<Edge *, Edge *>>, Graph *>> moves;

        for (Vertex *vertex1 : cycles[0].vertices)
        {
            for (Vertex *vertex2 : cycles[1].vertices)
            {
                moves.push_back(make_pair(make_pair(vertex1, vertex2), nullptr));
            }
        }

        for (Graph &graph : cycles)
        {
            for (Edge *edge1 : graph.edges)
            {
                for (Edge *edge2 : graph.edges)
                {
                    if (edge1 != edge2)
                    {
                        if (edge1->src != edge2->dest && edge1->dest == edge2->src)
                        {
                            moves.push_back(make_pair(make_pair(edge1, edge2), &graph));
                            moves.push_back(make_pair(make_pair(edge1, edge2), &graph));
                        }
                    }
                }
            }
        }

        // random_shuffle(moves.begin(), moves.end());

        pair<variant<pair<Vertex *, Vertex *>, pair<Edge *, Edge *>>, Graph *> bestMove;

        bestDelta = 0;

        for (auto move : moves)
        {
            int delta = 0;

            if (move.second == nullptr)
            {
                Vertex *vertex1 = get<0>(move.first).first;
                Vertex *vertex2 = get<0>(move.first).second;
                Graph *graphToSwap = move.second;

                delta = distanceMatrix[vertex1->neighbours[0]->id][vertex2->id] + distanceMatrix[vertex2->id][vertex1->neighbours[1]->id] + distanceMatrix[vertex2->neighbours[0]->id][vertex1->id] + distanceMatrix[vertex1->id][vertex2->neighbours[1]->id] - distanceMatrix[vertex1->neighbours[0]->id][vertex1->id] - distanceMatrix[vertex1->id][vertex1->neighbours[1]->id] - distanceMatrix[vertex2->neighbours[0]->id][vertex2->id] - distanceMatrix[vertex2->id][vertex2->neighbours[1]->id];
            }
            else
            {
                Edge *edge1 = get<1>(move.first).first;
                Edge *edge2 = get<1>(move.first).second;
                Graph *graphToSwap = move.second;

                Vertex *vertex11 = edge1->src;
                Vertex *vertex12 = edge1->dest;
                Vertex *vertex21 = edge2->src;
                Vertex *vertex22 = edge2->dest;

                delta = distanceMatrix[vertex11->id][vertex21->id] + distanceMatrix[vertex12->id][vertex22->id] - edge1->distance - edge2->distance;
            }

            if (delta < bestDelta)
            {
                bestMove = move;
                bestDelta = delta;
            }
        }

        if (bestDelta < 0)
        {
            if (bestMove.second == nullptr)
            {
                Vertex *vertex1 = get<0>(bestMove.first).first;
                Vertex *vertex2 = get<0>(bestMove.first).second;

                swapVerticesBetweenCycles(vertex1->id, vertex2->id, cycles, distanceMatrix);
            }
            else
            {
                Edge *edge1 = get<1>(bestMove.first).first;
                Edge *edge2 = get<1>(bestMove.first).second;

                Graph *graphToSwap = bestMove.second;

                Vertex *vertex11 = edge1->src;
                Vertex *vertex12 = edge1->dest;
                Vertex *vertex21 = edge2->src;
                Vertex *vertex22 = edge2->dest;

                graphToSwap->removeEdge(vertex11, vertex12);
                graphToSwap->removeEdge(vertex21, vertex22);

                graphToSwap->addEdge(vertex11, vertex21, distanceMatrix[vertex11->id][vertex21->id]);
                graphToSwap->addEdge(vertex12, vertex22, distanceMatrix[vertex12->id][vertex22->id]);
            }

            for(Graph &cycle : cycles)
            {
                cycle.normalizeGraph();
            }
        }
    }while(bestDelta < 0);

    return cycles;
}

vector<Graph> greedyLocalSearchNeighbourhood2(vector<Graph> &cycles, const vector<vector<int>> &distanceMatrix)
{
    int delta;
    do
    {
        vector<pair<variant<pair<Vertex *, Vertex *>, pair<Edge *, Edge *>>, Graph *>> moves;

        for (Vertex *vertex1 : cycles[0].vertices)
        {
            for (Vertex *vertex2 : cycles[1].vertices)
            {
                moves.push_back(make_pair(make_pair(vertex1, vertex2), nullptr));
            }
        }

        for (Graph &graph : cycles)
        {
            for (Edge *edge1 : graph.edges)
            {
                for (Edge *edge2 : graph.edges)
                {
                    if (edge1 != edge2)
                    {
                        if (edge1->src != edge2->dest && edge1->dest == edge2->src)
                        {
                            moves.push_back(make_pair(make_pair(edge1, edge2), &graph));
                            moves.push_back(make_pair(make_pair(edge1, edge2), &graph));
                        }
                    }
                }
            }
        }

        random_shuffle(moves.begin(), moves.end());

        for (auto move : moves)
        {
            delta = 0;

            if (move.second == nullptr)
            {
                Vertex *vertex1 = get<0>(move.first).first;
                Vertex *vertex2 = get<0>(move.first).second;
                Graph *graphToSwap = move.second;

                delta = distanceMatrix[vertex1->neighbours[0]->id][vertex2->id] + distanceMatrix[vertex2->id][vertex1->neighbours[1]->id] + distanceMatrix[vertex2->neighbours[0]->id][vertex1->id] + distanceMatrix[vertex1->id][vertex2->neighbours[1]->id] - distanceMatrix[vertex1->neighbours[0]->id][vertex1->id] - distanceMatrix[vertex1->id][vertex1->neighbours[1]->id] - distanceMatrix[vertex2->neighbours[0]->id][vertex2->id] - distanceMatrix[vertex2->id][vertex2->neighbours[1]->id];
            }
            else
            {
                Edge *edge1 = get<1>(move.first).first;
                Edge *edge2 = get<1>(move.first).second;
                Graph *graphToSwap = move.second;

                Vertex *vertex11 = edge1->src;
                Vertex *vertex12 = edge1->dest;
                Vertex *vertex21 = edge2->src;
                Vertex *vertex22 = edge2->dest;

                delta = distanceMatrix[vertex11->id][vertex21->id] + distanceMatrix[vertex12->id][vertex22->id] - edge1->distance - edge2->distance;
            }

            if (delta < 0)
            {
                if (move.second == nullptr)
                {
                    Vertex *vertex1 = get<0>(move.first).first;
                    Vertex *vertex2 = get<0>(move.first).second;

                    swapVerticesBetweenCycles(vertex1->id, vertex2->id, cycles, distanceMatrix);
                }
                else
                {
                    Edge *edge1 = get<1>(move.first).first;
                    Edge *edge2 = get<1>(move.first).second;

                    Graph *graphToSwap = move.second;

                    Vertex *vertex11 = edge1->src;
                    Vertex *vertex12 = edge1->dest;
                    Vertex *vertex21 = edge2->src;
                    Vertex *vertex22 = edge2->dest;

                    graphToSwap->removeEdge(vertex11, vertex12);
                    graphToSwap->removeEdge(vertex21, vertex22);

                    graphToSwap->addEdge(vertex11, vertex21, distanceMatrix[vertex11->id][vertex21->id]);
                    graphToSwap->addEdge(vertex12, vertex22, distanceMatrix[vertex12->id][vertex22->id]);
                }

                for(Graph &cycle : cycles)
                {
                    cycle.normalizeGraph();
                }

                break;
            }
        }
    } while(delta < 0);

    return cycles;
}

vector<Graph> randomWalkNeighbourhood1(vector<Graph> &cycles, const vector<vector<int>> &distanceMatrix)
{
    vector<Graph> bestCycles = cycles;
    for (int i = 0; i < 350; i++)
    {
        vector<pair<pair<Vertex *, Vertex *>, Graph *>> moves;

        for (Vertex *vertex1 : cycles[0].vertices)
        {
            for (Vertex *vertex2 : cycles[1].vertices)
            {
                moves.push_back(make_pair(make_pair(vertex1, vertex2), nullptr));
            }
        }

        for (Graph &cycle : cycles)
        {
            for (Vertex *vertex1 : cycle.vertices)
            {
                for (Vertex *vertex2 : cycle.vertices)
                {
                    if (vertex1 != vertex2)
                    {
                        moves.push_back(make_pair(make_pair(vertex1, vertex2), &cycle));
                    }
                }
            }
        }

        random_shuffle(moves.begin(), moves.end());

        auto move = moves[0];
        Vertex *vertex1 = move.first.first;
        Vertex *vertex2 = move.first.second;
        Graph *graphToSwap = move.second;

        if (graphToSwap == nullptr)
        {
            swapVerticesBetweenCycles(vertex1->id, vertex2->id, cycles, distanceMatrix);
        }
        else
        {
            swapVerticesInCycle(vertex1->id, vertex2->id, *graphToSwap, distanceMatrix);
        }

        if(cycles[0].distance + cycles[1].distance < bestCycles[0].distance + bestCycles[1].distance)
        {
            bestCycles = cycles;
        }

        for(Graph &cycle : cycles)
        {
            cycle.normalizeGraph();
        }

    }

    return bestCycles;
}

vector<Graph> randomWalkNeighbourhood2(vector<Graph> &cycles, const vector<vector<int>> &distanceMatrix)
{
    vector<Graph> bestCycles = cycles;
    for (int i = 0; i < 350; i++)
    {
        vector<pair<variant<pair<Vertex *, Vertex *>, pair<Edge *, Edge *>>, Graph *>> moves;

        for (Vertex *vertex1 : cycles[0].vertices)
        {
            for (Vertex *vertex2 : cycles[1].vertices)
            {
                moves.push_back(make_pair(make_pair(vertex1, vertex2), nullptr));
            }
        }

        for (Graph &graph : cycles)
        {
            for (Edge *edge1 : graph.edges)
            {
                for (Edge *edge2 : graph.edges)
                {
                    if (edge1 != edge2)
                    {
                        if (edge1->src != edge2->dest && edge1->dest == edge2->src)
                        {
                            moves.push_back(make_pair(make_pair(edge1, edge2), &graph));
                            moves.push_back(make_pair(make_pair(edge1, edge2), &graph));
                        }
                    }
                }
            }
        }

        random_shuffle(moves.begin(), moves.end());

        auto move = moves[0];

        if (move.second == nullptr)
        {
            Vertex *vertex1 = get<0>(move.first).first;
            Vertex *vertex2 = get<0>(move.first).second;

            swapVerticesBetweenCycles(vertex1->id, vertex2->id, cycles, distanceMatrix);
        }
        else
        {
            Edge *edge1 = get<1>(move.first).first;
            Edge *edge2 = get<1>(move.first).second;

            Graph *graphToSwap = move.second;

            Vertex *vertex11 = edge1->src;
            Vertex *vertex12 = edge1->dest;
            Vertex *vertex21 = edge2->src;
            Vertex *vertex22 = edge2->dest;

            graphToSwap->removeEdge(vertex11, vertex12);
            graphToSwap->removeEdge(vertex21, vertex22);

            graphToSwap->addEdge(vertex11, vertex21, distanceMatrix[vertex11->id][vertex21->id]);
            graphToSwap->addEdge(vertex12, vertex22, distanceMatrix[vertex12->id][vertex22->id]);
        }

        for(Graph &cycle : cycles)
        {
            cycle.normalizeGraph();
        }

        if(cycles[0].distance + cycles[1].distance < bestCycles[0].distance + bestCycles[1].distance)
        {
            bestCycles = cycles;
        }
    }

    return bestCycles;
}

int main()
{
    vector<vector<int>> verticesCoords = readKroaFile("kroB100.tsp");

    vector<vector<int>> distanceMatrix = createDistanceMatrix(verticesCoords);

    // srand(time(NULL));

    pair<int, int> bestRandomValue = {INT_MAX, -1};
    pair<int, int> bestGreedyValue = {INT_MAX, -1};

    pair<int, int> worstRandomValue = {-1, -1};
    pair<int, int> worstGreedyValue = {-1, -1};

    long averageRandomValue = 0;
    long averageGreedyValue = 0;

    pair<int, int> bestRandomTime = {INT_MAX, -1};
    pair<int, int> bestGreedyTime = {INT_MAX, -1};

    pair<int, int> worstRandomTime = {-1, -1};
    pair<int, int> worstGreedyTime = {-1, -1};

    long averageRandomTime = 0;
    long averageGreedyTime = 0;

    vector<Graph> bestRandomCyclesResult;
    vector<Graph> bestGreedyCyclesResult;

    for(int i = 0; i < 100; i++){

        chrono::steady_clock::time_point beginRandomCycles = chrono::steady_clock::now();

        vector<Graph> randomCyclesStart = randomCycles(distanceMatrix);

        vector<Graph> randomCyclesResult = randomWalkNeighbourhood2(randomCyclesStart, distanceMatrix);

        chrono::steady_clock::time_point endRandomCycles = chrono::steady_clock::now();

        // for(Edge *e : randomCyclesResult[0].edges){
        //     cout << e->src->id << " " << e->dest->id << endl;
        // }

        int elapsedRandomCycles = chrono::duration_cast<chrono::milliseconds>(endRandomCycles - beginRandomCycles).count();


        if(randomCyclesResult[0].distance + randomCyclesResult[1].distance < bestRandomValue.first){
            bestRandomValue = {randomCyclesResult[0].distance + randomCyclesResult[1].distance, i};
            bestRandomCyclesResult = randomCyclesResult;
        }

        if(randomCyclesResult[0].distance + randomCyclesResult[1].distance > worstRandomValue.first){
            worstRandomValue = {randomCyclesResult[0].distance + randomCyclesResult[1].distance, i};
        }

        averageRandomValue += randomCyclesResult[0].distance + randomCyclesResult[1].distance;

        if(elapsedRandomCycles > worstRandomTime.first){
            worstRandomTime = {elapsedRandomCycles, i};
        }

        if(elapsedRandomCycles < bestRandomTime.first){
            bestRandomTime = {elapsedRandomCycles, i};
        }

        averageRandomTime += elapsedRandomCycles;




        chrono::steady_clock::time_point beginGreedyCycles = chrono::steady_clock::now();

        vector<Graph> greedyCyclesStart = greedyCycles(distanceMatrix, i);

        vector<Graph> greedyCyclesResult = randomWalkNeighbourhood2(greedyCyclesStart, distanceMatrix);

        chrono::steady_clock::time_point endGreedyCycles = chrono::steady_clock::now();

        int elapsedGreedyCycles = chrono::duration_cast<chrono::milliseconds>(endGreedyCycles - beginGreedyCycles).count();


        if(greedyCyclesResult[0].distance + greedyCyclesResult[1].distance < bestGreedyValue.first){
            bestGreedyValue = {greedyCyclesResult[0].distance + greedyCyclesResult[1].distance, i};
            bestGreedyCyclesResult = greedyCyclesResult;
        }

        if(greedyCyclesResult[0].distance + greedyCyclesResult[1].distance > worstGreedyValue.first){
            worstGreedyValue = {greedyCyclesResult[0].distance + greedyCyclesResult[1].distance, i};
        }

        averageGreedyValue += greedyCyclesResult[0].distance + greedyCyclesResult[1].distance;

        if(elapsedGreedyCycles > worstGreedyTime.first){
            worstGreedyTime = {elapsedGreedyCycles, i};
        }

        if(elapsedGreedyCycles < bestGreedyTime.first){
            bestGreedyTime = {elapsedGreedyCycles, i};
        }

        averageGreedyTime += elapsedGreedyCycles;

        cout << i << endl;
    }

    averageRandomValue /= 100;
    averageGreedyValue /= 100;

    averageRandomTime /= 100;
    averageGreedyTime /= 100;

    cout << "greedyLocalSearchNeighbourhood2" << endl;

    cout << "RandomCycles" << endl;
    cout << averageRandomValue << " (" << bestRandomValue.first << " – " << worstRandomValue.first << ")" << endl;
    cout << averageRandomTime << " (" << bestRandomTime.first << " – " << worstRandomTime.first << ")" << endl;

    cout << "GreedyCycles" << endl;
    cout << averageGreedyValue << " (" << bestGreedyValue.first << " – " << worstGreedyValue.first << ")" << endl;
    cout << averageGreedyTime << " (" << bestGreedyTime.first << " – " << worstGreedyTime.first << ")" << endl;

    saveGraphs(bestRandomCyclesResult, "randomWalkNeighbourhood2RandomCyclesB.txt");
    saveGraphs(bestGreedyCyclesResult, "randomWalkNeighbourhood2GreedyCyclesB.txt");

    return 0;
}