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
    Vertex* prev;
    Vertex* next;

    Vertex(int _id) : id(_id), prev(nullptr), next(nullptr) {}
};

class Edge
{
public:
    Vertex *src;
    Vertex *dest;
    int distance;

    Edge(Vertex *_src, Vertex *_dest, int _distance) : src(_src), dest(_dest), distance(_distance) {}

    void remove(){
        delete this;
    }
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
        try{
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
                src->next = dest;
                dest->prev = src;
                this->distance += distance;
                edges.push_back(e);
            }
        }
        catch(exception &e){
            cout << e.what() << endl;
        }
    }

    void removeVertex(Vertex *v)
    {
        if (v)
        {
            removeEdge(v->prev, v);
            removeEdge(v, v->next);

            vertices.erase(remove(vertices.begin(), vertices.end(), v), vertices.end());
        }
    }

    void removeEdge(Vertex *src, Vertex *dest)
    {
        if (src && dest)
        {
            for (Edge *e : edges)
            {
                if (e->src == src && e->dest == dest)
                {
                    this->distance -= e->distance;
                    edges.erase(remove(edges.begin(), edges.end(), e), edges.end());
                    // e->remove();
                    break;
                }
            }
            src->next = nullptr;
            dest->prev = nullptr;
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

    vector<Edge*> findPath(Vertex *start, Vertex *end)
    {
        vector<Edge*> path;
        Vertex *currentVertex = start;

        do
        {
            for (Edge *e : edges)
            {
                if (e->src == currentVertex)
                {
                    path.push_back(e);
                    currentVertex = e->dest;
                    break;
                }
            }
        } while (currentVertex != end);

        return path;
    }
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
        cycle.addEdge(nearestNeighbour, vertex);
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

        // if (cycle.vertices.size() > 3)
        // {
            cycle.removeEdge(minEdge->src, minEdge->dest);
        // }

        cycle.addEdge(minEdge->src, newVertex, distanceMatrix[minEdge->src->id][vertexId]);
        cycle.addEdge(newVertex, minEdge->dest, distanceMatrix[newVertex->id][minEdge->dest->id]);

        if (cycle.vertices.size() == 3)
        {
            cycle.addEdge(minEdge->dest, minEdge->src, minEdge->distance);
        }

        visited[vertexId] = true;
    }

    return cycles;
}

void swapVerticesBetweenCycles(int i, int j, vector<Graph> &graphs, const vector<vector<int>> &distanceMatrix)
{
    Graph *graph1 = &graphs[0];
    Graph *graph2 = &graphs[1];

    Vertex *vertex1 = graph1->findVertex(i);
    Vertex *vertex2 = graph2->findVertex(j);

    Vertex *vertex1Prev = vertex1->prev;
    Vertex *vertex1Next = vertex1->next;

    Vertex *vertex2Prev = vertex2->prev;
    Vertex *vertex2Next = vertex2->next;

    graph1->removeVertex(vertex1);

    graph2->removeVertex(vertex2);

    graph1->addVertex(vertex2);
    graph1->addEdge(vertex1Prev, vertex2, distanceMatrix[vertex1Prev->id][vertex2->id]);
    graph1->addEdge(vertex2, vertex1Next, distanceMatrix[vertex2->id][vertex1Next->id]);

    graph2->addVertex(vertex1);
    graph2->addEdge(vertex2Prev, vertex1, distanceMatrix[vertex2Prev->id][vertex1->id]);
    graph2->addEdge(vertex1, vertex2Next, distanceMatrix[vertex1->id][vertex2Next->id]);
}

void swapVerticesInCycle(int i, int j, Graph &graph, const vector<vector<int>> &distanceMatrix)
{
    Vertex *vertex1 = graph.findVertex(i);
    Vertex *vertex2 = graph.findVertex(j);

    Vertex *vertex1Prev = vertex1->prev;
    Vertex *vertex1Next = vertex1->next;

    Vertex *vertex2Prev = vertex2->prev;
    Vertex *vertex2Next = vertex2->next;

    if (vertex1Prev == vertex2 || vertex1Next == vertex2)
    {
        if (vertex1Prev == vertex2 && vertex2Prev == vertex1)
        {
            graph.removeEdge(vertex1Next, vertex1);
            graph.removeEdge(vertex2, vertex2Next);

            graph.addEdge(vertex1, vertex2Next, distanceMatrix[vertex1->id][vertex2Next->id]);
            graph.addEdge(vertex1Next, vertex2, distanceMatrix[vertex1Next->id][vertex2->id]);
        }
        else if (vertex1Next == vertex2 && vertex2Prev == vertex1)
        {
            graph.removeEdge(vertex1Prev, vertex1);
            graph.removeEdge(vertex2, vertex2Next);

            graph.addEdge(vertex1, vertex2Next, distanceMatrix[vertex1->id][vertex2Next->id]);
            graph.addEdge(vertex1Prev, vertex2, distanceMatrix[vertex1Prev->id][vertex2->id]);
        }
        else if (vertex1Prev == vertex2 && vertex2Next == vertex1)
        {
            graph.removeEdge(vertex1Next, vertex1);
            graph.removeEdge(vertex2Prev, vertex2);

            graph.addEdge(vertex1, vertex2Prev, distanceMatrix[vertex1->id][vertex2Prev->id]);
            graph.addEdge(vertex1Next, vertex2, distanceMatrix[vertex1Next->id][vertex2->id]);
        }
        else if (vertex1Next == vertex2 && vertex2Next == vertex1)
        {
            graph.removeEdge(vertex1Prev, vertex1);
            graph.removeEdge(vertex2Prev, vertex2);

            graph.addEdge(vertex1, vertex2Prev, distanceMatrix[vertex1->id][vertex2Prev->id]);
            graph.addEdge(vertex1Prev, vertex2, distanceMatrix[vertex1Prev->id][vertex2->id]);
        }
    }
    else
    {
        graph.removeEdge(vertex1Prev, vertex1);
        graph.removeEdge(vertex1, vertex1Next);
        graph.removeEdge(vertex2Prev, vertex2);
        graph.removeEdge(vertex2, vertex2Next);

        graph.addEdge(vertex1Prev, vertex2, distanceMatrix[vertex1Prev->id][vertex2->id]);
        graph.addEdge(vertex2, vertex1Next, distanceMatrix[vertex2->id][vertex1Next->id]);
        graph.addEdge(vertex2Prev, vertex1, distanceMatrix[vertex2Prev->id][vertex1->id]);
        graph.addEdge(vertex1, vertex2Next, distanceMatrix[vertex1->id][vertex2Next->id]);
    }
}

vector<Graph> steepestLocalSearchNeighbourhood1(vector<Graph> &cycles, const vector<vector<int>> &distanceMatrix)
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

            if (vertex1->prev == vertex2 || vertex1->next == vertex2)
            {
                if (vertex1->prev == vertex2 && vertex2->prev == vertex1)
                {
                    delta = distanceMatrix[vertex1->next->id][vertex2->id] + distanceMatrix[vertex2->next->id][vertex1->id] - distanceMatrix[vertex1->next->id][vertex1->id] - distanceMatrix[vertex2->next->id][vertex2->id];
                }
                else if (vertex1->next == vertex2 && vertex2->next == vertex1)
                {
                    delta = distanceMatrix[vertex1->prev->id][vertex2->id] + distanceMatrix[vertex2->prev->id][vertex1->id] - distanceMatrix[vertex1->prev->id][vertex1->id] - distanceMatrix[vertex2->prev->id][vertex2->id];
                }
                else if (vertex1->prev == vertex2 && vertex2->next == vertex1)
                {
                    delta = distanceMatrix[vertex1->next->id][vertex2->id] + distanceMatrix[vertex2->prev->id][vertex1->id] - distanceMatrix[vertex1->next->id][vertex1->id] - distanceMatrix[vertex2->prev->id][vertex2->id];
                }
                else if (vertex1->next == vertex2 && vertex2->prev == vertex1)
                {
                    delta = distanceMatrix[vertex1->prev->id][vertex2->id] + distanceMatrix[vertex2->next->id][vertex1->id] - distanceMatrix[vertex1->prev->id][vertex1->id] - distanceMatrix[vertex2->next->id][vertex2->id];
                }
            }
            else
            {
                delta = distanceMatrix[vertex1->prev->id][vertex2->id] + distanceMatrix[vertex2->id][vertex1->next->id] + distanceMatrix[vertex2->prev->id][vertex1->id] + distanceMatrix[vertex1->id][vertex2->next->id] - distanceMatrix[vertex1->prev->id][vertex1->id] - distanceMatrix[vertex1->id][vertex1->next->id] - distanceMatrix[vertex2->prev->id][vertex2->id] - distanceMatrix[vertex2->id][vertex2->next->id];
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

            if (vertex1->prev == vertex2 && vertex2->prev == vertex1)
            {
                delta = distanceMatrix[vertex1->next->id][vertex2->id] + distanceMatrix[vertex2->next->id][vertex1->id] - distanceMatrix[vertex1->next->id][vertex1->id] - distanceMatrix[vertex2->next->id][vertex2->id];
            }
            else if (vertex1->next == vertex2 && vertex2->next == vertex1)
            {
                delta = distanceMatrix[vertex1->prev->id][vertex2->id] + distanceMatrix[vertex2->prev->id][vertex1->id] - distanceMatrix[vertex1->prev->id][vertex1->id] - distanceMatrix[vertex2->prev->id][vertex2->id];
            }
            else if (vertex1->prev == vertex2 && vertex2->next == vertex1)
            {
                delta = distanceMatrix[vertex1->next->id][vertex2->id] + distanceMatrix[vertex2->prev->id][vertex1->id] - distanceMatrix[vertex1->next->id][vertex1->id] - distanceMatrix[vertex2->prev->id][vertex2->id];
            }
            else if (vertex1->next == vertex2 && vertex2->prev == vertex1)
            {
                delta = distanceMatrix[vertex1->prev->id][vertex2->id] + distanceMatrix[vertex2->next->id][vertex1->id] - distanceMatrix[vertex1->prev->id][vertex1->id] - distanceMatrix[vertex2->next->id][vertex2->id];
            }
            else
            {
                delta = distanceMatrix[vertex1->prev->id][vertex2->id] + distanceMatrix[vertex2->id][vertex1->next->id] + distanceMatrix[vertex2->prev->id][vertex1->id] + distanceMatrix[vertex1->id][vertex2->next->id] - distanceMatrix[vertex1->prev->id][vertex1->id] - distanceMatrix[vertex1->id][vertex1->next->id] - distanceMatrix[vertex2->prev->id][vertex2->id] - distanceMatrix[vertex2->id][vertex2->next->id];
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

                break;
            }
        }
    } while(delta < 0);

    return cycles;
}

vector<Graph> steepestLocalSearchNeighbourhood2(vector<Graph> &cycles, const vector<vector<int>> &distanceMatrix)
{
    int bestDelta = 0;

    int i = 0;

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
                        Vertex *vertex11 = edge1->src;
                        Vertex *vertex12 = edge1->dest;
                        Vertex *vertex21 = edge2->src;
                        Vertex *vertex22 = edge2->dest;

                        if (vertex11->id != vertex12->id && vertex11->id != vertex21->id && vertex11->id != vertex22->id && vertex12->id != vertex21->id && vertex12->id != vertex22->id && vertex21->id != vertex22->id)
                        {
                            moves.push_back(make_pair(make_pair(edge1, edge2), &graph));
                        }
                    }
                }
            }
        }

        pair<variant<pair<Vertex *, Vertex *>, pair<Edge *, Edge *>>, Graph *> bestMove;

        bestDelta = 0;

        for (auto move : moves)
        {
            int delta = 0;

            if (move.second == nullptr)
            {
                Vertex *vertex1 = get<0>(move.first).first;
                Vertex *vertex2 = get<0>(move.first).second;

                delta = distanceMatrix[vertex1->prev->id][vertex2->id] + distanceMatrix[vertex2->id][vertex1->next->id] + distanceMatrix[vertex2->prev->id][vertex1->id] + distanceMatrix[vertex1->id][vertex2->next->id] - distanceMatrix[vertex1->prev->id][vertex1->id] - distanceMatrix[vertex1->id][vertex1->next->id] - distanceMatrix[vertex2->prev->id][vertex2->id] - distanceMatrix[vertex2->id][vertex2->next->id];
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

                vector <Edge*> pathToReverse = graphToSwap->findPath(vertex12, vertex21);

                for(Edge *edge: pathToReverse){
                    swap(edge->src, edge->dest);
                    swap(edge->dest->next, edge->dest->prev);
                    if(edge->src->next == nullptr){
                        swap(edge->src->next, edge->src->prev);
                    }
                }

                graphToSwap->addEdge(vertex11, vertex21, distanceMatrix[vertex11->id][vertex21->id]);
                graphToSwap->addEdge(vertex12, vertex22, distanceMatrix[vertex12->id][vertex22->id]);
            }

        }

        i++;
    }
    while(bestDelta < 0);

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
                        Vertex *vertex11 = edge1->src;
                        Vertex *vertex12 = edge1->dest;
                        Vertex *vertex21 = edge2->src;
                        Vertex *vertex22 = edge2->dest;

                        if (vertex11->id != vertex12->id && vertex11->id != vertex21->id && vertex11->id != vertex22->id && vertex12->id != vertex21->id && vertex12->id != vertex22->id && vertex21->id != vertex22->id)
                        {
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

                delta = distanceMatrix[vertex1->prev->id][vertex2->id] + distanceMatrix[vertex2->id][vertex1->next->id] + distanceMatrix[vertex2->prev->id][vertex1->id] + distanceMatrix[vertex1->id][vertex2->next->id] - distanceMatrix[vertex1->prev->id][vertex1->id] - distanceMatrix[vertex1->id][vertex1->next->id] - distanceMatrix[vertex2->prev->id][vertex2->id] - distanceMatrix[vertex2->id][vertex2->next->id];
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

                    vector <Edge*> pathToReverse = graphToSwap->findPath(vertex12, vertex21);

                    for(Edge *edge: pathToReverse){
                        swap(edge->src, edge->dest);
                        swap(edge->dest->next, edge->dest->prev);
                        if(edge->src->next == nullptr){
                            swap(edge->src->next, edge->src->prev);
                        }
                    }

                    graphToSwap->addEdge(vertex11, vertex21, distanceMatrix[vertex11->id][vertex21->id]);
                    graphToSwap->addEdge(vertex12, vertex22, distanceMatrix[vertex12->id][vertex22->id]);
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
                        Vertex *vertex11 = edge1->src;
                        Vertex *vertex12 = edge1->dest;
                        Vertex *vertex21 = edge2->src;
                        Vertex *vertex22 = edge2->dest;

                        if (vertex11->id != vertex12->id && vertex11->id != vertex21->id && vertex11->id != vertex22->id && vertex12->id != vertex21->id && vertex12->id != vertex22->id && vertex21->id != vertex22->id)
                        {
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

            vector <Edge*> pathToReverse = graphToSwap->findPath(vertex12, vertex21);

            for(Edge *edge: pathToReverse){
                swap(edge->src, edge->dest);
                swap(edge->dest->next, edge->dest->prev);
                if(edge->src->next == nullptr){
                    swap(edge->src->next, edge->src->prev);
                }
            }

            graphToSwap->addEdge(vertex11, vertex21, distanceMatrix[vertex11->id][vertex21->id]);
            graphToSwap->addEdge(vertex12, vertex22, distanceMatrix[vertex12->id][vertex22->id]);
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

        vector<Graph> randomCyclesResult = steepestLocalSearchNeighbourhood2(randomCyclesStart, distanceMatrix);

        chrono::steady_clock::time_point endRandomCycles = chrono::steady_clock::now();

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

        // for(Graph &graph: greedyCyclesStart){
        //     for(Edge *edge: graph.edges){
        //         cout << edge->src->id << " " << edge->dest->id << endl;
        //     }
        //     cout << endl;
        // }

        vector<Graph> greedyCyclesResult = steepestLocalSearchNeighbourhood2(greedyCyclesStart, distanceMatrix);

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

    cout << "SteepestLocalSearchNeighbourhood2" << endl;

    cout << "RandomCycles" << endl;
    cout << averageRandomValue << " (" << bestRandomValue.first << " – " << worstRandomValue.first << ")" << endl;
    cout << averageRandomTime << " (" << bestRandomTime.first << " – " << worstRandomTime.first << ")" << endl;

    cout << "GreedyCycles" << endl;
    cout << averageGreedyValue << " (" << bestGreedyValue.first << " – " << worstGreedyValue.first << ")" << endl;
    cout << averageGreedyTime << " (" << bestGreedyTime.first << " – " << worstGreedyTime.first << ")" << endl;

    // saveGraphs(bestRandomCyclesResult, "randomWalkNeighbourhood2RandomCyclesB.txt");
    // saveGraphs(bestGreedyCyclesResult, "randomWalkNeighbourhood2GreedyCyclesB.txt");

    return 0;
}