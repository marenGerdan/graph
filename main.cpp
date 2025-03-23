#include <iostream>
#include <vector>
#include <list>
#include <queue>
#include <stack>
#include <limits>
#include <cmath>
#include <unordered_map>
#include <set>
#include <iomanip>
#include <fstream>
#include <GL/glut.h>
#include <GL/glu.h>
#include <GL/gl.h>
#include <sstream>
#include <filesystem>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <climits>

const int WINDOW_WIDTH = 800;
const int WINDOW_HEIGHT = 800;
const int RADIUS = 5;
const int OFFSET = 50;

class Graph
{
private:
    int vertices;
    std::vector<std::vector<int>> capacity;
    std::vector<std::vector<int>> adjMatrix;
    std::vector<std::list<std::pair<int, int>>> adjList;
    std::vector<std::pair<GLfloat, GLfloat>> positions;
    std::vector<std::pair<int, int>> highlightedEdges;

    double heuristic(int a, int b)
    {
        return std::hypot(positions[a].first - positions[b].first,
                          positions[a].second - positions[b].second);
    }

    void drawCircle(GLfloat x, GLfloat y)
    {
        glBegin(GL_TRIANGLE_FAN);
        glVertex2f(x, y);
        for (int i = 0; i <= 72; i++)
        {
            float angle = 2.0f * M_PI * i / 72;
            glVertex2f(x + RADIUS * cos(angle), y + RADIUS * sin(angle));
        }
        glEnd();
    }

    void DFSUtil(int v, std::vector<bool> &visited, std::vector<std::pair<int, int>> &dfsEdges)
    {
        visited[v] = true;
        std::cout << v << " ";

        std::vector<std::pair<int, int>> neighbors(adjList[v].begin(), adjList[v].end());

        std::sort(neighbors.begin(), neighbors.end(),
                  [](const std::pair<int, int> &a, const std::pair<int, int> &b)
                  {
                      return a.first < b.first;
                  });

        for (const auto &[neighbor, _] : neighbors)
        {
            if (!visited[neighbor])
            {
                dfsEdges.emplace_back(v, neighbor); // додати ребро
                DFSUtil(neighbor, visited, dfsEdges);
            }
        }
    }

public:
    explicit Graph(int v) : vertices(v),
                            adjMatrix(v, std::vector<int>(v, 0)),
                            capacity(v, std::vector<int>(v, 0)),
                            adjList(v)
    {
        // Початкове випадкове розміщення вершин
        for (int i = 0; i < v; ++i)
        {
            positions.emplace_back(
                OFFSET + std::rand() % (WINDOW_WIDTH - 2 * OFFSET),
                OFFSET + std::rand() % (WINDOW_HEIGHT - 2 * OFFSET));
        }

        // Застосування алгоритму Fruchterman-Reingold для покращення розташування вершин
        arrangeLayeredBFS(0); // Починаємо з вершини 0

    }

    void loadFromFile(const std::string &filename)
    {
        std::ifstream file(filename);
        if (!file)
        {
            std::cerr << "Error opening file: " << filename << std::endl;
            return;
        }

        int v, e;
        file >> v >> e;
        if (v > vertices)
        {
            std::cerr << "Error: File contains more vertices than initialized graph!" << std::endl;
            return;
        }

        int node = 0;
        std::string line;
        std::getline(file, line);
        while (std::getline(file, line) && node < vertices)
        {
            std::istringstream iss(line);
            int neighbor;
            while (iss >> neighbor)
            {
                if (neighbor >= 0 && neighbor < vertices)
                {
                    addEdge(node, neighbor);
                }
                if (neighbor >= vertices)
                {
                    std::cerr << "Invalid neighbor index: " << neighbor << " for node " << node << std::endl;
                    continue;
                }
            }
            node++;
        }
    }

    void addEdge(int u, int v, int weight = 1)
    {
        adjMatrix[u][v] = weight;
        adjMatrix[v][u] = weight;
        adjList[u].emplace_back(v, weight);
        adjList[v].emplace_back(u, weight);

        // Ініціалізуємо пропускну здатність для Ford-Fulkerson
        capacity[u][v] = weight;
        capacity[v][u] = weight; // Для неорієнтованого графа
    }


    void BFS(int start)
    {
        std::vector<bool> visited(vertices, false);
        std::queue<int> q;
        std::vector<std::pair<int, int>> bfsEdges;

        std::cout << "BFS: ";

        q.push(start);
        visited[start] = true;

        while (!q.empty())
        {
            int v = q.front();
            q.pop();
            std::cout << v << " ";

            std::vector<std::pair<int, int>> neighbors(adjList[v].begin(), adjList[v].end());

            std::sort(neighbors.begin(), neighbors.end(),
                      [](const std::pair<int, int> &a, const std::pair<int, int> &b)
                      {
                          return a.first < b.first;
                      });

            for (const auto &[neighbor, _] : neighbors)
            {
                if (!visited[neighbor])
                {
                    visited[neighbor] = true;
                    q.push(neighbor);
                    bfsEdges.emplace_back(v, neighbor); // додати ребро
                }
            }
        }

        // підсвітити всі ребра, які були використані в BFS
        highlightEdges(bfsEdges);
        std::cout << std::endl;
    }

    void DFS(int start)
    {
        std::vector<bool> visited(vertices, false);
        std::vector<std::pair<int, int>> dfsEdges;
        std::cout << "DFS: ";
        DFSUtil(start, visited, dfsEdges);
        highlightEdges(dfsEdges);
        std::cout << std::endl;
    }

    void primMST()
    {
        std::vector<int> parent(vertices, -1);
        std::vector<int> key(vertices, INT_MAX);
        std::vector<bool> inMST(vertices, false);
        std::vector<std::pair<int, int>> mstEdges;
        key[0] = 0; // Почати з першої вершини

        for (int count = 0; count < vertices - 1; ++count)
        {
            int u = -1;

            // Знайти вершину з мінімальним ключем, яка ще не в MST
            for (int v = 0; v < vertices; ++v)
            {
                if (!inMST[v] && (u == -1 || key[v] < key[u]))
                {
                    u = v;
                }
            }

            inMST[u] = true;

            // Оновлюємо ключі сусідів вибраної вершини
            for (auto &[v, weight] : adjList[u])
            {
                if (!inMST[v] && weight < key[v])
                {
                    parent[v] = u;
                    key[v] = weight;
                    mstEdges.emplace_back(u, v); // Додаємо ребро до MST
                }
            }
        }

        // Виведення ребер MST у консоль
        std::cout << "Prim's MST: \n";
        for (const auto &[u, v] : mstEdges)
        {
            std::cout << u << " - " << v << std::endl;
        }

        // Підсвітити ребра MST
        highlightEdges(mstEdges);
    }

    void dijkstra(int start)
    {
        std::vector<int> dist(vertices, std::numeric_limits<int>::max());
        dist[start] = 0;
        std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<>> pq;
        pq.emplace(0, start);
        std::vector<std::pair<int, int>> dijkstraEdges;

        while (!pq.empty())
        {
            int u = pq.top().second;
            pq.pop();

            for (const auto &[v, weight] : adjList[u])
            {
                if (dist[u] + weight < dist[v])
                {
                    dist[v] = dist[u] + weight;
                    pq.emplace(dist[v], v);
                    dijkstraEdges.emplace_back(u, v); // додати ребро
                }
            }
        }

        std::cout << "Shortest Paths (Dijkstra's Algorithm) from " << start << ":\n";
        for (int i = 0; i < vertices; ++i)
        {
            std::cout << "To " << i << ": " << dist[i] << "\n";
        }

        // Підсвітити всі ребра, які були використані в Dijkstra
        highlightEdges(dijkstraEdges);
    }

    void aStar(int start, int goal)
    {
        std::vector<int> parent(vertices, -1);
        std::vector<double> gScore(vertices, std::numeric_limits<double>::max());
        std::vector<double> fScore(vertices, std::numeric_limits<double>::max());
        gScore[start] = 0;
        fScore[start] = heuristic(start, goal);

        std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<>> openSet;
        openSet.emplace(fScore[start], start);
        std::vector<std::pair<int, int>> aStarEdges;

        while (!openSet.empty())
        {
            int current = openSet.top().second;
            openSet.pop();

            if (current == goal)
            {
                std::cout << "A* Path: ";
                int node = goal;
                while (node != -1)
                {
                    std::cout << node << " ";
                    node = parent[node];
                }
                std::cout << std::endl;
                // Підсвітити ребра A* шляху
                highlightEdges(aStarEdges);
                return;
            }

            for (const auto &[neighbor, weight] : adjList[current])
            {
                double tentative_gScore = gScore[current] + weight;
                if (tentative_gScore < gScore[neighbor])
                {
                    parent[neighbor] = current;
                    gScore[neighbor] = tentative_gScore;
                    fScore[neighbor] = gScore[neighbor] + heuristic(neighbor, goal);
                    openSet.emplace(fScore[neighbor], neighbor);
                    aStarEdges.emplace_back(current, neighbor); // додати ребро
                }
            }
        }

        std::cout << "No path found using A*" << std::endl;
    }

    int bfsFordFulkerson(int source, int sink, std::vector<int> &parent)
    {
        std::fill(parent.begin(), parent.end(), -1);
        std::vector<int> maxFlow(vertices, 0);
        std::vector<bool> visited(vertices, false);

        // Пріоритетна черга: {максимальний потік до вершини, вершина}
        std::priority_queue<std::pair<int, int>> pq;
        pq.push({INT_MAX, source});
        maxFlow[source] = INT_MAX;
        visited[source] = true;

        while (!pq.empty())
        {
            int flow = pq.top().first;
            int u = pq.top().second;
            pq.pop();

            if (u == sink)
                return flow; // Досягли призначення, повертаємо максимальний потік для цього шляху

            for (int v = 0; v < vertices; v++)
            {
                if (!visited[v] && capacity[u][v] > 0)
                {
                    parent[v] = u;
                    int new_flow = std::min(flow, capacity[u][v]); // Мінімальне ребро у шляху
                    pq.push({new_flow, v});
                    maxFlow[v] = new_flow;
                    visited[v] = true;
                }
            }
        }

        return 0; // Якщо не знайшли шлях
    }

    int fordFulkerson(int source, int sink)
    {
        if (source == sink)
            return 0;

        int max_flow = 0;
        std::vector<int> parent(vertices);
        std::vector<std::pair<int, int>> bestPath;

        auto original_capacity = capacity; // Копія початкової матриці потужностей

        while (int flow = bfsFordFulkerson(source, sink, parent))
        {
            max_flow += flow;
            bestPath.clear();

            // Відновлення шляху та оновлення залишкової місткості
            for (int v = sink; v != source; v = parent[v])
            {
                int u = parent[v];
                capacity[u][v] -= flow;
                capacity[v][u] += flow;
                bestPath.push_back({u, v});
            }
        }

        capacity = original_capacity; // Відновлення початкових потужностей після виконання алгоритму
        highlightedEdges = bestPath;

        std::cout << "\nОптимальний шлях із максимальною пропускною здатністю (" << max_flow << "):\n";
        for (auto &edge : bestPath)
            std::cout << edge.first << " -> " << edge.second << "\n";

        glutPostRedisplay(); // Оновлення візуалізації
        return max_flow;
    }

    int bidirectionalDijkstra(int start, int goal)
    {
        std::vector<int> distFwd(vertices, std::numeric_limits<int>::max());
        std::vector<int> distBwd(vertices, std::numeric_limits<int>::max());
        std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<>> pqFwd, pqBwd;
        distFwd[start] = 0;
        distBwd[goal] = 0;
        pqFwd.emplace(0, start);
        pqBwd.emplace(0, goal);

        std::set<int> visitedFwd, visitedBwd;
        std::vector<std::pair<int, int>> bidirectionalEdges;

        while (!pqFwd.empty() && !pqBwd.empty())
        {
            // Forward search
            int u = pqFwd.top().second;
            pqFwd.pop();
            visitedFwd.insert(u);

            for (const auto &[v, weight] : adjList[u])
            {
                if (distFwd[u] + weight < distFwd[v])
                {
                    distFwd[v] = distFwd[u] + weight;
                    pqFwd.emplace(distFwd[v], v);
                    bidirectionalEdges.push_back({u, v}); // Додаємо ребро для підсвічування
                    highlightEdges(bidirectionalEdges);   // Підсвічуємо після кожного кроку
                }
            }

            // Backward search
            int w = pqBwd.top().second;
            pqBwd.pop();
            visitedBwd.insert(w);

            for (const auto &[v, weight] : adjList[w])
            {
                if (distBwd[w] + weight < distBwd[v])
                {
                    distBwd[v] = distBwd[w] + weight;
                    pqBwd.emplace(distBwd[v], v);
                    bidirectionalEdges.push_back({w, v}); // Додаємо ребро для підсвічування
                    highlightEdges(bidirectionalEdges);   // Підсвічуємо після кожного кроку
                }
            }

            // Перевіряємо, чи є спільні відвідані вершини
            if (visitedFwd.count(w) || visitedBwd.count(u))
            {
                highlightEdges(bidirectionalEdges); // Підсвічуємо відповідні ребра
                return distFwd[u] + distBwd[w];     // Повертаємо мінімальний шлях
            }
        }

        highlightEdges(bidirectionalEdges);     // Якщо не знайдено шляхів
        return std::numeric_limits<int>::max(); // Повертаємо невідому відстань
    }

    int bidirectionalAStar(int start, int goal)
    {
        std::vector<int> distFwd(vertices, std::numeric_limits<int>::max());
        std::vector<int> distBwd(vertices, std::numeric_limits<int>::max());
        std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<>> pqFwd, pqBwd;
        distFwd[start] = 0;
        distBwd[goal] = 0;
        pqFwd.emplace(heuristic(start, goal), start);
        pqBwd.emplace(heuristic(goal, start), goal);

        std::set<int> visitedFwd, visitedBwd;
        std::vector<std::pair<int, int>> bidirectionalEdges;

        while (!pqFwd.empty() && !pqBwd.empty())
        {
            // Forward search
            int u = pqFwd.top().second;
            pqFwd.pop();
            visitedFwd.insert(u);

            for (const auto &[v, weight] : adjList[u])
            {
                if (distFwd[u] + weight < distFwd[v])
                {
                    distFwd[v] = distFwd[u] + weight;
                    pqFwd.emplace(distFwd[v] + heuristic(v, goal), v);
                    bidirectionalEdges.push_back({u, v}); // Додаємо ребро для підсвічування
                }
            }

            // Backward search
            int w = pqBwd.top().second;
            pqBwd.pop();
            visitedBwd.insert(w);

            for (const auto &[v, weight] : adjList[w])
            {
                if (distBwd[w] + weight < distBwd[v])
                {
                    distBwd[v] = distBwd[w] + weight;
                    pqBwd.emplace(distBwd[v] + heuristic(v, start), v);
                    bidirectionalEdges.push_back({w, v}); // Додаємо ребро для підсвічування
                }
            }

            // Перевіряємо, чи є спільні відвідані вершини
            if (visitedFwd.count(w) || visitedBwd.count(u))
            {
                highlightEdges(bidirectionalEdges); // Підсвічуємо відповідні ребра
                return distFwd[u] + distBwd[w];     // Повертаємо мінімальний шлях
            }
        }

        highlightEdges(bidirectionalEdges);     // Якщо не знайдено шляхів
        return std::numeric_limits<int>::max(); // Повертаємо невідому відстань
    }

    void drawGraph()
    {
        glClear(GL_COLOR_BUFFER_BIT);

        glColor3f(0.7f, 0.7f, 0.7f);
        glBegin(GL_LINES);
        for (int u = 0; u < vertices; ++u)
        {
            for (const auto &[v, _] : adjList[u])
            {
                glVertex2f(positions[u].first, positions[u].second);
                glVertex2f(positions[v].first, positions[v].second);
            }
        }
        glEnd();

        glColor3f(1.0f, 0.0f, 0.0f);
        glBegin(GL_LINES);
        for (const auto &[u, v] : highlightedEdges)
        {
            glVertex2f(positions[u].first, positions[u].second);
            glVertex2f(positions[v].first, positions[v].second);
        }
        glEnd();

        // Малюємо вершини
        glColor3f(0.3f, 0.6f, 1.0f);
        for (int i = 0; i < vertices; ++i)
        {
            drawCircle(positions[i].first, positions[i].second);
        }

        // Відображаємо номери вершин
        glColor3f(1.0f, 1.0f, 1.0f);
        for (int i = 0; i < vertices; ++i)
        {
            float textX = std::min(std::max(positions[i].first + RADIUS, 5.0f), WINDOW_WIDTH - 15.0f);
            float textY = std::min(std::max(positions[i].second + RADIUS, 5.0f), WINDOW_HEIGHT - 15.0f);
            glRasterPos2f(positions[i].first + RADIUS, positions[i].second + RADIUS);
            std::string label = std::to_string(i);
            for (char c : label)
            {
                glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, c);
            }
        }

        glutSwapBuffers();
    }

    void drawEdge(int u, int v)
    {
        GLfloat x1 = positions[u].first;
        GLfloat y1 = positions[u].second;
        GLfloat x2 = positions[v].first;
        GLfloat y2 = positions[v].second;
        glBegin(GL_LINES);
        glVertex2f(x1, y1);
        glVertex2f(x2, y2);
        glEnd();
    }

    void highlightEdges(const std::vector<std::pair<int, int>> &edges)
    {
        highlightedEdges = edges;
        glutPostRedisplay();
    }

    void arrangeLayeredBFS(int start = 0)
    {
        std::vector<int> level(vertices, -1);
        std::queue<int> q;

        level[start] = 0;
        q.push(start);

        int maxLevel = 0;
        while (!q.empty())
        {
            int v = q.front();
            q.pop();

            for (const auto &[neighbor, _] : adjList[v])
            {
                if (level[neighbor] == -1) // Ще не відвідано
                {
                    level[neighbor] = level[v] + 1;
                    maxLevel = std::max(maxLevel, level[neighbor]);
                    q.push(neighbor);
                }
            }
        }

        // Визначаємо кількість вершин на кожному рівні
        std::vector<int> countPerLevel(maxLevel + 1, 0);
        for (int l : level)
        {
            if (l != -1)
                countPerLevel[l]++;
        }

        // Розташування вершин на рівнях
        std::vector<int> offsetPerLevel(maxLevel + 1, 0);
        for (int i = 1; i <= maxLevel; i++)
        {
            offsetPerLevel[i] = offsetPerLevel[i - 1] + countPerLevel[i - 1];
        }

        // Розміщуємо вершини в сітковому форматі
        const float levelSpacing = (WINDOW_HEIGHT - 2 * OFFSET) / (maxLevel + 1);
        std::vector<int> used(vertices, 0);

        for (int i = 0; i < vertices; i++)
        {
            if (level[i] != -1)
            {
                float x = OFFSET + (WINDOW_WIDTH - 2 * OFFSET) * (used[level[i]] + 1) / (countPerLevel[level[i]] + 1);
                float y = WINDOW_HEIGHT - OFFSET - level[i] * levelSpacing;
                positions[i] = {x, y};
                used[level[i]]++;
            }
        }
    }
};

Graph *g = nullptr;
// Graph *e = nullptr;

void display()
{
    g->drawGraph();
    // e->drawGraph();
}

struct WindowData
{
    Graph *graph;
    std::string name;
    int windowId;
};

std::vector<WindowData> windows;

void printSeparator()
{
    std::cout << "\n========================================\n";
}

// Функція відображення для конкретного вікна
void displayWindow(void)
{
    int windowId = glutGetWindow();
    for (const auto &window : windows)
    {
        if (window.windowId == windowId)
        {
            window.graph->drawGraph();
            break;
        }
    }
}

void initGL()
{
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, WINDOW_WIDTH, 0, WINDOW_HEIGHT);
}

void createWindow(Graph *graph, const std::string &name, int x, int y)
{
    glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
    glutInitWindowPosition(x, y);
    int windowId = glutCreateWindow(name.c_str());

    initGL();

    WindowData window;
    window.graph = graph;
    window.name = name;
    window.windowId = windowId;
    windows.push_back(window);

    glutDisplayFunc(displayWindow);
}

void runAllAlgorithms(const std::vector<Graph *> &graphs)
{
    // BFS
    std::cout << "\nRunning BFS from node 0:";
    printSeparator();
    graphs[0]->BFS(0);

    // DFS
    std::cout << "\nRunning DFS from node 0:";
    printSeparator();
    graphs[1]->DFS(0);

    // Prim's MST
    std::cout << "\nRunning Prim's MST Algorithm:";
    printSeparator();
    graphs[2]->primMST();

    // Dijkstra
    std::cout << "\nRunning Dijkstra's Algorithm from node 0:";
    printSeparator();
    graphs[3]->dijkstra(0);

    // A*
    std::cout << "\nRunning A* Algorithm from node 0 to node 3:";
    printSeparator();
    graphs[4]->aStar(0, 3);

    // Ford-Fulkerson
    std::cout << "\nRunning Ford-Fulkerson Algorithm:";
    printSeparator();
    int maxFlow = graphs[5]->fordFulkerson(0, 4);
    std::cout << "Maximum Flow: " << maxFlow << std::endl;

    // Bidirectional Dijkstra
    std::cout << "\nRunning Bidirectional Dijkstra from node 0 to node 50:";
    printSeparator();
    int biDijkstraDistance = graphs[6]->bidirectionalDijkstra(0, 50);
    std::cout << "Shortest path length: " << biDijkstraDistance << std::endl;

    // Bidirectional A*
    std::cout << "\nRunning Bidirectional A* from node 0 to node 50:";
    printSeparator();
    int biAstarDistance = graphs[7]->bidirectionalAStar(0, 50);
    std::cout << "Shortest path length: " << biAstarDistance << std::endl;
}

int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);

    // Створюємо вектор для зберігання всіх графів
    std::vector<Graph *> graphs;

    // Створюємо окремі графи для кожного алгоритму
    graphs.push_back(new Graph(199)); // BFS
    graphs.push_back(new Graph(199)); // DFS
    graphs.push_back(new Graph(199)); // Prim
    graphs.push_back(new Graph(199)); // Dijkstra
    graphs.push_back(new Graph(199)); // A*
    graphs.push_back(new Graph(199)); // Ford-Fulkerson
    graphs.push_back(new Graph(199)); // Bi-Dijkstra
    graphs.push_back(new Graph(199)); // Bi-A*

    // Завантажуємо дані в кожен граф
    std::string filename = "/home/sasha/Документи/suka/jazz.graph";
    for (auto *graph : graphs)
    {
        graph->loadFromFile(filename);
    }

    // Створюємо вікна для кожного алгоритму
    std::vector<std::string> windowNames = {
        "BFS Visualization",
        "DFS Visualization",
        "Prim's MST Visualization",
        "Dijkstra Visualization",
        "A* Visualization",
        "Ford-Fulkerson Visualization",
        "Bidirectional Dijkstra",
        "Bidirectional A*"};

    // Створюємо вікна в сітці 4x2
    for (size_t i = 0; i < graphs.size(); ++i)
    {
        int x = (i % 4) * WINDOW_WIDTH;
        int y = (i / 4) * WINDOW_HEIGHT;
        createWindow(graphs[i], windowNames[i], x, y);
    }

    std::cout << "Starting algorithm visualizations...\n";
    printSeparator();

    // Запускаємо всі алгоритми з виведенням в консоль
    runAllAlgorithms(graphs);

    glutMainLoop();

    // Очищення пам'яті
    for (auto *graph : graphs)
    {
        delete graph;
    }

    return 0;
}
