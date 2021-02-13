#include <iostream>
#include "csv.h"
#include <vector>
#include "string"
#include <cmath>
#include <algorithm>
#include <tuple>
#include <chrono>
#include "omp.h"
#include "string"
#include <set>

using namespace io;

std::vector<std::string> tests = {
    "_5",
    "_8",
    "_10",
    "_14",
    "_15",
    // "_16",
};

std::vector<std::vector<float>> coords_from_csv(std::string filepath)
{
    CSVReader<2, trim_chars<' '>, no_quote_escape<';'>> in(filepath);
    float x;
    float y;
    std::vector<float> city_coords;
    std::vector<std::vector<float>> coords;
    while (in.read_row(x, y))
    {
        city_coords = {x, y};
        coords.emplace_back(city_coords);
    }
    return coords;
}

std::vector<int> load_results(std::string filepath)
{
    CSVReader<1, trim_chars<>> in(filepath);
    int x;
    std::vector<int> results;
    while (in.read_row(x))
    {
        results.emplace_back(x);
    }
    return results;
}

float euclidean_distance(std::vector<float> coords1, std::vector<float> coords2)
{
    float x1 = coords1[0];
    float y1 = coords1[1];

    float x2 = coords2[0];
    float y2 = coords2[1];

    float distance = std::sqrt(std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2));
    return distance;
}

std::vector<std::vector<float>> get_distance_matrix(std::vector<std::vector<float>> coords)
{
    auto n_cities = coords.size();
    std::vector<std::vector<float>> matrix;
    std::cout << "Count of cities: " << n_cities << std::endl;
    for (int i = 0; i < n_cities; i++)
    {
        matrix.emplace_back();
        for (int j = 0; j < n_cities; j++)
        {
            if (i != j)
            {
                auto distance = euclidean_distance(coords[i], coords[j]);
                matrix[i].push_back(distance);
            }
            else
            {
                matrix[i].push_back(0.0);
            }
        }
    }

    return matrix;
}

int find_min_idx(const std::vector<float> &vec)
{
    size_t min_idx = 0;
    size_t size = vec.size();
    float min_value = INFINITY;
    float curr_val;

    for (size_t i = 0; i < size; i++)
    {
        curr_val = vec[i];
        if (curr_val < min_value)
        {
            min_idx = i;
            min_value = curr_val;
        }
    }
    return min_idx;
}

// Функция нахождения кратчайщего пути с учетом
// заданного начального перехода 0->first_step_destination
std::tuple<float, std::vector<int>> find_shortest(const std::vector<std::vector<float>> matrix,
                                                  int first_step_destination)
{
    auto n_cities = matrix.size();

    std::vector<int> valid_cities;
    for (int i = 1; i < n_cities; ++i)
    {
        if (i != first_step_destination)
        {
            valid_cities.push_back(i);
        }
    }

    std::vector<int> best_path;
    int prev_city;
    float min_distance = INFINITY;
    float path_distance;

    do
    {
        path_distance = matrix[0][first_step_destination]; // исходное расстояние
        prev_city = first_step_destination;
        for (int city_idx : valid_cities)
        {
            path_distance += matrix[prev_city][city_idx];
            prev_city = city_idx;
        }

        path_distance += matrix[0][prev_city];
        if (path_distance < min_distance)
        {
            min_distance = path_distance;
            best_path = valid_cities;
        }

    } while (std::next_permutation(valid_cities.begin(), valid_cities.end()));

    // добавление в лучший путь начальных точек
    best_path.insert(best_path.begin(), first_step_destination);
    best_path.insert(best_path.begin(), 0);
    return {min_distance, best_path};
}

std::vector<int> test(std::string test_num)
{
    std::string filename = "tests\\test" + test_num + ".csv";
    auto coords = coords_from_csv(filename);
    std::vector<std::vector<float>> matrix = get_distance_matrix(coords);

    std::vector<std::vector<int>> paths;
    std::vector<float> distances;
    size_t matrix_size = matrix.size();
    for (int i = 0; i < matrix_size - 1; i++)
    {
        distances.emplace_back(0);
        paths.emplace_back(0);
    }

    auto start = std::chrono::high_resolution_clock::now();
    auto wall_timer = omp_get_wtime();

    int i;
    float distance;
    std::vector<int> path;
#pragma omp parallel for
    for (i = 0; i < matrix_size - 1; i++)
    {
        std::tie(distance, path) = find_shortest(matrix, i + 1);
        distances[i] = distance;
        paths[i] = path;
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    int min_idx = find_min_idx(std::ref(distances));
    auto min_distance = distances[min_idx];
    auto best_path = paths[min_idx];

    std::cout << "MIN DIST: " << min_distance << std::endl;
    std::cout << "PATH: ";
    for (auto i : best_path)
    {
        std::cout << i << " ";
    }

    std::cout << "\nTime taken by function: "
              << omp_get_wtime() - wall_timer << " Seconds\n";
    return best_path;
}

int main()
{
    for (const auto &filenum : tests)
    {
        auto path = test(filenum);
        std::string filename = "tests\\result" + filenum + ".csv";
        auto resultCsv = load_results(filename);
        if (path.empty() || (resultCsv.size() > path.size()))
        {
            std::cout << "test" << test << ": FAIL\n";
            continue;
        }
        std::set<std::pair<int, int>> edges;
        auto prev = path.back();
        for (int i = 0; i < resultCsv.size(); ++i)
        {
            edges.emplace(prev, path[i]);
            prev = path[i];
        }
        auto ok = true;
        prev = resultCsv[resultCsv.size() - 1];
        for (int i = 0; i < resultCsv.size(); ++i)
        {
            if (edges.find(std::pair<int, int>(prev, resultCsv[i])) == edges.end() &&
                edges.find(std::pair<int, int>(resultCsv[i], prev)) == edges.end())
            {
                std::cout << "test" << filenum << ": FAIL\n";
                ok = false;
                break;
            }
            prev = resultCsv[i];
        }
        if (ok)
        {
            std::cout << "test" << filenum << ": OK\n\n";
        }
    }

    return 0;
}