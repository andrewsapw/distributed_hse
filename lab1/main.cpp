#include <iostream>
#include "csv.h"
#include <vector>
#include "string"
#include <cmath>
#include <algorithm>
#include <tuple>
#include <chrono>
#include "omp.h"

using namespace io;

std::vector<std::vector<double>> coords_from_csv(std::string filepath)
{
    CSVReader<2, trim_chars<' '>, no_quote_escape<';'>> in(filepath);
    double x;
    double y;

    std::vector<std::vector<double>> coords;
    while (in.read_row(x, y))
    {
        std::vector<double> city_coords = {x, y};
        coords.emplace_back(city_coords);
    }

    return coords;
}

double euclidean_distance(std::vector<double> coords1, std::vector<double> coords2)
{
    double x1 = coords1[0];
    double y1 = coords1[1];

    double x2 = coords2[0];
    double y2 = coords2[1];

    double distance = std::sqrt(std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2));
    return distance;
}

std::vector<std::vector<double>> get_distance_matrix(std::vector<std::vector<double>> coords)
{
    auto n_cities = coords.size();
    std::vector<std::vector<double>> matrix;
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

int find_min_idx(const std::vector<double> &vec)
{
    size_t min_idx = 0;
    size_t size = vec.size();
    double min_value = INFINITY;
    double curr_val;

    for (int i = 0; i < size; i++)
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

double cals_path_distance(std::vector<int> path, const std::vector<std::vector<double>> &matrix, int first_point)
{
    double path_distance, distance = 0;
    auto prev_city = first_point;
    for (auto city_idx : path)
    {
        distance = matrix[city_idx][prev_city];
        path_distance += distance;
        prev_city = city_idx;
    }
    return path_distance;
}

std::tuple<double, std::vector<int>> find_shortest(const std::vector<std::vector<double>> matrix,
                                                   int first_step_destination, int &count)
{
    auto n_cities = matrix.size();

    std::vector<int> valid_cities;
    for (int i = 0; i < n_cities; i++)
    {
        if (i != first_step_destination)
        {
            valid_cities.push_back(i);
        }
    }

    double path_distance;
    std::vector<int> best_path;
    bool useless_path;
    int prev_city;
    int idx = 0;
    double min_distance = INFINITY;

    do
    {
        if (valid_cities.front() < valid_cities.back())
        {
            continue;
        }
#pragma omp atomic
        count++;
        path_distance = 0.0;
        useless_path = false;

        prev_city = first_step_destination;
        for (int city_idx : valid_cities)
        {
            path_distance += matrix[prev_city][city_idx];
            prev_city = city_idx;

            if (path_distance >= min_distance)
            {
                useless_path = true;
                break;
            }
        }
        if (!useless_path)
        {
            path_distance += matrix[prev_city][first_step_destination];
            if (path_distance < min_distance)
            {
                min_distance = path_distance;
                best_path = valid_cities;
            }
        }

    } while (std::next_permutation(valid_cities.begin(), valid_cities.end()));

    best_path.insert(best_path.begin(), first_step_destination);
    return {min_distance, best_path};
}

std::tuple<double, std::vector<int>> find_shortest_2(const std::vector<std::vector<double>> matrix,
                                                     int first_step_destination, int &count)
{
    auto n_cities = matrix.size();

    std::vector<int> valid_cities;
    for (int i = 1; i < n_cities; i++)
    {
        if (i != first_step_destination)
        {
            valid_cities.push_back(i);
        }
    }

    double path_distance;
    std::vector<int> best_path;
    bool useless_path;
    int prev_city;
    int idx = 0;
    double min_distance = INFINITY;

    do
    {
#pragma omp atomic
        count++;
        path_distance = matrix[0][first_step_destination];
        useless_path = false;

        prev_city = first_step_destination;
        for (int city_idx : valid_cities)
        {
            path_distance += matrix[prev_city][city_idx];
            prev_city = city_idx;

            if (path_distance >= min_distance)
            {
                useless_path = true;
                break;
            }
        }
        if (!useless_path)
        {
            // path_distance += matrix[prev_city][first_step_destination];
            if (path_distance < min_distance)
            {
                min_distance = path_distance;
                best_path = valid_cities;
            }
        }

    } while (std::next_permutation(valid_cities.begin(), valid_cities.end()));

    best_path.insert(best_path.begin(), first_step_destination);
    best_path.insert(best_path.begin(), 0);
    return {min_distance, best_path};
}

int main()
{
    auto coords = coords_from_csv("tests\\test_14.csv");
    std::vector<std::vector<double>> matrix = get_distance_matrix(coords);

    std::vector<std::vector<int>> paths;
    std::vector<double> distances;
    double min_distance = INFINITY;
    int count;

    size_t matrix_size = matrix.size();
    for (int i = 0; i < matrix_size - 1; i++)
    {
        distances.emplace_back();
        paths.emplace_back();
    }
    auto start = std::chrono::high_resolution_clock::now();

#pragma omp parallel for num_threads(4) shared(matrix, count)
    for (size_t i = 0; i < matrix_size - 1; i++)
    {
        auto [distance, path] = find_shortest_2(matrix, i + 1, std::ref(count));
        distances[i] = distance;
        paths[i] = path;
    }

    // #pragma omp parallel for num_threads(4) shared(matrix, count)
    //     for (size_t i = 0; i < matrix_size; i++)
    //     {
    //         auto [distance, path] = find_shortest(matrix, i, std::ref(count));
    //         distances[i] = distance;
    //         paths[i] = path;
    //     }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    int min_idx = find_min_idx(std::ref(distances));
    min_distance = distances[min_idx];
    auto best_path = paths[min_idx];

    std::cout << "MIN DIST: " << min_distance << std::endl;
    std::cout << "PATH: ";
    for (auto i : best_path)
    {
        std::cout << i << " ";
    }

    std::cout << "\nTime taken by function: "
              << duration.count() << " milliseconds"
              << "(" << count << ")" << std::endl;
}