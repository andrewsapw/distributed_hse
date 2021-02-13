#include <iostream>
#include "csv.h"
#include <vector>
#include "string"
#include <algorithm>
#include <tuple>
#include <chrono>
#include "omp.h"
#include <set>
#include <atomic>
#include <stdio.h>
#include <cstdlib>

using namespace io;

// прочитать данные из .csv
std::tuple<std::vector<float>, std::vector<float>> read_csv(std::string filepath)
{
    CSVReader<2, trim_chars<' '>, no_quote_escape<';'>> in(filepath);
    float weight;
    float price;
    std::vector<float> prices;
    std::vector<float> weigths;
    std::vector<std::vector<float>> res;
    while (in.read_row(weight, price))
    {
        weigths.emplace_back(weight);
        prices.emplace_back(price);
    }
    return {weigths, prices};
}

/*
 * Нахождение лучшнй комбинации предметов.
 * Возвращает лучшую комбинацию.
 * 
 * Аргументы:
 * N : общее количество предметов.
 * prices : вектор цен предметов.
 * weigths : вектор весов предметов.
 * capacity : вместимость рюкзака.
*/
std::vector<int> solve(int N, std::vector<float> prices, std::vector<float> weigths, int capacity)
{
    float best_value = 0;
    std::vector<int> best_comb;
    std::atomic<size_t> counter{0};
#pragma omp parallel for schedule(static)
    for (int K = 1; K < N; ++K)
    {
        float best_K_value = 0;
        std::vector<int> best_K_comb;

        std::string bitmask(K, 1);
        bitmask.resize(N, 0);

        // все комбинации
        do
        {
            float current_weigth = 0; // вес комбинации предметов, shape=(K)
            float current_value = 0;  // стоимость комбинации предметов
            std::vector<int> current_comb;
            counter++;

            // для каждого предмета из комбинации
            for (int i = 0; i < N; ++i)
            {
                if (bitmask[i])
                {
                    current_weigth += weigths[i];

                    // если вес учтенных на данный момент предметов
                    // превышает вместимость - пропускаем
                    if (current_weigth > capacity)
                    {
                        break;
                    }
                    current_value += prices[i];
                    current_comb.push_back(i);
                }
            }

            if (current_value > best_K_value)
            {
                best_K_value = current_value;
                best_K_comb = current_comb;
            }

        } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

        // обновляем лучшее цену
        if (best_K_value > best_value)
        {
            best_value = best_K_value;
            best_comb = best_K_comb;
        }
    }

    std::cout << "Best value: " << best_value << " " << counter << std::endl;

    return best_comb;
}

// сумма элементов в векторе
template <typename T>
T sum(std::vector<T> vec)
{
    float result_sum = 0;
    for (auto x : vec)
    {
        result_sum += x;
    }
    return result_sum;
}

std::vector<int> read_solution(const int test_num)
{
    std::string filepath = "tests/bpresult_" + std::to_string(test_num) + ".csv";
    LineReader in(filepath);

    std::vector<int> solution;
    int counter = 0;
    while (char *line = in.next_line())
    {
        auto line_length = strlen(line);

        for (int i = 0; i < line_length; ++i)
        {
            if (line[i] != ';')
            {
                if (line[i] == '1')
                {
                    solution.push_back(counter);
                }
                counter += 1;
            }
        }
    }
    return solution;
}

void check_solution(std::vector<int> best_path)
{
    return;
}

int main(int argc, char const *argv[])
{
    if (argc < 3)
    {
        printf("Not enough arguments");
        return 1;
    }

    int N = atoi(argv[1]);
    int N_THREADS = atoi(argv[2]);
    omp_set_num_threads(N_THREADS);

    std::vector<float> prices;
    std::vector<float> weigths;
    std::string file_num = std::to_string(N);

    std::string filename = "tests\\test_" + file_num + ".csv";
    std::tie(weigths, prices) = read_csv(filename);

    float sum_weigths = sum(weigths);
    int capacity = static_cast<int>(sum_weigths / 2);

    auto wall_timer = omp_get_wtime(); // начало отсчета времени
    std::vector<int> best_comb = solve(N, prices, weigths, capacity);
    auto elapsed_time = omp_get_wtime() - wall_timer;

    int weight = 0;
    std::cout << "Capacity: " << capacity << " Num threads: " << N_THREADS << std::endl;
    printf_s("Wiegth: %d\n", weight);
    printf_s("Time taken by function: %f seconds", elapsed_time);

    printf("\nBest comb: ");
    for (auto x : best_comb)
    {
        weight += weigths[x];
        printf_s("%d ", x);
    }
    printf("\nSolution: ");
    auto solution = read_solution(N);
    for (auto x : solution)
    {
        std::cout << x << ' ';
    }
    if (solution == best_comb)
    {
        printf_s("\nTest %d: OK", N);
    }
    else
    {
        printf_s("\nTest %d: FAIL", N);
    }

    return 0;
}
