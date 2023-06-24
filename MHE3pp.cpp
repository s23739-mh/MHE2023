#include <vector>
#include <iostream>
#include <numeric>
#include <random>
#include <math.h>
#include <list>
#include <functional>
#include <iostream>
#include <string>
#include <fstream>

std::mt19937 rgen(109908093657865);
using problem_t = std::vector<int>;



class Solution : public std::vector<int>
{

public:
    std::shared_ptr<problem_t> problem;

    Solution(std::shared_ptr<problem_t> _problem) : problem(_problem) {
    }

    auto tripletAverage() {
        return std::accumulate(this->begin(), this->end(), 0.0) * 3 / this->size();
    }


    auto score() -> double {
        std::vector<int> tripletSumsDiff;

        auto sum = 0.0;

        for (auto i = 0; i < size(); i += 3) {
            sum += (abs((at(i) + at(i + 1) + at(i + 2) - tripletAverage())));
        }

        return sum * 3 / size();
    }


    void print() {
        std::cout << "\n{";
        for (auto i = 0; i < size(); i++) {
            std::cout << at(i) << ',';
        }
        std::cout << "}\n";
        for (auto i = 0; i < size(); i += 3) {
            std::cout << " sum[" << i / 3 << "]:" << at(i) + at(i + 1) + at(i + 2);
        }
        std::cout << std::endl;
    }
};



auto generateNeighbors(Solution current) -> std::vector<Solution> {
    std::vector<Solution> neighbors;
    auto maximum = current.score();
    for (auto i = 0; i < current.size(); i++) {
        for (auto j = 0; j < current.size(); j++) {
            auto copy = current;
            std::swap(copy[i], copy[j]);
            if (copy.score() < maximum) {
                neighbors.push_back(copy);
            }
        }
    }
    return neighbors;
}

auto findMinimumElement(std::vector<Solution> solutions) {

    auto lowest = solutions.front().score();
    auto lowest_sol = solutions.front();
    for (auto solution : solutions) {
        auto solutionScore = solution.score();
        if (solutionScore < lowest) {
            lowest = solutionScore;
            lowest_sol = solution;
        }
    }
    return lowest_sol;
}



auto tabu_search(Solution current, int max_size) -> Solution {
    auto lowest = current.score();
    auto lowest_sol = current;
    std::list<Solution> tabu_list = { current };
    std::vector<Solution> neighbors;
    auto iter = 0;
    do {
        neighbors = generateNeighbors(tabu_list.back());
        for (auto tabu_item : tabu_list) {
            auto found = std::find(neighbors.begin(), neighbors.end(), tabu_item);
            if (found != neighbors.end()) {
                neighbors.erase(found);
                if (neighbors.size() == 0) {
                    return lowest_sol;
                }
            }

            if (neighbors.empty()) { return lowest_sol; }
        }
        
        auto lowest_neighbor = findMinimumElement(neighbors);
        if (lowest_neighbor.score() < lowest) {
            lowest = lowest_neighbor.score();
            lowest_sol = lowest_neighbor;
        }
        if (tabu_list.size() < max_size) { tabu_list.push_back(lowest_neighbor); }
        else {
            return lowest_sol;
        }
        
        iter++;
    } while (!neighbors.empty()&&iter<540);
    return lowest_sol;
}

Solution better_modify(Solution current) {
    auto neighbors = generateNeighbors(current);
    if (neighbors.empty()) { return current; }
    return findMinimumElement(neighbors);
}

Solution random_modify(Solution current, std::mt19937 rgen) {
    std::uniform_int_distribution<int> distr(0, current.size() - 1);
    int a = distr(rgen);
    std::swap(current[a], current[(a + 1) % current.size()]);
    return current;
}

Solution sim_annealing(Solution current, std::function<double(int)> T, std::mt19937 rgen) {
    //for
    auto lowest_sol = current;
    auto current_sol = current;
    for (auto i = 0; i < 1000; i++) {
        auto new_sol = random_modify(lowest_sol, rgen);
        if (new_sol.score() <= current_sol.score()) {
            current_sol = new_sol;
            if (new_sol.score() <= lowest_sol.score()) {
                lowest_sol = current_sol;
            }
        }
        else
        {
            std::uniform_real_distribution<double> u(0.0, 1.0);
            if (u(rgen) < std::exp(-std::abs(new_sol.score() - current_sol.score()) / T(i))) {
                current_sol = new_sol;
            }

        }

    }
    return lowest_sol;
}

problem_t generate_problem(int size, std::mt19937 rgen) {
    std::uniform_int_distribution<int> number(0, 50);
    problem_t problem;
    for (auto i = 0; i < size; i++) {
        problem.push_back(number(rgen));
    }
    return problem;
}

Solution random_solution(problem_t problem, std::mt19937 rgen) {
    Solution solution(std::make_shared<problem_t>(problem));
    for (int i = 0; i < problem.size(); i++) {
        solution.push_back(problem.at(i));
    }
    std::shuffle(solution.begin(), solution.end(), rgen);
    return solution;
}


Solution brute_force(Solution current) {
    do
    {
        if (current.score() == 0) { return current; }
    } while (std::next_permutation(current.begin(), current.end()));
    return current;
}

template<class T>
class generic_algorithm_config_t {
public:
    int population_size = 100;
    T sample_solution=0;
    std::string crossoverConf = "randomCross"; //pointCross
    std::string mutationConf = "worstMut"; //randomMut
    std::string conditionConf = "maxPopul"; //goalAchieved

    std::mt19937 rgen;
    virtual bool terminal_condition1() = 0;
    virtual bool terminal_condition2(std::vector<T>) = 0;
    virtual std::vector<T> get_initial_population() = 0;
    virtual double fitness(T) = 0;
    virtual std::vector<T> selection(std::vector<double>, std::vector<T>, std::mt19937) = 0;
    virtual std::vector<T> crossover(std::vector<T>, std::mt19937) = 0;
    virtual std::pair<Solution,Solution> crossover1(Solution, Solution , std::mt19937) = 0;
    virtual std::pair<Solution, Solution> crossover2(Solution, Solution, std::mt19937) = 0;
    virtual std::vector<T> mutation1(std::vector<T>, std::mt19937) = 0;
    virtual std::vector<T> mutation2(std::vector<T>, std::mt19937) = 0;
};

template<class T>
Solution generic_algorithm(generic_algorithm_config_t<T>& config, std::mt19937 rgen) {
    auto population = config.get_initial_population();

    auto condition = true;
    while (condition) {
        std::vector<double> fitnesses;
        for (auto i = 0; i < population.size(); i++) {
            fitnesses.push_back(config.fitness(population[i]));
        }
        auto parents = config.selection(fitnesses, population,rgen);
        std::vector<Solution> offspring;
        
        offspring = config.crossover(parents,rgen);
        
        if (config.mutationConf == "worstMut") {
            offspring = config.mutation1(offspring,rgen);
        }
        else {
            offspring = config.mutation2(offspring,rgen);
        }
        
        population = offspring;
        if (config.conditionConf == "maxPopul") {
            condition = config.terminal_condition1();
        }
        else {
            condition = config.terminal_condition2(population);
        }
        fitnesses.clear();
    }
    return *std::max_element(population.begin(), population.end(), [&](T l, T r) {return config.fitness(l) > config.fitness(r); });
}

class tsp_config : public generic_algorithm_config_t<Solution> {
public:
    int iteration=0;
    int max_iterations;
    problem_t problem;
    
    
    tsp_config(int iter, int pop_size, std::string crossover, std::string mutation, std::string condition, problem_t problem_, std::mt19937 rgenn) {
        iteration = 0;
        max_iterations = iter;
        population_size = pop_size;
        crossoverConf = crossover;
        mutationConf = mutation;
        conditionConf = condition;
        problem = problem_;
        rgen = rgenn;
    }
    

    virtual bool terminal_condition1() {
        iteration++;
        return iteration <= max_iterations;
    }

    virtual bool terminal_condition2(std::vector<Solution> sols) {
        for (auto i : sols) {
            if (fitness(i) == 1) {
                return true;
            }
        }
        return false;
    }

    virtual std::vector<Solution> get_initial_population() {
        std::vector<Solution> vect;
        for (auto i = 0; i < population_size; i++) {
            vect.push_back(random_solution(problem, rgen));
        }
        return vect;
    }
    virtual double fitness(Solution solution) {
        return 1 / (1 + solution.score());
    }

    virtual std::vector<Solution> selection(std::vector<double> fitnesses, std::vector<Solution> population, std::mt19937 rgen) {
        std::vector<Solution> vect;
        while (vect.size() < population.size()) {
            std::uniform_int_distribution<int> dist(0, population.size() - 1);
            int a_idx = dist(rgen);
            int b_idx = dist(rgen);
            if (fitnesses[a_idx] >= fitnesses[b_idx]) {
                vect.push_back(population[a_idx]);
            }
            else {
                vect.push_back(population[b_idx]);
            }
        }
        return vect;
    }
    std::pair<Solution, Solution> crossover1(Solution a, Solution b, std::mt19937 rgen) {
        std::uniform_int_distribution<int> distr(0, 1);
        for (auto j = 0; j < a.size(); j++) {
            if (distr(rgen)) {
                std::swap(a[j], b[j]);
            }
        }
        return { a,b };
    }
       
    std::pair<Solution, Solution> crossover2(Solution a, Solution b, std::mt19937 rgen) {
        std::uniform_int_distribution<int> pivot(1, b.size() - 2);
        int piv = pivot(rgen);
        for (auto i = 0; i < piv; i++) {
            std::swap(a[i], b[i]);
        }
        return { a,b };
    }           

    virtual std::vector<Solution> crossover(std::vector<Solution> pop, std::mt19937 rgen) {

        std::vector<Solution> offspring;
        for (auto i = 0; i < pop.size() - 1; i += 2) {
            if (crossoverConf == "randomCross") {
                auto ab = crossover1(pop.at(i), pop.at(i + 1), rgen);
                offspring.push_back(ab.first);
                offspring.push_back(ab.second);
            }
            else {
                auto ab = crossover2(pop.at(i), pop.at(i + 1), rgen);
                offspring.push_back(ab.first);
                offspring.push_back(ab.second);
            }
            
        }
        return offspring;
    }

    virtual std::vector<Solution> mutation1(std::vector<Solution> solution, std::mt19937 rgen) {
        std::vector<Solution> pop;
        for (auto sol : solution) {
            std::uniform_int_distribution<int> distr(0, 1);
            if (distr(rgen)) {
                sol = random_modify(sol, rgen);
            }
            pop.push_back(sol);
        }
        return pop;
    }

    virtual std::vector<Solution> mutation2(std::vector<Solution> solution, std::mt19937 rgen) {
        std::vector<Solution> pop;
        std::vector<double> scores;
        for (auto i = 0; i < solution.size(); i++) {
            scores.push_back(fitness(solution[i]));
            pop.push_back(solution[i]);
        }
        auto lowest = scores[0];
        auto index = 0;
        for (auto j = 1; j < scores.size(); j++) {
            if (scores[j] < lowest) {

                lowest = scores[j];
                index = j;
            }
        }
        pop[index] = random_modify(pop[index], rgen);
        return pop;
    }
};




int main(int argc, char** argv)
{
    //auto help = arg(argc, argv, "help", false, "help screen");
    std::random_device rd;
    std::mt19937 rgen(rd());
    for (auto i = 0; i < argc; i++) {
        std::cout << argv[i]<<std::endl;
    }
    if (argc < 3|| std::string(argv[1])=="help") {
        std::cout << "random/file | size/filename | method | params";
        std::cout << "genetic params: population | iterations | terminal condition | crossing method | mutation method";
    }
    
    if (std::string(argv[1]) == "random") {
        if (std::atoi(argv[2]) % 3 != 0) { return 0; }
        
        problem_t random_problem = generate_problem(std::atoi(argv[2]), rgen);
        Solution random_sol = random_solution(random_problem, rgen);
        switch (std::atoi(argv[3])) {
        case 1:
            
     std::cout << "Hillclimb:\n";
     for (auto i = 0; i < (random_problem.size()*random_problem.size()); i++) {

         auto new_solution = random_modify(random_sol,rgen);
         if (new_solution.score() < random_sol.score()) {
             random_sol = new_solution;
             random_sol.print();
         }
         
     }
            break;
        case 2:
            
      std::cout << "Better Hillclimb:\n";
     for (auto i = 0; i < 200; i++) {
         auto new_solution = better_modify(random_sol);
         if (new_solution.score() < random_sol.score()) {
             random_sol = new_solution;
             random_sol.print();
         }
     }

            break;
        case 3:
            std::cout << "Tabu Search\n";
            random_sol = tabu_search(random_sol, std::atoi(argv[4]));
            random_sol.print();
            break;

        case 4:
            std::cout << "Sim annealinng\n";
            switch (std::stoi(argv[4])) {
            case 1: {
                std::cout << "1000/k";
                random_sol = sim_annealing(random_sol, [](int k) {return 1000 / (k + 1); }, rgen);
                break;
            }
            case 2:
            {
                std::cout << "a^k";
                random_sol = sim_annealing(random_sol, [](int k) {return pow(5, k); }, rgen);
                break;
            }
            case 3: {
                std::cout << "1/log(k)";
                random_sol = sim_annealing(random_sol, [](int k) {return 1 / (log(k)); }, rgen);
                break;
            }
            }
            
            random_sol.print();
            break;

        case 5:
        {
            std::cout << "Genetic\n";
            std::cout << argc << std::endl;
            int population = std::atoi(argv[4]);
            int iterations = std::atoi(argv[5]);
            std::string terminal = std::string(argv[6]);
            std::string cross = std::string(argv[7]);
            std::string mutation = std::string(argv[8]);
            tsp_config conf(population, iterations, cross, mutation, terminal, random_problem, rgen);

            //tsp_config conf(30, 10, "randomCross", "worstMut", "maxPopul", random_problem, rgen);
            auto sol = generic_algorithm<Solution>(conf, rgen);
            sol.print();
            break;
        }
        default:
            std::cout << "Error!";
            return 0;
            break;
        }
       
    }

    if (std::string(argv[1]) == "file") {
        std::string filename = std::string(argv[2]);
        std::ifstream file(filename);
        problem_t numbers;
        while (!file.eof()) {
            std::string text;
            std::getline(file, text);
            numbers.push_back(std::stoi(text));
        }
        Solution file_sol(std::make_shared<problem_t>(numbers));
        if (file_sol.size() % 3 != 0) { return 0; }
        switch (std::atoi(argv[3])) {
        case 1:

            std::cout << "Hillclimb:\n";
            for (auto i = 0; i < 1000; i++) {

                auto new_solution = random_modify(file_sol, rgen);
                if (new_solution.score() < file_sol.score()) {
                    file_sol = new_solution;
                    file_sol.print();
                }
            }
            
            break;
        case 2:

            std::cout << "Better Hillclimb:\n";
            for (auto i = 0; i < 200; i++) {
                auto new_solution = better_modify(file_sol);
                if (new_solution.score() < file_sol.score()) {
                    file_sol = new_solution;
                    file_sol.print();
                }
            }

            break;
        case 3:
            std::cout << "Tabu Search\n";
            file_sol = tabu_search(file_sol, std::atoi(argv[4]));
            file_sol.print();
            break;

        case 4:
            std::cout << "Sim annealinng\n";
            file_sol = sim_annealing(file_sol, [](int k) {return 100 / (k + 1); }, rgen);
            file_sol.print();
            break;

        case 5: {
            
            break;
        }
        default:
            std::cout << "Error!";
            return 0;
            break;
        }
    }
    
    //problem_t tp_problem = { 20,23,25,30,49,45,27,30,30,40,22,19 };
    
    problem_t tp_problem = generate_problem(12, rgen);
    Solution tp_solution = random_solution(tp_problem, rgen);
    
    if (tp_problem.size() % 3 != 0) {
        //error1
    }
    if (tp_problem.size() < 6) {
        //error2
    }
    
    /* std::cout << "standard deviation" << tp_solution.goal();
     tp_solution = random_modify(tp_solution);
     tp_solution.print();
     std::cout << "standard eviation" << tp_solution.goal();*/

     /*
     std::cout << "\nmean:" << std::accumulate(tp_problem.begin(), tp_problem.end(), 0.0) * 3 / tp_problem.size();
     std::cout << "Hillclimb:\n";
     for (auto i = 0; i < 200; i++) {
         auto new_solution = random_modify(tp_solution);
         if (new_solution.score() < tp_solution.score()) {
             tp_solution = new_solution;
             tp_solution.print();
         }
     }*/

     /*
     std::cout << "Better Hillclimb:\n";
     for (auto i = 0; i < 200; i++) {
         auto new_solution = better_modify(tp_solution);
         if (new_solution.score() < tp_solution.score()) {
             tp_solution = new_solution;
             tp_solution.print();
         }
     }
     */
     /*
     std::cout << "Tabu Search\n";
     tp_solution = tabu_search(tp_solution,50);
     tp_solution.print();
     */
    /*
    tsp_config conf(30,10,"randomCross","worstMut","maxPopul",tp_problem,rgen);
    auto sol = generic_algorithm<Solution>(conf,rgen);
    sol.print();
    */
    /*
    std::cout << "Sim annealinng\n";
    tp_solution = sim_annealing(tp_solution, [](int k) {return 1000 / (k + 1); }, rgen);
    tp_solution.print();
    */
    /*
    std::cout << "Brute Force:\n";
    brute_force(tp_solution).print();
    */
}
//parametr hillclimb
//wybor temperatura schemat