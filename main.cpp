#include <bits/stdc++.h>
#define POPULATION_SIZE 1000
#define MAX_NUMBER 2147483647

using namespace std;

const int N = 1e3 + 2;

// Global variables declaration
map<int, int> taskFrequency;
set<int> jobSet;
vector<vector<pair<int, int>>> jobs; // sequence of operations for each job
vector<vector<pair<int, int>>> machineSchedules; // the schedule for each machine
int numberOfMachines;
map<string, double> dp;
long long maxJob ;
double averageTime ;

// This function reads the test input from a file
// input: file path
// output : the sequence of operations (vector of pairs) for each job (vector of vectors of pairs)
vector<vector<pair<int, int>>> readJobs(const string& filePath) {
    ifstream file(filePath);
    if (!file.is_open()) {
        cerr << "Failed to open file." << endl;
        throw runtime_error("File could not be opened.");
    }

    file >> numberOfMachines;
    vector<vector<pair<int, int>>> jobList;
    string line;
    getline(file, line);

    while (getline(file, line)) {
        getline(file, line);
        istringstream iss(line);
        vector<pair<int, int>> job;

        string task;
        while (getline(iss, task, ',')) {
            istringstream taskStream(task);
            int machineNumber, time;
            taskStream >> machineNumber >> time;
            job.emplace_back(machineNumber, time);
        }

        jobList.push_back(job);
    }

    file.close();
    return jobList;
}

void initializeMachines() {
    machineSchedules.clear();
    for (int i = 0; i < numberOfMachines; i++) {
        machineSchedules.push_back({});
    }
}

double calculateFitness(const vector<int>& chromosome) {
    int numberOfJobs = jobs.size();
    double average = 0;
    double averageFactor = 1.0 / numberOfJobs;
    double maxJobTime = 0;

    vector<int> currentTaskNumber(numberOfJobs, 0);
    vector<long long> finishTime(numberOfJobs, 0);
    vector<long long> machineFinishTime(numberOfMachines, 0);

    initializeMachines();

    for (int i = 0; i < chromosome.size(); i++) {
        int jobNumber = chromosome[i];
        auto task = jobs[jobNumber - 1][currentTaskNumber[jobNumber - 1]];
        currentTaskNumber[jobNumber - 1]++;

        int machineNumber = task.first;
        long long idleTime = max(0LL, finishTime[jobNumber - 1] - machineFinishTime[machineNumber - 1]);

        if (idleTime > 0) {
            machineSchedules[machineNumber - 1].push_back({-1, idleTime});
        }

        machineSchedules[machineNumber - 1].push_back({jobNumber, task.second});
        average -= averageFactor * finishTime[jobNumber - 1];
        finishTime[jobNumber - 1] = idleTime + task.second + machineFinishTime[machineNumber - 1];
        machineFinishTime[machineNumber - 1] = finishTime[jobNumber - 1];
        average += averageFactor * finishTime[jobNumber - 1];

        if (finishTime[jobNumber - 1] > maxJobTime) {
            maxJobTime = finishTime[jobNumber - 1];
        }
    }
    maxJob = maxJobTime ;
    averageTime = average ;
    return 1/maxJobTime ;
}

// This struct acts as a comparator to compare fitness, it is used for sorting the priority queue
struct ChromosomeComparator {
    bool operator()(const vector<int>& a, const vector<int>& b) {
        return calculateFitness(a) < calculateFitness(b);
    }
};

void createInitialPopulation(priority_queue<vector<int>, vector<vector<int>>, ChromosomeComparator>& population) {
    if (jobs.empty()) {
        cout << "No jobs to create initial population." << endl;
        return;
    }

    vector<int> initialPopulation;

    for (const auto& job : jobs) {
        int numberOfJobInstances = job.size();
        while (numberOfJobInstances-- > 0) {
            initialPopulation.push_back(&job - &jobs[0] + 1);
        }
    }

    for (int i = 0; i < POPULATION_SIZE; i++) {
        vector<int> chromosome = initialPopulation;
        random_shuffle(chromosome.begin(), chromosome.end());
        population.push(chromosome);
    }
}

void addRandomChromosomes(vector<vector<int>>& newGeneration, const vector<vector<int>>& generation, double percent) {
    vector<int> chromosome = generation[0];

    for (int i = 0; i < static_cast<int>(POPULATION_SIZE * percent); i++) {
        random_shuffle(chromosome.begin(), chromosome.end());
        newGeneration.push_back(chromosome);
    }
}

static random_device rd;
static mt19937 gen(rd());

int generateRandom(int maxNumber) {
    return (rand() % maxNumber) + 1;
}

int generateWeightedIndex(int maxIndex) {
    vector<double> weights;
    for (int i = 1; i <= maxIndex; ++i) {
        weights.push_back(1.0 / i);
    }

    discrete_distribution<int> distrib(weights.begin(), weights.end());
    return distrib(gen) + 1;
}

long double generateProbability() {
    uniform_real_distribution<> dis(0.0, 1.0);
    return dis(gen);
}

vector<int> validateChromosome(vector<int> currentState) {
    map<int, int> currentFrequency;
    for (const auto& task : currentState) {
        currentFrequency[task]++;
    }

    vector<int> neededTasks, overloadedTasks;
    map<int, vector<int>> indices;
    for (int i = 0; i < currentState.size(); i++) {
        indices[currentState[i]].push_back(i);
    }

    for (const auto& job : jobSet) {
        if (currentFrequency[job] < taskFrequency[job]) {
            neededTasks.push_back(job);
        } else if (currentFrequency[job] > taskFrequency[job]) {
            overloadedTasks.push_back(job);
        }
    }

    while (!neededTasks.empty() || !overloadedTasks.empty()) {
        int neededTaskIndex = generateRandom(neededTasks.size()) - 1;
        int overloadedTaskIndex = generateRandom(overloadedTasks.size()) - 1;
        int overloadedTaskPos = generateRandom(indices[overloadedTasks[overloadedTaskIndex]].size()) - 1;

        currentState[indices[overloadedTasks[overloadedTaskIndex]][overloadedTaskPos]] = neededTasks[neededTaskIndex];
        indices[overloadedTasks[overloadedTaskIndex]].erase(indices[overloadedTasks[overloadedTaskIndex]].begin() + overloadedTaskPos);
        currentFrequency[neededTasks[neededTaskIndex]]++;
        currentFrequency[overloadedTasks[overloadedTaskIndex]]--;

        if (currentFrequency[neededTasks[neededTaskIndex]] == taskFrequency[neededTasks[neededTaskIndex]]) {
            neededTasks.erase(neededTasks.begin() + neededTaskIndex);
        }
        if (currentFrequency[overloadedTasks[overloadedTaskIndex]] == taskFrequency[overloadedTasks[overloadedTaskIndex]]) {
            overloadedTasks.erase(overloadedTasks.begin() + overloadedTaskIndex);
        }
    }

    return currentState;
}

vector<int> combineChromosomes(const vector<int>& first, const vector<int>& second, int cutOff) {
    vector<int> result;
    for (int i = 0; i < first.size(); i++) {
        if (i < cutOff) {
            result.push_back(first[i]);
        } else {
            result.push_back(second[i]);
        }
    }
    return result;
}

vector<int> crossoverChromosomes(const vector<int>& first, const vector<int>& second) {
    vector<int> result;
    int maxFitness = max(calculateFitness(first), calculateFitness(second));

    if (maxFitness == calculateFitness(first)) {
        result = first;
    } else {
        result = second;
    }

    int randomIndex = generateRandom(first.size());
    vector<int> currentState = combineChromosomes(first, second, randomIndex);
    currentState = validateChromosome(currentState);

    if (calculateFitness(currentState) > maxFitness) {
        result = currentState;
    }

    return result;
}

bool shouldTake(int index, int iterationNumber, int populationSize) {
    int adjustedIndex = generateWeightedIndex(populationSize);
    return adjustedIndex <= iterationNumber;
}

vector<double> generateExponentialDistribution(int size, double lambda, bool fromBack) {
    vector<double> distribution(size);
    if (fromBack) {
        for (int i = 0; i < size; i++) {
            distribution[i] = exp(-lambda * i);
        }
    } else {
        for (int i = 0; i < size; i++) {
            distribution[i] = exp(lambda * i);
        }
    }

    double sum = accumulate(distribution.begin(), distribution.end(), 0.0);
    for (auto& value : distribution) {
        value /= sum;
    }

    return distribution;
}

template <typename T>
vector<T> selectElements(const vector<T>& previousGeneration, int count, bool fromBack) {
    vector<T> result;
    vector<double> distribution = generateExponentialDistribution(previousGeneration.size(), 0.5, fromBack);
    vector<bool> taken(previousGeneration.size(), false);

    default_random_engine generator(random_device{}());
    discrete_distribution<int> dist(distribution.begin(), distribution.end());

    for (int i = 0; i < count; ++i) {
        int index = dist(generator);
        while (taken[index]) {
            index = (index + 1) % previousGeneration.size();
        }

        result.push_back(previousGeneration[index]);
        taken[index] = true;
    }

    return result;
}

vector<vector<int>> chooseRandomly(const vector<vector<int>>& previousGeneration, double percent, bool fromBack) {
    return selectElements(previousGeneration, static_cast<int>(previousGeneration.size() * percent), fromBack);
}

void mutateChromosome(vector<int>& chromosome) {
    int times = generateRandom(max(1LL, static_cast<long long>(chromosome.size()) / 2 - 1LL)) - 1;
    for (int i = 0; i < times; i++) {
        int l = 0, r = 0;
        while (l == r || chromosome[l] == chromosome[r]) {
            l = generateRandom(chromosome.size()) - 1;
            r = generateRandom(chromosome.size()) - 1;
        }

        if (generateProbability() >= generateProbability()) {
            swap(chromosome[l], chromosome[r]);
        }
    }
}

void applyCrossover(vector<vector<int>>& newGeneration, const vector<vector<int>>& population, long double probability) {
    vector<vector<int>> tempGeneration;

    for (int i = 1; i < static_cast<int>(population.size() * 0.3); i++) {
        tempGeneration.push_back(crossoverChromosomes(population[i], population[i - 1]));
    }

    tempGeneration.push_back(crossoverChromosomes(population[0], population.back()));

    for (int i = 0; i < static_cast<int>(population.size() * 0.7); i++) {
        int l = 0, r = 0;
        while (l == r) {
            l = generateRandom(population.size()) - 1;
            r = generateRandom(population.size()) - 1;
        }

        tempGeneration.push_back(crossoverChromosomes(population[l], population[r]));
    }

    vector<vector<int>> selectedChromosomes = chooseRandomly(tempGeneration, probability, false);

    for (const auto& chromosome : selectedChromosomes) {
        newGeneration.push_back(chromosome);
    }
}

void applyMutation(vector<vector<int>>& newGeneration, const vector<vector<int>>& population, long double probability) {
    vector<vector<int>> selectedChromosomes = chooseRandomly(population, probability, true);

    for (auto& chromosome : selectedChromosomes) {
        mutateChromosome(chromosome);
        newGeneration.push_back(chromosome);
    }
}

void applyReproduction(vector<vector<int>>& newGeneration, const vector<vector<int>>& population, long double probability) {
    vector<vector<int>> selectedChromosomes = chooseRandomly(population, probability, false);

    for (int i = 0; i < static_cast<int>(0.1 * population.size()); i++) {
        newGeneration.push_back(population[i]);
    }
}

vector<vector<int>> copyPopulation(priority_queue<vector<int>, vector<vector<int>>, ChromosomeComparator>& population) {
    vector<vector<int>> result;
    priority_queue<vector<int>, vector<vector<int>>, ChromosomeComparator> temp;

    while (!population.empty()) {
        vector<int> chromosome = population.top();
        population.pop();
        temp.push(chromosome);
        result.push_back(chromosome);
    }

    while (!temp.empty()) {
        population.push(temp.top());
        temp.pop();
    }

    return result;
}

void processGeneration(priority_queue<vector<int>, vector<vector<int>>, ChromosomeComparator>& previousPopulation) {
    vector<vector<int>> newGeneration;
    vector<vector<int>> population = copyPopulation(previousPopulation);

    applyCrossover(newGeneration, population, 0.50);
    applyMutation(newGeneration, population, 0.20);
    addRandomChromosomes(newGeneration, population, 0.10);
    applyReproduction(newGeneration, population, 0.20);

    while (!previousPopulation.empty()) {
        previousPopulation.pop();
    }

    for (const auto& chromosome : newGeneration) {
        previousPopulation.push(chromosome);
    }
}

int main() {
    const string filePath = "jobs.txt";
    jobs = readJobs(filePath);

    priority_queue<vector<int>, vector<vector<int>>, ChromosomeComparator> population;

    for (int i = 0; i < jobs.size(); i++) {
        taskFrequency[i + 1] = jobs[i].size();
        jobSet.insert(i + 1);
    }

    createInitialPopulation(population);
    int iterations =  100 ;
    set<vector<int>> uniqueChromosomes;

    while (iterations--) {
        processGeneration(population);
        cout << calculateFitness(population.top()) << endl;
        cout << "this is max job time: " << maxJob << endl ;
        cout << "this is average job time:" << averageTime << endl ;
        vector<vector<int>> currentPopulation = copyPopulation(population);

        for (const auto& chromosome : currentPopulation) {
            uniqueChromosomes.insert(chromosome);
        }
    }

    cout << "best fitness score:" << calculateFitness(population.top()) << endl;

    for (int i = 0 ; i < machineSchedules.size() ; i ++ ){
        cout << "--------------------------M" << i + 1 << "--------------------------\n" ;
        int time = 0;
        cout << time << " -- ";
        for (int j = 0 ; j < machineSchedules[i].size() ; j ++ ){
            time += machineSchedules[i][j].second;
            // -1 represents the idle state
            if(machineSchedules[i][j].first != -1)
                cout << "J" << machineSchedules[i][j].first << " -- " << time << " -- " ;
            else
                cout << "idle" << " -- " << time << " -- " ;

            if(j % 10 == 0 && j != 0)
                cout << endl;
        }
        cout << endl ;
        cout << endl ;
    }

    return 0;
}
