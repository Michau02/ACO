#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <fstream>
#include <sstream>
#include <random>
#include <algorithm>
#include <chrono>

using namespace std;

chrono::high_resolution_clock::time_point start;
chrono::duration<double, milli> timePassed;
chrono::duration<double, milli> timePassed2;

int n, m, generations; // n miast, m mrowek, ilosc/liczba (dalej nie wiem) generacji
double alpha, beta, ro, Q; // alfa, beta, wspolczynnik parowania, sta³a feromonu

vector<vector<double>> distances;
vector<vector<double>> pheromones;

void loadFromFile(string fileName) {
    ifstream file(fileName);

    if (file.is_open()) {
        file >> n;
        m = n;

        // Dynamiczna alokacja pamiêci dla macierzy distances i pheromones
        distances.resize(n, vector<double>(n, 0));
        pheromones.resize(n, vector<double>(n, 1.0));

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                file >> distances[i][j];
            }
        }

        file.close();
    }
    else {
        cout << "Something went wrong with file opening." << endl;
    }
}

int calculateTotalDistance(const vector<int>& tour) {
    int totalDistance = 0;
    for (int i = 0; i < n - 1; ++i) {
        totalDistance += distances[tour[i]][tour[i + 1]];
    }
    totalDistance += distances[tour[n - 1]][tour[0]]; // Powrót do pierwszego miasta
    return totalDistance;
}

double estimateInitialPheromones() {
    vector<bool> visited(n, false);
    vector<int> tour;
    int startingCity = rand() % n;
    tour.push_back(startingCity);
    visited[startingCity] = true;

    int currentCity = tour.back();
    double minDistance = numeric_limits<double>::infinity();
    int nextCity = -1;

    for (int i = 1; i < n; ++i) {
        currentCity = tour.back();
        minDistance = numeric_limits<double>::infinity();
        nextCity = -1;

        for (int j = 0; j < n; ++j) {
            if (!visited[j] && distances[currentCity][j] < minDistance) {
                minDistance = distances[currentCity][j];
                nextCity = j;
            }
        }

        tour.push_back(nextCity);
        visited[nextCity] = true;
    }

    return static_cast<double>(m) / calculateTotalDistance(tour);
}

void initializePheromones() {
    double initialPheromoneValue = estimateInitialPheromones();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            pheromones[i][j] = initialPheromoneValue;
        }
    }
}

void updatePheromones(const vector<vector<int>>& antTours) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            pheromones[i][j] *= (1.0 - ro);
        }
    }

    for (int k = 0; k < m; ++k) {
        double tourLength = calculateTotalDistance(antTours[k]);
        for (int i = 0; i < n - 1; ++i) {
            int from = antTours[k][i];
            int to = antTours[k][i + 1];
            pheromones[from][to] += Q / tourLength;
        }
    }
}

vector<int> antColony() {
    vector<vector<int>> antTours(m, vector<int>(n, -1));

    // Lista dostêpnych miast dla ka¿dej mrówki
    vector<int> availableCities(n);
    for (int i = 0; i < n; i++) {
        availableCities[i] = i;
    }
    std::random_device rd;
    std::mt19937 g(rd());  // U¿yj mt19937 do generowania liczb losowych
    std::shuffle(availableCities.begin(), availableCities.end(), g);

    double epsilon = 1e-9;

    for (int i = 0; i < generations; i++) {
        // Generowanie tras dla mrówek
        cout << "|";
        for (int k = 0; k < m; k++) {
            vector<bool> visited(n, false);
            antTours[k][0] = availableCities[k]; // Losowy startowy punkt
            visited[antTours[k][0]] = true;

            for (int step = 1; step < n; step++) {
                
                int currentCity = antTours[k][step - 1];
                double totalProbability = 0.0;

                for (int nextCity = 0; nextCity < n; nextCity++) {
                    if (!visited[nextCity]) {
                        double probability = pow(pheromones[currentCity][nextCity], alpha)
                            * pow((1.0 / (distances[currentCity][nextCity] + epsilon)), beta);
                        totalProbability += probability;
                    }
                }

                double randomValue = ((double)rand() / RAND_MAX) * totalProbability;
                double cumulativeProbability = 0.0;

                for (int nextCity = 0; nextCity < n; nextCity++) {
                    if (!visited[nextCity]) {
                        double probability = pow(pheromones[currentCity][nextCity], alpha)
                            * pow((1.0 / (distances[currentCity][nextCity] + epsilon)), beta);
                        cumulativeProbability += probability;

                        if (cumulativeProbability >= randomValue) {
                            antTours[k][step] = nextCity;
                            visited[nextCity] = true;
                            break;
                        }
                    }
                }
            }
        }
        
        // Aktualizacja feromonu
        updatePheromones(antTours);
        
    }

    // Wybór najlepszej trasy
    double minDistance = numeric_limits<double>::infinity();
    vector<int> bestTour;

    for (int k = 0; k < m; ++k) {
        double tourLength = calculateTotalDistance(antTours[k]);
        if (tourLength < minDistance) {
            minDistance = tourLength;
            bestTour = antTours[k];
        }
    }

    return bestTour;
}

void saveToCSV(const vector<vector<int>>& tours, const vector<double>& distances, int iterations, string fileName, double averageCost, double time, int BESTDISTANCE) {
    ofstream outputFile("output.csv", ios::app);  // Otwórz w trybie dopisywania

    outputFile << "Nazwa pliku: " << fileName << "; " << "Iloœæ testow: " << iterations << ";" << "Ilosc generacji: " << generations << endl << "Znalezione sciezki: " << "; " << "Koszt sciezki" << endl;

    for (size_t i = 0; i < tours.size(); ++i) {
        for (size_t j = 0; j < tours[i].size(); ++j) {
            outputFile << tours[i][j];
            if (j < tours[i].size() - 1) {
                outputFile << " ";
            }
        }

        // Oddziel trasê od odleg³oœci œrednikiem
        outputFile << "; " << distances[i] << endl;
    }

    outputFile << "Srednia dlugosc: " << ";" << averageCost << endl << "Najkrotsza znaleziona trasa : " << BESTDISTANCE << endl << "Sredni czas[ms] : " << "; " << time << endl << endl;

    // Zamknij plik
    outputFile.close();
    cout << "Wyniki dzialania programu zostaly zapisane w pliku outputFile.csv" << endl;
}

void loadParametersFromFile(string iniFileName, double& alpha, double& beta, double& ro, double& Q, int& generations, string& txtFileName, bool& excel, int &numOfTests) {
    ifstream file(iniFileName);

    if (file.is_open()) {
        string line;

        while (getline(file, line)) {
            istringstream iss(line);
            string key, value;
            if (getline(iss, key, '=')) {
                getline(iss, value);
                if (key == "alpha") {
                    alpha = stod(value);
                }
                else if (key == "beta") {
                    beta = stod(value);
                }
                else if (key == "ro") {
                    ro = stod(value);
                }
                else if (key == "Q") {
                    Q = stod(value);
                }
                else if (key == "generations") {
                    generations = stoi(value);
                }
                else if (key == "fileName") {
                    txtFileName = value;
                }
                else if (key == "saveToExcel") {
                    if (value == "1") excel = true;
                    else if (value == "0") excel = false;
                    else {
                        excel = false;
                        cout << endl << endl << "Zostala podana zla wartosc w polu saveToExcel. Prosze podac 0 (false - wypis do konsoli) lub 1 (zapisanie danych w pliku.csv)" << endl << endl;
                    }
                }
                else if (key == "numOfTests") {
                    numOfTests = stoi(value);
                }
            }
        }

        file.close();
    }
    else {
        cout << "Something went wrong with file opening." << endl;
    }
}

void showResult(bool &excel, vector<int> &BESTSOLUTION, int BESTDISTANCE, vector<double> & allDistances, int &numOfTests, vector<vector<int>> &allTours, string &fileName, double time) {
    if (!excel) {
        cout << "Algorytm wykonano " << numOfTests << " razy." << endl << "Najlepsza znaleziona trasa: ";
        for (int city : BESTSOLUTION) {
            cout << city << " ";
        }
        cout << BESTSOLUTION[0] << endl << "Koszt znalezionej trasy: " << BESTDISTANCE << endl << "Average time: " << time << " ms" << endl;
    }
    else {
        double average = 0;
        for (double dist : allDistances) {
            average += dist;
        }
        average /= numOfTests;
        saveToCSV(allTours, allDistances, numOfTests, fileName, average, time, BESTDISTANCE);
    }
}

int main() {
    srand(static_cast<unsigned int>(time(NULL)));

    string fileName; 
    bool excel;
    int numOfTests;
    loadParametersFromFile("parametry.ini", alpha, beta, ro, Q, generations, fileName, excel, numOfTests);
    loadFromFile(fileName);
    cout << fileName << endl;
    vector<vector<int>> allTours;
    vector<double> allDistances;

    vector<int> BESTSOLUTION, bestTour;
    int distance, BESTDISTANCE;

    start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < numOfTests; i++) {
        initializePheromones();
        bestTour = antColony();
        distance = calculateTotalDistance(bestTour);
        allTours.push_back(bestTour);
        allDistances.push_back(distance);

        if (i == 0) {
            BESTDISTANCE = distance;
            BESTSOLUTION = bestTour;
        }
        else {
            if (distance < BESTDISTANCE) {
                BESTDISTANCE = distance;
                BESTSOLUTION = bestTour;
            }
        }
    }
    timePassed = std::chrono::high_resolution_clock::now() - start;
    system("cls");
    showResult(excel, BESTSOLUTION, BESTDISTANCE, allDistances, numOfTests, allTours, fileName, (ceil((timePassed.count()) * 1000) / 1000) / numOfTests);
    
    return 0;
}
