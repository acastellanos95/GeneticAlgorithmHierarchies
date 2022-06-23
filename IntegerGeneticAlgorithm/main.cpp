#include <iostream>
#include <limits>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include "Utils.h"

int main(int argc, char *argv[]) {
  double pm, pc;
  size_t Gmax, populationSize;
  std::string filename;
  if (argc == 1) {
    std::cout << "Introduzca la probabilidad de mutación (0.01 <= pm <= 0.1): ";
    std::cin >> pm;
    std::cout << "Introduzca la probabilidad de cruza (pc aprox. 0.5): ";
    std::cin >> pc;
    std::cout << "Introduzca el número máximo de generaciones (>=100): ";
    std::cin >> Gmax;
    std::cout << "Introduzca el tamaño máximo de población (>=50): ";
    std::cin >> populationSize;
    std::cout << "Nombrar archivo con reporte: ";
    std::cin >> filename;
  } /*else if (argc == 7) {
    pm = std::strtod(argv[1], nullptr);
    pc = std::strtod(argv[2], nullptr);
    Gmax = std::strtoul(argv[3], nullptr, 0);
    populationSize = std::strtoul(argv[4], nullptr, 0);
    filename = std::string(argv[5]);
    Rseed = std::strtof(argv[6], nullptr);
  } */else {
    throw std::runtime_error("programa no admite entradas en línea de comandos");
  }

  /* Semilla hardcoded ya que la tarea indica que solo se pide las variables de arriba */
    Rseed = 0.1045;

  // Verificación de condiciones (petición de clase)
  if (pm > 0.1)
    std::cout << "Advertencia!!!! probabilidad de mutación mayor que la sugerida (0.01 <= pm <= 0.1)\n";
  if (pm > 1.0 || pm < 0.0)
    throw std::runtime_error("Probabilidad de mutación no puede ser < 0.0 ó > 1.0");
  if (pc != 0.5)
    std::cout << "Advertencia!!!! probabilidad de cruza menor que la sugerida (aprox. 0.5)\n";
  if (pm > 1.0 || pm < 0.0)
    throw std::runtime_error("Probabilidad de cruza no puede ser < 0.0 ó > 1.0");
  if (Gmax < 100)
    std::cout << "Advertencia!!!! número de generaciones menor que la sugerida (>=100)\n";
  if (populationSize < 50)
    std::cout << "Advertencia!!!! tamaño de población menor que la sugerida (>=50)\n";

  // Inicializa población y generador de números aleatorios
  std::ofstream file(filename);
  randomize();
  std::vector<individual> oldPopulation(populationSize);
  for (auto &individuo: oldPopulation) {
    // Inicializar cromosoma
    std::string chrom(numVariables * numDecVariable, '0');
    for (size_t alleleIndex = 0; alleleIndex < numVariables * numDecVariable; ++alleleIndex) {
      auto value = rnd(0, 9);
      chrom[alleleIndex] = std::to_string(value)[0];
    }
    individuo.chromosome = chrom;

    // Decodificar
    individuo.x = decode(individuo.chromosome);
  }
  std::vector<individual> newPopulation(populationSize);

  // Reporte Inicial

  file << std::fixed;
  file << std::setprecision(5);
  file << "probabilidad de mutación: " << pm << ", probabilidad de cruza: " << pc << ", número de generaciones: "
       << Gmax << ", tamaño de población: " << populationSize << ", semilla: " << Rseed << '\n';
  file << std::setprecision(3);

  // Generaciones y variables globales para generación
  std::pair<size_t, individual> bestIndividual;
  // Gráfica de convergencia y jeraruias
  std::vector<std::pair<size_t, double>> convergence;
  double maxFx = -std::numeric_limits<double>::max();
  const double MaxHierarchySelection = 1.1;
  const double MinHierarchySelection = 2.0 - MaxHierarchySelection;
  for (size_t generationIndex = 1; generationIndex <= Gmax; ++generationIndex) {
    file << "--------------------- Generación " << generationIndex << " ---------------------\n";

    // Cálculo de aptitud
    std::vector<std::pair<size_t, double>> fitness(populationSize, {0, 0.0});
    individual generationBestIndividual;

    // Valores de f(x) a fitness y buscar el máximo valor
    for (size_t popIndex = 0; popIndex < populationSize; ++popIndex) {
      fitness[popIndex].second = function(oldPopulation[popIndex].x);
      fitness[popIndex].first = popIndex;
    }

    auto maxIndividualFx = std::max_element(fitness.begin(), fitness.end(),
                                            [](std::pair<size_t, double> i, std::pair<size_t, double> j) {
                                              return i.second < j.second;
                                            });
    maxFx = maxIndividualFx->second > maxFx ? maxIndividualFx->second : maxFx;

    for (size_t popIndex = 0; popIndex < populationSize; ++popIndex) {
      fitness[popIndex].second = maxFx - fitness[popIndex].second + 0.1 * maxFx;
      oldPopulation[popIndex].fitness = fitness[popIndex].second;
    }

    double avgfitness = avgPairVector(fitness);
    size_t maxElementIndex = std::max_element(fitness.begin(), fitness.end(),
                                              [](std::pair<size_t, double> i, std::pair<size_t, double> j) {
                                                return i.second < j.second;
                                              }) - fitness.begin();
    double maxfitness = std::max_element(fitness.begin(), fitness.end(),
                                         [](std::pair<size_t, double> i, std::pair<size_t, double> j) {
                                           return i.second < j.second;
                                         })->second;
    double minfitness = std::min_element(fitness.begin(), fitness.end(),
                                         [](std::pair<size_t, double> i, std::pair<size_t, double> j) {
                                           return i.second < j.second;
                                         })->second;

    // Mejor individuo de la generación
    generationBestIndividual = oldPopulation[maxElementIndex];
    if (bestIndividual.second.chromosome.empty()) {
      bestIndividual.second = generationBestIndividual;
      bestIndividual.first = generationIndex;
    } else {
      // Actualizar mejor individuo global para comparar en caso de que maxFx haya cambiado
      bestIndividual.second.fitness = maxFx - function(bestIndividual.second.x) + 0.1 * maxFx;
      if (bestIndividual.second.fitness < generationBestIndividual.fitness)
        bestIndividual.second = generationBestIndividual;
      bestIndividual.first = generationIndex;
    }

    file << std::setprecision(5) << "Media de aptitud de población: " << avgfitness << '\n';
    file << "Aptitud máxima de población: " << maxfitness << '\n';
    file << "Aptitud mínima de población: " << minfitness << '\n';
    file << "Mejor individuo de la generación: " << generationBestIndividual.chromosome << ", ("
         << std::setprecision(3) << generationBestIndividual.x[0] << "," << generationBestIndividual.x[1] << "), fitness: "
         << std::setprecision(5) << generationBestIndividual.fitness << '\n';
    file << "Mejor individuo global: " << bestIndividual.second.chromosome << ", ("
         << std::setprecision(3) << bestIndividual.second.x[0] << "," << bestIndividual.second.x[1] << "), fitness: "
         << std::setprecision(5) << bestIndividual.second.fitness << ", generación: " << bestIndividual.first << '\n';
    // Gráfica de convergencia
    convergence.emplace_back(generationIndex, generationBestIndividual.fitness);

    // Jerarquias
    std::sort(fitness.begin(), fitness.end(),
              [](std::pair<size_t, double> i, std::pair<size_t, double> j) { return i.second < j.second; });
    for (size_t hierarchyIndex = 0; hierarchyIndex < fitness.size(); ++hierarchyIndex) {
      size_t popIndex = fitness[hierarchyIndex].first;

      // jerarquia(i,t) - 1 = hierarchyIndex + 1 - 1 = hierarchyIndex
      oldPopulation[popIndex].expectedValue = MinHierarchySelection +
          (MaxHierarchySelection - MinHierarchySelection) *
              ((round((double) hierarchyIndex)) /
                  (round((double) populationSize) - 1.0));
    }

    double sumExpected = 0.0;
    for (auto &individuo: oldPopulation) {
      sumExpected += individuo.expectedValue;
    }

    // Seleccionamos pares para cruza (universal estocástica)
    std::vector<size_t> preselected;
    std::vector<std::pair<size_t, size_t>> selected;
    std::vector<size_t> freqSelected(populationSize, 0);
    double r1 = randomperc();
    size_t i = 0;
    for (double sum = 0.0; i < populationSize; i++)
      for (sum += oldPopulation[i].expectedValue; sum > r1; r1 += 1.0)
        preselected.push_back(i);

    // construimos el arreglo de selección
    for (size_t selectedIndex = 0; selectedIndex < populationSize; selectedIndex += 2) {
      selected.emplace_back(preselected[selectedIndex], preselected[selectedIndex + 1]);
    }

    // Cruza de un punto
    uniformCrossover(selected, oldPopulation, newPopulation, pc, file);

    std::vector<individual> newPopulationBeforeMutation = newPopulation;

    // Mutación
    mutations(newPopulation, pm, file);

    // Elitismo
    elitism(newPopulation, generationBestIndividual);

    // Actualizar x en la nueva población
    for (auto &individuo: newPopulation) {
      // Decodificar
      individuo.x = decode(individuo.chromosome);
    }

    auto tmpPopulation = oldPopulation;
    oldPopulation = newPopulation;
    newPopulation = tmpPopulation;
  }

  // Reporte final
  file << std::setprecision(10) << "Final f(x): " << function(bestIndividual.second.x) << '\n';
  file << std::setprecision(3) << "$(" << bestIndividual.second.x[0] << "," << bestIndividual.second.x[1] << ")$ & $"
       << std::setprecision(10) << function(bestIndividual.second.x) << "$ & $" << bestIndividual.second.fitness
       << "$\\\\\n";
  // Punto 5 convergencia de la mediana
//  std::ofstream fileConvergence("convergence.dat");
//  for (auto &pair: convergence) {
//    fileConvergence << pair.first << ", " << std::setprecision(10) << std::to_string(pair.second) << '\n';
//  }
//  fileConvergence.close();
  file.close();
  return 0;
}
