//
// Created by andre on 6/9/22.
//

#ifndef BINARYGENETICALGORITHM__INDIVIDUAL_H_
#define BINARYGENETICALGORITHM__INDIVIDUAL_H_

#include <vector>
#include <string>

/**
 * Estructura para representar un individuo
 */
struct individual {
  std::string chromosome; /* Cromosoma */
  std::vector<double> x; /* Cromosoma decodificado */
  double fitness; /* Valor de aptitud */
  double expectedValue; /* Valor esperado */
  size_t crossoverPoint; /* Punto de cruza */
  bool hasCrossover; /* Punto de cruza */
  size_t numberMutations; /* Número de mutaciones en el cromosoma */
  std::vector<unsigned long> parents = std::vector<unsigned long>(2, 0); /* Padres de individuo */
};

#endif //BINARYGENETICALGORITHM__INDIVIDUAL_H_
