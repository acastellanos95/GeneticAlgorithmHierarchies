//
// Created by andre on 6/9/22.
//

#ifndef BINARYGENETICALGORITHM__UTILS_H_
#define BINARYGENETICALGORITHM__UTILS_H_

#include <vector>
#include <bitset>
#include <cmath>
#include "Individual.h"
#include "lib/RandomGA.h"

/* Definiciones de problema */
const std::vector<std::vector<double>> a{{-32.0, -16.0, 0.0,   16.0,  32.0,  -32.0, -16.0, 0.0,   16.0,  32.0,  -32.0, -16.0, 0.0, 16.0, 32.0, -32.0, -16.0, 0.0,  16.0, 32.0, -32.0, -16.0, 0.0,  16.0, 32.0},
                                         {-32.0, -32.0, -32.0, -32.0, -32.0, -16.0, -16.0, -16.0, -16.0, -16.0, 0.0,   0.0,   0.0, 0.0,  0.0,  16.0,  16.0,  16.0, 16.0, 16.0, 32.0,  32.0,  32.0, 32.0, 32.0}};
const size_t numVariables = 2;
const size_t numBitsVariable = 17; /* log2([65.536-(-65.536)x10^3]) = 17 */
const double l_sup = 65.536;
const double l_inf = -65.536;

/// Función a minimizar F5 de De Jong
/// \param x
/// \return f(x)
double function(const std::vector<double> &x) {
    double fjSum = 0.0;
    for (size_t j = 0; j < 25; ++j) {
        double tmp = round((double) (j + 1)) + pow((x[0] - a[0][j]), 6.0) + pow((x[1] - a[1][j]), 6.0);
        fjSum += 1.0 / tmp;
    }
    double tmp = 1 / 500.0;
    tmp += fjSum;

    double result = 1 / tmp;
    return result;
}

/// Gray a binario de clase
/// \param grayGene: cadena representa gen en binario
/// \return gen en código de gray
std::string grayToBinary(const std::string &grayGene) {
    std::string binary(grayGene.size(), '0');

    binary[0] = grayGene[0];
    auto valor = grayGene[0];

    for (size_t i = 1; i < grayGene.length(); i++) {
        if (grayGene[i] == '1')
            valor = valor == '1' ? '0' : '1';
        binary[i] = valor;
    }

    return binary;
}

/// Decodificar cromosoma a x O(numVariables)
/// \param cromosoma: Cadena binaria que representa
/// \return Vector x
std::vector<double> decode(const std::string &cromosoma) {
    // Resultado
    std::vector<double> result(numVariables);
    // Vector de bitset de cada variable x_i
    std::vector<std::bitset<numBitsVariable>> variableGenes(numVariables);

    // Copiar de cadena al vector de bitset
    for (size_t i = 0; i < numVariables; ++i) {
        std::string geneString = grayToBinary(cromosoma.substr(i * numBitsVariable, numBitsVariable));
        variableGenes[i] = std::bitset<numBitsVariable>(geneString);
    }

    // Decodificación
    for (size_t i = 0; i < numVariables; ++i) {
        double x_i =
                l_inf + ((double) variableGenes[i].to_ullong()) * ((l_sup - l_inf) / (pow(2.0, numBitsVariable) - 1.0));
        result[i] = x_i;
    }

    return result;
}

/// Calcular el promedio de un vector de pares usando el segundo
/// \param x
/// \return avg(x)
double avgPairVector(const std::vector<std::pair<size_t, double>> &x) {
    double sum = 0.0;
    for (auto &x_i: x) {
        sum += x_i.second;
    }
    return sum / (double) x.size();
}

/// Cruza uniforme
/// \param selected: arreglo de par de índices individuos en la población seleccionados para cruza
/// \param oldPopulation: población de padres
/// \param newPopulation: población de hijos
/// \param pc: probabilidad de cruza
/// \param file: archivo
void uniformCrossover(const std::vector<std::pair<size_t, size_t>> &selected, const std::vector<individual> &oldPopulation,
                 std::vector<individual> &newPopulation, const double pc, std::ofstream &file) {
    size_t numberCrossover = 0;
    size_t newPopIndex = 0;
    for (auto &selectPair: selected) {
        newPopulation[newPopIndex].chromosome = oldPopulation[newPopIndex].chromosome;
        newPopulation[newPopIndex + 1].chromosome = oldPopulation[newPopIndex].chromosome;
        bool flag = false;
        for (size_t alleleIndex = 0; alleleIndex < numVariables * numBitsVariable; ++alleleIndex) {
            // si flip es verdadero entonces el alelo del primer hijo es del primer padre
            if (flip(pc)) {
                newPopulation[newPopIndex].chromosome[alleleIndex] = oldPopulation[selectPair.first].chromosome[alleleIndex];
                newPopulation[newPopIndex + 1].chromosome[alleleIndex] = oldPopulation[selectPair.second].chromosome[alleleIndex];
            } else {
                /* El hijo tiene un alelo de otro padre (es una cruza), si no entrara aquí nunca entonces los hijos pasan como los padres */
                flag = true;
                newPopulation[newPopIndex + 1].chromosome[alleleIndex] = oldPopulation[selectPair.first].chromosome[alleleIndex];
                newPopulation[newPopIndex].chromosome[alleleIndex] = oldPopulation[selectPair.second].chromosome[alleleIndex];
            }
        }
        if (flag)
            ++numberCrossover;
        newPopIndex += 2;
    }

    file << "Número de cruzas: " << numberCrossover << '\n';
}

/// Muta los alelos de cada cromosoma de un individuo
/// \param newPopulation: población de hijos a mutar
/// \param pm: probabilidad de mutación
/// \return lista de mutaciones
void mutations(std::vector<individual> &newPopulation, const double pm, std::ofstream &file) {
    size_t numberMutations = 0;
    for (auto &individuo: newPopulation) {
        size_t individualNumberMutations = 0;
        for (auto &allele: individuo.chromosome) {
            if (flip(pm)) {
                ++numberMutations;
                ++individualNumberMutations;
                allele = allele == '0' ? '1' : '0';
            }
        }
        individuo.numberMutations = individualNumberMutations;
    }

    file << "Número de mutaciones totales: " << numberMutations << '\n';
}

/// Tomar el más apto y ponerlo en el primer hijo
/// \param newPopulation: población de hijos
/// \param bestIndividual: individuo más apto de la población de padres
void elitism(std::vector<individual> &newPopulation, individual &bestIndividual) {
    newPopulation[0] = bestIndividual;
}

#endif //BINARYGENETICALGORITHM__UTILS_H_
