/*

Copyright (C) 2018  A. Niño

This file is part of NetworkCommunities.

NetworkCommunities is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

NetworkCommunities is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with NetworkCommunities.  If not, see <https://www.gnu.org/licenses/>.

*/


/*****************************************************************************
*
* Authors: SciCom research group-E.S.I. Universidad de Castilla-La Mancha
*          Paseo de la Universidad 4, 13004-Ciudad Real. SPAIN
*
* Release date: September, 2018
*
* Purpose: Definition (header) and implementation of class FuzzyCommunities
*
*****************************************************************************/


/**
 * @file
 * Definition (header) and implementation file for class FuzzyCommunities.
 * This class is part of the package networkCommunities.
*/
#ifndef FUZZYCOMMUNITIES_H
#define FUZZYCOMMUNITIES_H

#include <iostream>
#include <fstream>
#include <math.h>
#include <omp.h>
#include <chrono>
#include <random>
#include <sstream>
#include <ctype.h>
#include <iomanip>

#include "SparseArray.h"

using namespace std;

namespace NetworkCommunities {

  /**
   * \brief Defines a class for determining fuzzy communities in networks. Here,
   *        the first node is labeled as 0.

   * @author SciCom research group-E.S.I. Universidad de Castilla-La Mancha.
   *         Paseo de la Universidad 4, 13004-Ciudad Real. SPAIN
   * @date November, 2018
   */
  

  class FuzzyCommunities {
    private:
      // Declaring constants common to all objects of the class
      static constexpr double ZERO = 0.0;
      static constexpr double ONE  = 1.0;
      static constexpr double TWO  = 2.0;
      static constexpr double FOUR = 4.0;
      static constexpr double HALF = 0.5;

      // Declaring member variables
      SparseArray* A;      // Adjacency matrix
      int N;               // Number of vertices
      int c;               // Number of communities minus 1.
      int M;               // Number of edges in the network;
      int* degree;         // Degree array
      double** U;          // Membership array
      double** grad;       // Gradient array
      int nIterations;     // Number of iterations in the minimization
      int nF;              // Number of function evaluations in the line search
      int nThreads;        // Number of threads to use in OpemMP
                            

    public:      
       
      /**
       * Constructor method. Reading from a file with the network given as an 
       * edge list. In this case, the two first data in the file must be the
       * number of nodes and the number of edges. In each edge definition line 
       * the two first data are mandatory and represent the origin and 
       * destination nodes. The third datum (if it exists)
       * in each edge line is the edge weight. If no third datum exists in the 
       * edge lines a default weight of one is set.
       * 
       * IMPORTANT: A symmetric adjacency matrix is assumed.
       * 
       * @param fileName Name of the file containing the network
       * @param ci Number of communities
       */
       FuzzyCommunities(string fileName, int ci){
        ifstream file;
        this->c = ci - 1;  // Communities minus 1
        nThreads = 1;      // Default value
        
        // Reading edge list
        try {
          file.open(fileName);
          file >> N;  // Reading number of nodes
          file >> M;  // Reading number of edges
          readFuzzyCommunities(file, true); // Auxiliary reading method
          file.close();
          
        } catch (ifstream::failure f) {        
          cout << "File " + fileName + " not found";
          exit (EXIT_FAILURE);
        }
      }
       

      /**
       * Constructor method. Reading from a file with the network given as an 
       * edge list. In this case, the number of nodes and the number of edges
       * are given in the constructor. In each edge definition line 
       * the two first data are mandatory and represent the origin and 
       * destination nodes. The third datum (if it exists)
       * in each edge line is the edge weight. If no third datum exists in the 
       * edge lines a default weight of one is set.
       * 
       * IMPORTANT: A symmetric adjacency matrix is assumed.
       * 
       * @param fileName Name of the file containing the network
       * @param ci Number of communities
       * @param nNodes Number of nodes in the network
       * @param nEdges Number of edges in the network
       */
       FuzzyCommunities(string fileName, int ci, int nNodes, int nEdges){
        ifstream file;

        this->c = ci - 1;  // Communities minus 1
        nThreads = 1;      // Default value

        N = nNodes; // Setting number of nodes
        M = nEdges; // Setting number of edges

        // Reading edge list
        try {
          file.open(fileName);
          readFuzzyCommunities(file, false);   // Auxiliary reading method
          file.close();
          
        } catch (ifstream::failure f) {        
          cout << "File " + fileName + " not found";
          exit (EXIT_FAILURE);
        }
      }
       
       
       
      /**
       * Determines the fuzzified modularity of the network.
       * @param U  The current U matrix
       * @return Q, the fuzzified modularity
       */
      double modularity(){
        int i, j, ki;
        double aux, aux1, aux2;
        
        // Contribution from pairs of nodes
        aux1 = aux2 = ZERO;
        for(i = 0; i < N; i++){
          ki = degree[i];
          aux1 += ki * ki * sii(i);
          for(j = i+1; j < N; j++){
            aux2 += ki * degree[j] * sij(i, j);
          }
        }
        
        // Additional contribution from connected nodes (edges)
        aux = ZERO;
        for (j = 0; j < M; j++){
          aux += sij(A->r(j), A->c(j)) * A->v(j);
        }
        
        return (aux / M -(aux1 / TWO + aux2) / (TWO * M * M));
      }
      
      
      
      /**
       * Finding communities using conjugate gradient in a symmetric network.
       * Applies the conjugate gradient (CG) algorithm to calculate the fuzzy 
       * communities membership of each vertex in the network. The local 
       * line search is performed using a Newton-Raphson method.
       */
      void findCommunities(){
        bool continueIteration = true, stopIter;
        int n_max_iterations, i, j, nStop;
        double** g = new double* [N];
        double** h = new double* [N];
        double alpha, gNum, gDen, gamma,
               aux, maxUij, uLimit, fn, fOld, fLimit;
        
        
        
        // Allocating memory for the matrices
        for (i = 0; i < N; i++){
          g[i] = new double [c];
          h[i] = new double [c];
        }

        
        // General constants
        n_max_iterations =  300;  // Max number of iterations allowed
        uLimit = 0.001;           // Limit for variation of U matrix elements
        fLimit = 0.1;             // Limit for D function variation
        
        fOld = ZERO;              // Initial value for previous function
     
        
        // Initializing U
        initializeU();        
        

        nF = 0; // Number of function evaluations in line searches
        nIterations = 0;
        fn = function();

        // Computing gradients for a symmetric network
        gradient(); 
        
        
        // Initializing CG
        for (i = 0; i < N; i++){
          for (j =0; j < c; j++) {
            g[i][j] = h[i][j] = -grad[i][j]; // Search direction
          }
        }
        

        stopIter = false; 
        nStop = 0;
        while(continueIteration){      

          // Computing minimum along the search direction h
          alpha = alpha_NR(h);
         
          
          // Updating coordinates
          maxUij = ZERO;
          for(i = 0; i < N; i++){
            for(j = 0; j < c; j++){
              aux = U[i][j];
              U[i][j] +=  alpha * h[i][j];
              aux = fabs(U[i][j] - aux);
              if (aux > maxUij) maxUij = aux;
            }
          }

          if (maxUij < uLimit){
            nStop++;
            if (nStop == 2) stopIter = true;
          } else{
            nStop = 0;
          }


          if(stopIter || (fabs(fOld - fn) < fLimit) || 
            nIterations >= n_max_iterations){   // End conditions

             continueIteration = false;
          } else {            
            fOld = fn;
            fn = function();
            
                        
            // Computing gradients
            gradient(); 
  
            
            // Computing gamma_i (Polak-Ribière)
            gNum = gDen = ZERO;
            for (i = 0; i < N; i++){
              for (j = 0; j < c; j++) {
                gNum -= (-grad[i][j] - g[i][j]) * grad[i][j];
                gDen += g[i][j] * g[i][j];
              }
            }
            gamma = gNum / gDen; 
            
            
            for (i = 0; i < N; i++){
              for (j = 0; j < c; j++) {
                g[i][j] = -grad[i][j];                // Updating g
                h[i][j] = g[i][j] + gamma * h[i][j];  // Updating h
              }
            }
            
            // Increasing counter
            nIterations++; 
          }
        }
          
        
        // Deallocating auxiliary arrays
        for (i = 0; i < N; i++){
          delete [] g[i];
          delete [] h[i];
        }
        delete [] g;
        delete [] h;
      }    
 
 
      
      /**
       * Finding communities using conjugate gradient in a symmetric network.
       * Applies the conjugate gradient (CG) algorithm to calculate the fuzzy 
       * communities membership of each vertex in the network. The local 
       * line search is performed using the Brent's method.
       */
      void findCommunities_Brent(){
        bool continueIteration = true, stopIter;
        int n_max_iterations, i, j, nStop;
        double** g = new double* [N];
        double** h = new double* [N];
        double alpha, gNum, gDen, gamma, 
               aux, maxUij, uLimit, fn, fOld, fLimit;
        
        
        
        // Allocating memory for the matrices
        for (i = 0; i < N; i++){
          g[i] = new double [c];
          h[i] = new double [c];
        }

        
        // General constants
        n_max_iterations =  300;  // Max number of iterations allowed
        uLimit = 0.001;           // Limit for variation of U matrix elements
        fLimit = 0.1;             // Limit for D function variation
        
        fOld = ZERO;              // Initial value for previous function
     
        
        // Initializing U
        initializeU();
        

        nF = 0; // Number of function evaluations in line searches
        nIterations = 0;
        fn = function();

        // Computing gradients for a symmetric network
        gradient(); 
        
        
        // Initializing CG
        for (i = 0; i < N; i++){
          for (j =0; j < c; j++) {
            g[i][j] = h[i][j] = -grad[i][j]; // Search direction
          }
        }
        

        stopIter = false; 
        nStop = 0;
        while(continueIteration){      

          // Computing minimum along the search direction h
          alpha = alpha_Brent(h);
          
          
          // Updating coordinates
          maxUij = ZERO;
          for(i = 0; i < N; i++){
            for(j = 0; j < c; j++){
              aux = U[i][j];
              U[i][j] +=  alpha * h[i][j];
              aux = fabs(U[i][j] - aux);
              if (aux > maxUij) maxUij = aux;
            }
          }
  

          if (maxUij < uLimit){
            nStop++;
            if (nStop == 2) stopIter = true;
          } else{
            nStop = 0;
          }


          if(stopIter || (fabs(fOld - fn) < fLimit) || 
            nIterations >= n_max_iterations){   // End conditions

             continueIteration = false;
          } else {            
            fOld = fn;
            fn = function();
            
                        
            // Computing gradients
            gradient(); 
  
            
            // Computing gamma_i (Polak-Ribière)
            gNum = gDen = ZERO;
            for (i = 0; i < N; i++){
              for (j = 0; j < c; j++) {
                gNum -= (-grad[i][j] - g[i][j]) * grad[i][j];
                gDen += g[i][j] * g[i][j];
              }
            }
            gamma = gNum / gDen; 
            
            
            
            for (i = 0; i < N; i++){
              for (j = 0; j < c; j++) {
                g[i][j] = -grad[i][j];                // Updating g
                h[i][j] = g[i][j] + gamma * h[i][j];  // Updating h
              }
            }

            // Increasing counter
            nIterations++; 
          }
        }
          
        
        // Deallocating auxiliary arrays
        for (i = 0; i < N; i++){
          delete [] g[i];
          delete [] h[i];
        }
        delete [] g;
        delete [] h;
      
      }    
 
      
      /**
       * Finding communities using steepest descent in a symmetric network.
       * Applies the steepest descent (SD) algorithm to calculate the fuzzy 
       * communities membership of each vertex in the network. The local 
       * line search is performed using a Newton-Raphson method.
       */
      void findCommunities_SD(){
        bool continueIteration = true, stopIter;
        int n_max_iterations, i, j, nStop;
        double** h = new double* [N];
        double alpha, aux, maxUij, uLimit, fn, fOld, fLimit;        
        
        
        // Allocating memory for the h matrix
        for (i = 0; i < N; i++){
          h[i] = new double [c];
        }

        
        // General constants
        n_max_iterations =  300;  // Max number of iterations allowed
        uLimit = 0.001;           // Limit for variation of U matrix elements
        fLimit = 0.1;             // Limit for D function variation
        
        fOld = ZERO;              // Initial value for previous function
     
        
        // Initializing U
        initializeU();        
        

        nF = 0; // Number of function evaluations in line searches
        nIterations = 0;
        fn = function();

        // Computing gradients for a symmetric network
        gradient(); 
        
        
        // Initializing SD
        for (i = 0; i < N; i++){
          for (j =0; j < c; j++) {
            h[i][j] = -grad[i][j]; // Search direction
          }
        }
        

        stopIter = false; 
        nStop = 0;
        while(continueIteration){      

          // Computing minimum along the search direction h
          alpha = alpha_NR(h);
         
          
          // Updating coordinates
          maxUij = ZERO;
          for(i = 0; i < N; i++){
            for(j = 0; j < c; j++){
              aux = U[i][j];
              U[i][j] +=  alpha * h[i][j];
              aux = fabs(U[i][j] - aux);
              if (aux > maxUij) maxUij = aux;
            }
          }

          if (maxUij < uLimit){
            nStop++;
            if (nStop == 2) stopIter = true;
          } else{
            nStop = 0;
          }


          if(stopIter || (fabs(fOld - fn) < fLimit) || 
            nIterations >= n_max_iterations){   // End conditions

             continueIteration = false;
          } else {            
            fOld = fn;
            fn = function();
            
                        
            // Computing gradients
            gradient(); 
                        
            for (i = 0; i < N; i++){
              for (j = 0; j < c; j++) {
                h[i][j] = -grad[i][j];    // Updating h
              }
            }
            
            // Increasing counter
            nIterations++; 
          }
        }
          
        
        // Deallocating auxiliary array
        for (i = 0; i < N; i++){
          delete [] h[i];
        }
        delete [] h;
      }    
 
 
      
      /**
       * Finding communities using steepest descent in a symmetric network.
       * Applies the steepest descent (SD) algorithm to calculate the fuzzy 
       * communities membership of each vertex in the network. The local 
       * line search is performed using the Brent's method.
       */
      void findCommunities_SD_Brent(){
        bool continueIteration = true, stopIter;
        int n_max_iterations, i, j, nStop;
        double** h = new double* [N];
        double alpha, aux, maxUij, uLimit, fn, fOld, fLimit;        
        
        
        // Allocating memory for the h matrix
        for (i = 0; i < N; i++){
          h[i] = new double [c];
        }

        
        // General constants
        n_max_iterations =  300;  // Max number of iterations allowed
        uLimit = 0.001;           // Limit for variation of U matrix elements
        fLimit = 0.1;             // Limit for D function variation
        
        fOld = ZERO;              // Initial value for previous function
     
        
        // Initializing U
        initializeU();
        

        nF = 0; // Number of function evaluations in line searches
        nIterations = 0;
        fn = function();

        // Computing gradients for a symmetric network
        gradient(); 
        
        
        // Initializing SD
        for (i = 0; i < N; i++){
          for (j =0; j < c; j++) {
            h[i][j] = -grad[i][j]; // Search direction
          }
        }
        

        stopIter = false; 
        nStop = 0;
        while(continueIteration){      

          // Computing minimum along the search direction h
          alpha = alpha_Brent(h);
          
          
          // Updating coordinates
          maxUij = ZERO;
          for(i = 0; i < N; i++){
            for(j = 0; j < c; j++){
              aux = U[i][j];
              U[i][j] +=  alpha * h[i][j];
              aux = fabs(U[i][j] - aux);
              if (aux > maxUij) maxUij = aux;
            }
          }
  

          if (maxUij < uLimit){
            nStop++;
            if (nStop == 2) stopIter = true;
          } else{
            nStop = 0;
          }


          if(stopIter || (fabs(fOld - fn) < fLimit) || 
            nIterations >= n_max_iterations){   // End conditions

             continueIteration = false;
          } else {            
            fOld = fn;
            fn = function();
            
                        
            // Computing gradients
            gradient(); 
            
            for (i = 0; i < N; i++){
              for (j = 0; j < c; j++) {
                h[i][j] = -grad[i][j];    // Updating h
              }
            }

            // Increasing counter
            nIterations++; 
          }
        }
          
        
        // Deallocating auxiliary array
        for (i = 0; i < N; i++){
          delete [] h[i];
        }
        delete [] h;
      
      }   
      
      
      
      /**
       * Finding communities in a BA-benchmark network using 
       * conjugate gradient in a symmetric network. The method 
       * applies the conjugate gradient (CG) algorithm to calculate the fuzzy 
       * communities membership of each vertex in the network. The local 
       * line search is performed using the Newton-Raphson method. The method
       * differs from the other find communities methods in that it applies a 
       * taylored intializeU_Benchmark() method.
       */
      void findCommunities_benchmark(){
        bool continueIteration = true, stopIter;
        int n_max_iterations, i, j, nStop;
        double** g = new double* [N];
        double** h = new double* [N];
        double alpha, gNum, gDen, gamma, 
               aux, maxUij, uLimit, fn, fOld, fLimit;
        
        
        // Allocating memory for the matrices
        for (i = 0; i < N; i++){
          g[i] = new double [c];
          h[i] = new double [c];
        }

        
        // General constants
        n_max_iterations =  300;  // Max number of iterations allowed
        uLimit = 0.001;           // Limit for variation of U matrix elements
        fLimit = 0.1;             // Limit for D function variation
        
        fOld = ZERO;              // Initial value for previous function
     
        
        // Initializing U
        initializeU_benchmark();        
        

        nF = 0; // Number of function evaluations in line searches
        nIterations = 0;
        fn = function();

        // Computing gradients for a symmetric network
        gradient(); 
        
        
        // Initializing CG
        for (i = 0; i < N; i++){
          for (j =0; j < c; j++) {
            g[i][j] = h[i][j] = -grad[i][j]; // Search direction
          }
        }
        

        stopIter = false; 
        nStop = 0;
        while(continueIteration){      

          // Computing minimum along the search direction h
          alpha = alpha_NR(h);
          
          
          // Updating coordinates
          maxUij = ZERO;
          for(i = 0; i < N; i++){
            for(j = 0; j < c; j++){
              aux = U[i][j];
              U[i][j] +=  alpha * h[i][j];
              aux = fabs(U[i][j] - aux);
              if (aux > maxUij) maxUij = aux;
            }
          }
  

          if (maxUij < uLimit){
            nStop++;
            if (nStop == 2) stopIter = true;
          } else{
            nStop = 0;
          }


          if(stopIter || (fabs(fOld - fn) < fLimit) || 
            nIterations >= n_max_iterations){   // End conditions

             continueIteration = false;
          } else {            
            fOld = fn;
            fn = function();
            
                        
            // Computing gradients
            gradient(); 
  
            
            // Computing gamma_i (Polak-Ribière)
            gNum = gDen = ZERO;
            for (i = 0; i < N; i++){
              for (j = 0; j < c; j++) {
                gNum -= (-grad[i][j] - g[i][j]) * grad[i][j];
                gDen += g[i][j] * g[i][j];
              }
            }
            gamma = gNum / gDen; 
            
            
            for (i = 0; i < N; i++){
              for (j = 0; j < c; j++) {
                g[i][j] = -grad[i][j];                // Updating g
                h[i][j] = g[i][j] + gamma * h[i][j];  // Updating h
              }
            }
            
            // Increasing counter
            nIterations++; 
          }
        }
          
        
        // Deallocating auxiliary arrays
        for (i = 0; i < N; i++){
          delete [] g[i];
          delete [] h[i];
        }
        delete [] g;
        delete [] h;
      
      }    
 

      
      /**
       * Finding communities in a BA-benchmark network using 
       * conjugate gradient in a symmetric network. The method 
       * applies the conjugate gradient (CG) algorithm to calculate the fuzzy 
       * communities membership of each vertex in the network. The local 
       * line search is performed using the Brent's method.
       */
      void findCommunities_benchmark_Brent(){
        bool continueIteration = true, stopIter;
        int n_max_iterations, i, j, nStop;
        double** g = new double* [N];
        double** h = new double* [N];
        double alpha, gNum, gDen, gamma, 
               aux, maxUij, uLimit, fn, fOld, fLimit;
        
        
        
        // Allocating memory for the matrices
        for (i = 0; i < N; i++){
          g[i] = new double [c];
          h[i] = new double [c];
        }

        
        // General constants
        n_max_iterations =  300;  // Max number of iterations allowed
        uLimit = 0.001;           // Limit for variation of U matrix elements
        fLimit = 0.1;             // Limit for D function variation
        
        fOld = ZERO;              // Initial value for previous function
     
        
        // Initializing U
        initializeU_benchmark();        
        

        nF = 0; // Number of function evaluations in line searches
        nIterations = 0;
        fn = function();

        
        // Computing gradients for a similarity based network       
        gradient(); 
        
        
        // Initializing CG
        for (i = 0; i < N; i++){
          for (j =0; j < c; j++) {
            g[i][j] = h[i][j] = -grad[i][j]; // Search direction
            cout << grad[i][j]<<endl;
          }
        }
        

        stopIter = false; 
        nStop = 0;
        while(continueIteration){      

          // Computing minimum along the search direction h
          alpha = alpha_Brent(h);
          
          
          // Updating coordinates
          maxUij = ZERO;
          for(i = 0; i < N; i++){
            for(j = 0; j < c; j++){
              aux = U[i][j];
              U[i][j] +=  alpha * h[i][j];
              aux = fabs(U[i][j] - aux);
              if (aux > maxUij) maxUij = aux;
            }
          }
  

          if (maxUij < uLimit){
            nStop++;
            if (nStop == 2) stopIter = true;
          } else{
            nStop = 0;
          }


          if(stopIter || (fabs(fOld - fn) < fLimit) || 
            nIterations >= n_max_iterations){   // End conditions

             continueIteration = false;
          } else {            
            fOld = fn;
            fn = function();
            
                        
            // Computing gradients
            gradient(); 
  
            
            // Computing gamma_i (Polak-Ribière)
            gNum = gDen = ZERO;
            for (i = 0; i < N; i++){
              for (j = 0; j < c; j++) {
                gNum -= (-grad[i][j] - g[i][j]) * grad[i][j];
                gDen += g[i][j] * g[i][j];
              }
            }
            gamma = gNum / gDen; 
            
            
            for (i = 0; i < N; i++){
              for (j = 0; j < c; j++) {
                g[i][j] = -grad[i][j];                // Updating g
                h[i][j] = g[i][j] + gamma * h[i][j];  // Updating h
              }
            }
            
            // Increasing counter
            nIterations++; 
          }
        }
          
        
        // Deallocating auxiliary arrays
        for (i = 0; i < N; i++){
          delete [] g[i];
          delete [] h[i];
        }
        delete [] g;
        delete [] h;
      
      }    
 
      
       /**
       * Finding communities in parallel using conjugate gradient in a symmetric
       * network. The method applies the conjugate gradient (CG) algorithm 
       * to calculate the fuzzy communities membership of each vertex in the
       *  network. The local line search is performed using the 
       * Newton-Raphson method.
       *  
       * @param n Number of threads to use in shared memory (OpenMP) 
       *        parallelism
       */
      void findCommunities_parallel(int n){
        bool continueIteration = true, stopIter;
        int n_max_iterations, i, j;
        double** g = new double* [N];
        double** h = new double* [N];
        double alpha, gNum, gDen, gamma, 
               aux, maxUij, uLimit, fn, fOld, fLimit;
        
        nThreads = n;  // Number of threads a usar
        
        // Allocating memory for the matrices
        #pragma omp parallel for num_threads (nThreads)
        for (i = 0; i < N; i++){
          g[i] = new double [c];
          h[i] = new double [c];
        }

        // General constants
        n_max_iterations =  300;  // Max number of iterations allowed
        uLimit = 0.001;           // Limit for variation of U matrix elements
        fLimit = 0.1;             // Limit for D function variation
        
        fOld = ZERO;              // Initial value for previous function

        // Initializing U
        initializeU_parallel();

        nF = 0; // Number of function evaluations in line searches
        nIterations = 0;
        fn = function_parallel();

        
        // Initializing Conjugate Gradient
        
        #pragma omp parallel for num_threads (nThreads)
        for (i = 0; i < N; i++){
          for (j =0; j < c; j++) {
            aux = dulk(i, j); // Computing gradient on the fly
            g[i][j] = h[i][j] = -aux; // Search direction
          }
        }
        

        stopIter = false; 
        while(continueIteration){      

          // Computing minimum along the search direction h
          alpha = alpha_NR_parallel(h);
          

          // Updating coordinates
          maxUij = ZERO;
          
          #pragma omp parallel for num_threads (nThreads) \
                  private (j, aux) reduction (max: maxUij)
          for(i = 0; i < N; i++){
            for(j = 0; j < c; j++){
              aux = U[i][j];
              U[i][j] +=  alpha * h[i][j];
              aux = fabs(U[i][j] - aux);
              if (aux > maxUij) maxUij = aux;
            }
          }
  
          if (maxUij < uLimit) stopIter = true;

          if(stopIter || (fabs(fOld - fn) < fLimit) || 
            nIterations >= n_max_iterations){   // End conditions
             continueIteration = false;
          } else {      
            fOld = fn;
            fn = function_parallel();


            // Computing gamma_i (Polak-Ribière)
            gNum = gDen = ZERO;
            
            #pragma omp parallel num_threads (nThreads)
            {
              #pragma omp for private (j) reduction(+: gDen) reduction(-: gNum)
              for (i = 0; i < N; i++){
                for (j = 0; j < c; j++) {
                  grad[i][j] = dulk(i, j); // Computing gradient on the fly
                  gNum -= (-grad[i][j] - g[i][j]) * grad[i][j];
                  gDen += g[i][j] * g[i][j];
                }
              }

              #pragma omp single
              gamma = gNum / gDen; 

              #pragma omp for private (j)
              for (i = 0; i < N; i++){
                for (j = 0; j < c; j++) {
                  g[i][j] = -grad[i][j];                // Updating g
                  h[i][j] = g[i][j] + gamma * h[i][j];  // Updating h
                }
              }
          }
            // Increasing counter
            nIterations++; 
          }
        }
        
        // Deallocating auxiliary arrays
        #pragma omp parallel for
        for (i = 0; i < N; i++){
          delete [] g[i];
          delete [] h[i];
        }
        delete [] g;
        delete [] h;
      
      }    

      
       /**
       * Finding communities in parallel using conjugate gradient in a symmetric
       * network. The method applies the conjugate gradient (CG) algorithm 
       * to calculate the fuzzy communities membership of each vertex in the
       * network. The local line search is performed using the Brent's method.
       *  
       * @param n Number of threads to use in shared memory (OpenMP) 
       *        parallelism
       */
      void findCommunities_Brent_parallel(int n){
        bool continueIteration = true, stopIter;
        int n_max_iterations, i, j;
        double** g = new double* [N];
        double** h = new double* [N];
        double alpha, gNum, gDen, gamma, 
               aux, maxUij, uLimit, fn, fOld, fLimit;
        
        nThreads = n;  // Number of threads a usar
        
        // Allocating memory for the matrices
        #pragma omp parallel for num_threads (nThreads)
        for (i = 0; i < N; i++){
          g[i] = new double [c];
          h[i] = new double [c];
        }

        // General constants
        n_max_iterations =  300;  // Max number of iterations allowed
        uLimit = 0.001;           // Limit for variation of U matrix elements
        fLimit = 0.1;             // Limit for D function variation
        
        fOld = ZERO;              // Initial value for previous function

        // Initializing U
        initializeU_parallel();

        nF = 0; // Number of function evaluations in line searches
        nIterations = 0;
        fn = function_parallel();

        
        // Initializing Conjugate Gradient
        
        #pragma omp parallel for num_threads (nThreads)
        for (i = 0; i < N; i++){
          for (j =0; j < c; j++) {
            aux = dulk(i, j); // Computing gradient on the fly
            g[i][j] = h[i][j] = -aux; // Search direction
          }
        }
        

        stopIter = false; 
        while(continueIteration){      

          // Computing minimum along the search direction h
          alpha = alpha_Brent_parallel(h);
          

          // Updating coordinates
          maxUij = ZERO;
          
          #pragma omp parallel for num_threads (nThreads) \
                  private (j, aux) reduction (max: maxUij)
          for(i = 0; i < N; i++){
            for(j = 0; j < c; j++){
              aux = U[i][j];
              U[i][j] +=  alpha * h[i][j];
              aux = fabs(U[i][j] - aux);
              if (aux > maxUij) maxUij = aux;
            }
          }
  
          if (maxUij < uLimit) stopIter = true;

          if(stopIter || (fabs(fOld - fn) < fLimit) || 
            nIterations >= n_max_iterations){   // End conditions
             continueIteration = false;
          } else {      
            fOld = fn;
            fn = function_parallel();


            // Computing gamma_i (Polak-Ribière)
            gNum = gDen = ZERO;
            
            #pragma omp parallel num_threads (nThreads)
            {
              #pragma omp for private (j) reduction(+: gDen) reduction(-: gNum)
              for (i = 0; i < N; i++){
                for (j = 0; j < c; j++) {
                  grad[i][j] = dulk(i, j); // Computing gradient on the fly
                  gNum -= (-grad[i][j] - g[i][j]) * grad[i][j];
                  gDen += g[i][j] * g[i][j];
                }
              }

              #pragma omp single
              gamma = gNum / gDen; 

              #pragma omp for private (j)
              for (i = 0; i < N; i++){
                for (j = 0; j < c; j++) {
                  g[i][j] = -grad[i][j];                // Updating g
                  h[i][j] = g[i][j] + gamma * h[i][j];  // Updating h
                }
              }
          }
            // Increasing counter
            nIterations++; 
          }
        }
        
        // Deallocating auxiliary arrays
        #pragma omp parallel for
        for (i = 0; i < N; i++){
          delete [] g[i];
          delete [] h[i];
        }
        delete [] g;
        delete [] h;
      
      }    

       
       /**
       * Finding communities in parallel using steepest descent in a symmetric
       * network. The method applies the steepest descent (SD) algorithm 
       * to calculate the fuzzy communities membership of each vertex in the
       * network. The local line search is performed using the 
       * Newton-Raphson method.
       *  
       * @param n Number of threads to use in shared memory (OpenMP) 
       *        parallelism
       */
      void findCommunities_SD_parallel(int n){
        bool continueIteration = true, stopIter;
        int n_max_iterations, i, j;
        double** h = new double* [N];
        double alpha, aux, maxUij, uLimit, fn, fOld, fLimit;
        
        nThreads = n;  // Number of threads a usar
        
        // Allocating memory for the matrices
        #pragma omp parallel for num_threads (nThreads)
        for (i = 0; i < N; i++){
          h[i] = new double [c];
        }

        // General constants
        n_max_iterations =  300;  // Max number of iterations allowed
        uLimit = 0.001;           // Limit for variation of U matrix elements
        fLimit = 0.1;             // Limit for D function variation
        
        fOld = ZERO;              // Initial value for previous function

        // Initializing U
        initializeU_parallel();

        nF = 0; // Number of function evaluations in line searches
        nIterations = 0;
        fn = function_parallel();

        
        // Initializing steepest descent
        
        #pragma omp parallel for num_threads (nThreads)
        for (i = 0; i < N; i++){
          for (j =0; j < c; j++) {
            h[i][j] = -dulk(i, j); // Computing gradient on the fly
          }
        }
        

        stopIter = false; 
        while(continueIteration){      

          // Computing minimum along the search direction h
          alpha = alpha_NR_parallel(h);
          

          // Updating coordinates
          maxUij = ZERO;
          
          #pragma omp parallel for num_threads (nThreads) \
                  private (j, aux) reduction (max: maxUij)
          for(i = 0; i < N; i++){
            for(j = 0; j < c; j++){
              aux = U[i][j];
              U[i][j] +=  alpha * h[i][j];
              aux = fabs(U[i][j] - aux);
              if (aux > maxUij) maxUij = aux;
            }
          }
  
          if (maxUij < uLimit) stopIter = true;

          if(stopIter || (fabs(fOld - fn) < fLimit) || 
            nIterations >= n_max_iterations){   // End conditions
             continueIteration = false;
          } else {      
            fOld = fn;
            fn = function_parallel();

            
            #pragma omp parallel for num_threads (nThreads) private (j)
              for (i = 0; i < N; i++){
                for (j = 0; j < c; j++) {
                  h[i][j] = -dulk(i, j); // Computing gradient on the fly
                }
              }


            // Increasing counter
            nIterations++; 
          }
        }
        
        // Deallocating auxiliary array
        #pragma omp parallel for
        for (i = 0; i < N; i++){
          delete [] h[i];
        }
        delete [] h;
      
      }    

      
       /**
       * Finding communities in parallel using steepest descent in a symmetric
       * network. The method applies the steepest descent (SD) algorithm 
       * to calculate the fuzzy communities membership of each vertex in the
       * network. The local line search is performed using the Brent's method.
       *  
       * @param n Number of threads to use in shared memory (OpenMP) 
       *        parallelism
       */
      void findCommunities_SD_Brent_parallel(int n){
        bool continueIteration = true, stopIter;
        int n_max_iterations, i, j;
        double** h = new double* [N];
        double alpha, aux, maxUij, uLimit, fn, fOld, fLimit;
        
        nThreads = n;  // Number of threads a usar
        
        // Allocating memory for the matrices
        #pragma omp parallel for num_threads (nThreads)
        for (i = 0; i < N; i++){
          h[i] = new double [c];
        }

        // General constants
        n_max_iterations =  300;  // Max number of iterations allowed
        uLimit = 0.001;           // Limit for variation of U matrix elements
        fLimit = 0.1;             // Limit for D function variation
        
        fOld = ZERO;              // Initial value for previous function

        // Initializing U
        initializeU_parallel();

        nF = 0; // Number of function evaluations in line searches
        nIterations = 0;
        fn = function_parallel();

        
        // Initializing Conjugate Gradient
        
        #pragma omp parallel for num_threads (nThreads)
        for (i = 0; i < N; i++){
          for (j =0; j < c; j++) {
            h[i][j] = -dulk(i, j); // Computing gradient on the fly
          }
        }
        

        stopIter = false; 
        while(continueIteration){      

          // Computing minimum along the search direction h
          alpha = alpha_Brent_parallel(h);
          

          // Updating coordinates
          maxUij = ZERO;
          
          #pragma omp parallel for num_threads (nThreads) \
                  private (j, aux) reduction (max: maxUij)
          for(i = 0; i < N; i++){
            for(j = 0; j < c; j++){
              aux = U[i][j];
              U[i][j] +=  alpha * h[i][j];
              aux = fabs(U[i][j] - aux);
              if (aux > maxUij) maxUij = aux;
            }
          }
  
          if (maxUij < uLimit) stopIter = true;

          if(stopIter || (fabs(fOld - fn) < fLimit) || 
            nIterations >= n_max_iterations){   // End conditions
             continueIteration = false;
          } else {      
            fOld = fn;
            fn = function_parallel();


              #pragma omp parallel for num_threads (nThreads) private (j) 
              for (i = 0; i < N; i++){
                for (j = 0; j < c; j++) {
                  h[i][j] = -dulk(i, j); // Computing gradient on the fly
                }
              }

            // Increasing counter
            nIterations++; 
          }
        }
        
        // Deallocating auxiliary array
        #pragma omp parallel for
        for (i = 0; i < N; i++){
          delete [] h[i];
        }
        delete [] h;
      
      }    

      
      
       /**
       * Finding communities in parallel using conjugate gradient in a 
       * BA-benchmark network. The method applies the conjugate gradient (CG) 
       * algorithm to calculate the fuzzy communities membership of each vertex 
       * in the network. The local line search is performed using the 
       * Brent's method.
       *  
       * @param n Number of threads to use in shared memory (OpenMP) 
       *        parallelism
       */
      void findCommunities_benchmark_Brent_parallel(int n){
        bool continueIteration = true, stopIter;
        int n_max_iterations, i, j;
        double** g = new double* [N];
        double** h = new double* [N];
        double alpha, gNum, gDen, gamma, 
               aux, maxUij, uLimit, fn, fOld, fLimit;
        
        nThreads = n;  // Number of threads a usar
        
        // Allocating memory for the matrices
        #pragma omp parallel for num_threads (nThreads)
        for (i = 0; i < N; i++){
          g[i] = new double [c];
          h[i] = new double [c];
        }

        // General constants
        n_max_iterations =  300;  // Max number of iterations allowed
        uLimit = 0.001;           // Limit for variation of U matrix elements
        fLimit = 0.1;             // Limit for D function variation
        
        fOld = ZERO;              // Initial value for previous function

        // Initializing U
        initializeU_benchmark_parallel();

        nF = 0; // Number of function evaluations in line searches
        nIterations = 0;
        fn = function_parallel();

        
        // Initializing Conjugate Gradient
        
        #pragma omp parallel for num_threads (nThreads)
        for (i = 0; i < N; i++){
          for (j =0; j < c; j++) {
            aux = dulk(i, j); // Computing gradient on the fly
            g[i][j] = h[i][j] = -aux; // Search direction
          }
        }
        

        stopIter = false; 
        while(continueIteration){      

          // Computing minimum along the search direction h
          alpha = alpha_Brent_parallel(h);
          

          // Updating coordinates
          maxUij = ZERO;
          
          #pragma omp parallel for num_threads (nThreads) \
                  private (j, aux) reduction (max: maxUij)
          for(i = 0; i < N; i++){
            for(j = 0; j < c; j++){
              aux = U[i][j];
              U[i][j] +=  alpha * h[i][j];
              aux = fabs(U[i][j] - aux);
              if (aux > maxUij) maxUij = aux;
            }
          }
  
          if (maxUij < uLimit) stopIter = true;

          if(stopIter || (fabs(fOld - fn) < fLimit) || 
            nIterations >= n_max_iterations){   // End conditions
             continueIteration = false;
          } else {      
            fOld = fn;
            fn = function_parallel();


            // Computing gamma_i (Polak-Ribière)
            gNum = gDen = ZERO;
            
            #pragma omp parallel num_threads (nThreads)
            {
              #pragma omp for private (j) reduction(+: gDen) reduction(-: gNum)
              for (i = 0; i < N; i++){
                for (j = 0; j < c; j++) {
                  grad[i][j] = dulk(i, j); // Computing gradient on the fly
                  gNum -= (-grad[i][j] - g[i][j]) * grad[i][j];
                  gDen += g[i][j] * g[i][j];
                }
              }

              #pragma omp single
              gamma = gNum / gDen; 

              #pragma omp for private (j)
              for (i = 0; i < N; i++){
                for (j = 0; j < c; j++) {
                  g[i][j] = -grad[i][j];                // Updating g
                  h[i][j] = g[i][j] + gamma * h[i][j];  // Updating h
                }
              }
          }
            // Increasing counter
            nIterations++; 
          }
        }
        
        // Deallocating auxiliary arrays
        #pragma omp parallel for
        for (i = 0; i < N; i++){
          delete [] g[i];
          delete [] h[i];
        }
        delete [] g;
        delete [] h;
      
      }    

      
      /**
       * This method checks the correctness of results when using the
       * BA-communities network benchmark
       * @return True if the result is correct, false otherwise
       */
      bool checkU_benchmark() {         
        bool correct = true;
        int i, j, k, m, sizeBA;
        double aux, aux1, m1, m2, limit0, limit1;
        limit0 = 0.01;  // Limit for the central node
        limit1 = 0.01;  // Limit1 is the error in memberships


        // Checking that the memberships of central node equal 1/(c + 1) 
        // Remember: c is number of communities -1
        if (correct) {
          for (i = 0; i < c && correct; i++) {
            correct = fabs(U[N-1][i] - ONE /(c + 1)) <= limit0;
          }
        }
      
        // Checking the equivalent nodes in different communities have the same
        // membership
        double a1, a2;

        if (correct){
          sizeBA = (N - 1) / (c + 1); // Size of each one of the c BA subnetworks
          for (i = 0; i < c && correct; i++){
            for (j = i + 1; j < c && correct; j++){
              for (k = 0; k < sizeBA && correct; k++){
                a1 = a2 = ZERO;
                m1 = m2 = ZERO;
                for (m = 0; m < c && correct; m++){
                  aux = U[k + i * sizeBA][m];
                  a1 += aux;
                  aux1 = U[k + j * sizeBA][m];
                  a2 += aux1;
                  m1 += aux  * aux;
                  m2 += aux1 * aux1;
                }
                m1 = m1 + (ONE - a1) * (ONE - a1);
                m2 = m2 + (ONE - a2) * (ONE - a2);
                m1 = sqrt(m1/(c + 1));  //RMS of the whole vector
                m2 = sqrt(m2/(c + 1));  //RMS of the whole vector
                correct = fabs((m1-m2)) <= limit1; 
                
              }
            }
          }
        }
        return correct;
      }


      /**
       * Writes the U matrix to file fileName. The first number in the first 
       * line of the file is the number of nodes in the network. The second 
       * number in the first line is the number of communities. Each remaining 
       * line contains the ID number of the corresponding node and the 
       * memberships of that node to the communities in the network.
       * 
       * @param fileName The file where the U matrix is written.
       */
      void writeU(string fileName){
        int i, j;
        double aux, uij;
        ofstream file;

        file.open(fileName);
        file << N << " " << (c + 1) << endl;  // Writting number of nodes and
                                              // number of communities
        
        // Setting 6 significant digits after the decimal point for 
        // output of floating-point values
        file.precision(6);
        file << fixed;
        
        // Writing a line per node
        for (i = 0; i < N; i++){
          aux = ZERO;
          file << setw (10) << i << " ";
          for (j = 0; j < c; j++) { // Here c equals number of communities - 1
            uij = U[i][j];
            file << setw(10) << uij << " ";
            aux += uij;
          }
          file << setw(10) 
               << (ONE - aux) << endl; // Writing last component of vector ui
        }
        file.close();
        
      }
    
      
      
      /**
       * Returns the membership matrix
       * @return The membership matrix
       */
      double** getU(){
        return U;
      }
      

      /**
       * Returns the number of function calls made in the alpha_exact method
       * @return 
       */
      int numberFunctionCalls() {
        return nF;
      }
      

      /**
       * Returns the number of iterations in the minimization method
       * @return 
       */
      int numberIterations() {
        return nIterations;
      }
      

      /**
       * Accessor method returning the number of nodes (vertices) in the network
       * @return  The number of nodes (vertices) in the network
       */
      int getN(){
        return N;
      }
      
  
      /**
       * Accessor method returning the number of edges in the network
       * @return  The number of edges in the network
       */
      int getM(){
        return M;
      }
      
      
      /**
       * Mutator method to update the number of communities
       * @param nc The new number of communities
       */
      void setNComm(int nc) {
        
        // Deallocating auxiliary arrays
        for (int i = 0; i < N; i++){
          delete [] U[i];
          delete [] grad[i];
        }
        delete [] U;
        delete [] grad;
      
        c = nc - 1; // Updating number of communities
      
        // Allocating memory for the nuew communities
        U = new double* [N];
        grad = new double* [N];
        for (int i = 0; i < N; i++){
          U[i] = new double [c];
          grad[i] = new double [c];
        }
      }
      
      
      /**
       * Destructor method. Releases memory from arrays A, degree, U and grad.
       */
      ~FuzzyCommunities(){
        // Deallocating memory
        delete [] A;
        delete [] degree;
        
        // Deallocating auxiliary arrays
        for (int i = 0; i < N; i++){
          delete [] U[i];
          delete [] grad[i];
        }
        
        delete [] U;
        delete [] grad;
      }
      
      
      
    private:

      /**
       * Auxiliary method to read a symmetric network from a file stream as an
       * edge list. The method determines if the number of tokens in the edge 
       * lines is two (origin and destination nodes) or three (origin node, 
       * destination node, and edge weight). If only two tokens are found
       * (origin and destination nodes) the default weight for the edges is set
       * to 1.0.
       * 
       * @param file The file stream connected to the data file
       * @param flag Flag to indicate the kind of reading. flag = true when the 
       *        file contains number of nodes and edges in the first line and
       *        flag = false if it does not.
       */
      void readFuzzyCommunities(ifstream &file, bool flag){
        int i, j, nt;
        double s;
        string line;
        string tokens [3];

        // Reading edge list
        A = new SparseArray(M);
        degree = new int[N];
        for (int i = 0; i < N; i++) degree[i] = 0; // Initializing degrees

        // Reading edges FOR A SYMMETRIC ADJACENCY MATRIX!!!!!!!!
        if (flag) getline(file, line); // Reading till end of first line (
                                       // the one with N and M)
        getline(file, line); // Reading the first edge
        stringstream readL(line);
        s = 1.0;             // Default weight value
        nt = 0;              // Number of tokens
        
        // Determining the number of tokens using the first edge.        
        // WARNING. The following loop uses short-circuit evaluation of the
        // loop condition. This means that the nt<=2 MUST be the first
        // condition to evaluate
        while(nt <= 2 && readL >> tokens[nt]) {  // A maximum of three 
          nt++;                                  // tokens are read
        }
        i = stoi(tokens[0]);              // String to int in C++11
        j = stoi(tokens[1]);              // String to int in C++11
        if (nt == 3) s = stod(tokens[2]); // String to double in C++11
        
        if (i >= N || j >= N) {
          cout << "Error. Reading edge " << i << "--" << j << endl
               << " There are " << N << " nodes in the network and the"
               << " first one is labeled as 0." 
               << "No node can have an index > N-1.";
          exit(1);
        }
        
        // Adding edge to the edge list
        if (i <= j) {     
          A->put(i, j, s);
        } else {
          A->put(j, i, s);
        }
        
        degree[i]++;
        degree[j]++;
        
        // Reading the remaining edges
        for (int k = 1; k < M; k++){
          file >> i;   // First node
          file >> j;   // Second node
          
          if (nt == 3) file >> s;   // Reading weight of the edge
          
          if (i >= N || j >= N) {
            cout << "Error. Reading edge " << i << "--" << j << endl
                 << " There are " << N << " nodes in the network and the"
                 << " first one is labeled as 0." 
                 << "No node can have an index > N-1.";
            exit(1);
          }
          if (i <= j) {     
            A->put(i, j, s);
          } else {
            A->put(j, i, s);
          }
          degree[i]++;
          degree[j]++;
        }

        // Allocating memory for the matrices when the number of communities 
        // have been provided in the constructor.
        if (c != 0) {
          U = new double* [N];
          grad = new double* [N];
          for (int i = 0; i < N; i++){
            U[i] = new double [c];
            grad[i] = new double [c];
          }
        }
      }
           
       
      /**
       * Calculates the initial U membership matrix used in the algorithm using 
       * a Dirichlet distribution. This method uses a Gamma distribution 
       * to generate c (total number of communities) random numbers with 
       * shape and scale parameters equal to 1.0. Then, each result is divided 
       * by the sum of the numbers to generate a normalized Dirichlet 
       * distribution with order c and alpha parameters equal to 1.0. For 
       * details see: 
       * 
       * - Bela A. Frigyik, Amol Kapila, and Maya R. Gupta. Introduction to the 
       * Dirichlet Distribution and Related Processes. Department of Electrical 
       * Engineering University of Washington. UWEE Technical Report Number 
       * UWEETR-2010-0006. December 2010. 
       *       https: isocpp.org/files/papers/n3551.pdf
       */
      void initializeU(){
        int i, j;
        double r, aux;
        
        // Generating seed for the random number generator
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        
        default_random_engine g(seed);              // Random numbers engine
        gamma_distribution<double> gDist (1.0,1.0); // Gamma distribution
        
        for (i = 0; i < N; i++) {
          aux = ZERO;
          for (j = 0; j < c; j++) {  // Remember that c here equals number of 
            r = gDist(g);            // communities minus one
            U[i][j] = r;
            aux += r;
          }
          aux += gDist(g); // Additional random number for the last community
          aux = ONE / aux; // Computing inverse just once
          
          // Normalizing the c-1 independent coefficients
          for (j = 0; j < c; j++) {
            U[i][j] *= aux; 
          }
        }
      }
    
   
      /**
       * Calculates the initial membership matrix U used in the algorithm using 
       * a Dirichlet distribution. This method uses a Gamma distribution 
       * to generate c (total number of communities) random numbers with 
       * shape and scale parameters equal to 1.0. Then, each result is divided 
       * by the sum of the numbers to generate a 
       * normalized Dirichlet distribution with order c and alpha parameters 
       * equal to 1.0. For details see: 
       * 
       * - Bela A. Frigyik, Amol Kapila, and Maya R. Gupta. Introduction to the 
       * Dirichlet Distribution and Related Processes. Department of Electrical 
       * Engineering University of Washington. UWEE Technical Report Number 
       * UWEETR-2010-0006. December 2010. 
       *       https: isocpp.org/files/papers/n3551.pdf
       */
      void initializeU_parallel(){
        int i, j;
        double r, aux;
        
        // Generating seed for the random number generator
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        
        default_random_engine g(seed);              // Random numbers engine
        gamma_distribution<double> gDist (1.0,1.0); // Gamma distribution
        
        
        #pragma omp parallel for num_threads (nThreads)
        for (i = 0; i < N; i++) {
          aux = ZERO;
          for (j = 0; j < c; j++) {  // Remember that c here equals number of 
            r = gDist(g);            // communities minus one
            U[i][j] = r;
            aux += r;
          }
          aux += gDist(g); // Additional random number for the last community
          
          // Normalizing the c-1 independent coefficients
          for (j = 0; j < c; j++) {
            U[i][j] /= aux; 
          }
        }
      }

      
      /**
       * Calculates the initial membership matrix U used in the algorithm for a 
       * BA-benchmark network with c BA_subnetworks (communities).
       * The membership of the central (last) node is set to a=0.5 for the first 
       * community and to (1-a)/(c-1) for the rest. In the BA-subnetworks we use 
       * uik = a for all the nodes, i, in subnetwork k and (1-a)/(c-1) for the 
       * remaining memberships. 
       */
      void initializeU_benchmark(){
        int i, j, k, size;
        double aux, aux1;
        
        // Initializing central node
        aux = 0.5;
        aux1 = (ONE - aux) / c; // Remember: c is number of communities -1
        U[N -1][0] = aux;
        for (i = 1; i < c; i++) U[N-1][i] = aux1;

        // Initializing rest of nodes
        size = (N - 1) / (c + 1);      // Number of nodes per community
        for (i = 0; i < c; i++){       // Running on the first c communities
          for (j = 0; j < size; j++){  // Running on the nodes of the community
            U[j + size * i][i] = aux;
            for (k = 0; k < c; k++){
              if (k != i) U[j + size * i][k] = aux1;
            }
          }
        }

        for (j = 0; j < size; j++){// Running on the nodes of the last community
          for (k = 0; k < c; k++){
            if (k != i) U[j + size * c][k] = aux1;
          }
        }

      }
      
      
       /**
       * Calculates in parallel the initial membership matrix U used in the 
       * algorithm for a BA-benchmark network with c BA_subnetworks 
       * (communities). The membership of the central (last) node is set to 
       * a=0.5 for the first community and to (1-a)/(c-1) for the rest. In 
       * the BA-subnetworks we use uik = a for all the nodes, i, in 
       * subnetwork k and (1-a)/(c-1) for the remaining memberships. 
       */
      void initializeU_benchmark_parallel(){
        int i, j, k, size;
        double aux, aux1;
        
        // Initializing central node
        aux = 0.5;
        aux1 = (ONE - aux) / c; // Remember: c is number of communities -1
        U[N -1][0] = aux;
        for (i = 1; i < c; i++) U[N-1][i] = aux1;

        // Initializing rest of nodes
        size = (N - 1) / (c + 1);      // Number of nodes per community
        
        #pragma omp parallel num_threads (nThreads)
        {
          #pragma omp for private (j, k) collapse (2)
          for (i = 0; i < c; i++){       // Running on the first c communities
            for (j = 0; j < size; j++){  // Running on the nodes of the community
              U[j + size * i][i] = aux;
              for (k = 0; k < c; k++){
                if (k != i) U[j + size * i][k] = aux1;
              }
            }
          }

          #pragma omp for private (k) collapse (2)
          for (j = 0; j < size; j++){// Running on the nodes of the last community
            for (k = 0; k < c; k++){
              if (k != i) U[j + size * c][k] = aux1;
            }
          }
        }
      }
      
       
       /**
        * Calculates alpha constant (exact line search) along a search 
        * direction in networks using Brent's method. Alpha is a small 
        * step size constant used to calculate the next membership matrix 
        * in the iteration. This method uses a line minimization implementing 
        * Brent's method which brackets the minimum by the golden rule and 
        * locates the minimum  using inverse parabolic interpolation or, if 
        * it does not work, golden ratio search.
        * 
        * @param h     The current search direction
        * @return      The alpha parameter minimizing the function along the 
        *              search direction
        */
      double alpha_Brent(double** h){
        double ax, bx, cx, x, fa, fb, fc, fx, aux, tol, golden = 1.6180339887,
               CGOLD, d, e, eps, xm, p, q, r, tol1, t2, u, v, w, fu, fv , fw, 
               tol3, a, b;

        tol = 3e-8; //Square root of machine double precision
        int t, tmax = 100;
        CGOLD = HALF * (3.0 - sqrt(5.0));
        d = ZERO;

        /* Bracketing the minimum using golden ratio */
        // Initializing
        ax = ZERO;
        bx = ONE;
        fa = function_s(ax, h);
        fb = function_s(bx, h);

        if (fb > fa) {  // We always want fb < fa (we always 
          aux = ax;     // search in the decreasing direction)
          ax = bx;
          bx = aux;
          aux = fa;
          fa = fb;
          fb = aux;
        }

        cx = bx + golden * (bx-ax); // Increasing by the golden ratio of the
                                    // a-b segment
        fc = function_s(cx, h) ;

        while (fc < fb) { // Searching until the function increases
          ax = bx;
          fa = fb;
          bx = cx;
          fb = fc;
          cx = bx + golden * (bx-ax);
          fc = function_s(cx, h);
        }

        /* Minimum bracketed */

        // Computing minimum applying Brent's approach

        a = (ax < cx ? ax : cx);
        b = (ax > cx ? ax : cx);

        eps = 1.2e-16;       //Machine precision
        tol1 = eps + ONE;
        eps = sqrt(eps);

        x = w = v = bx;
        e = ZERO;
        fw = fv = fx = function_s(x, h);
        tol3 = tol / 3.0;

        xm = HALF * (a + b);
        tol1 = eps * fabs(x) + tol3;
        t2 = TWO * tol1;

        // Main loop
        for(t = 0; fabs(x - xm) > (t2 - HALF * (b - a)) && t < tmax; t++){
          p = q = r = ZERO;

          if(fabs(e) > tol1){ // Fitting a parabola
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = TWO * (q - r);
            if(q > ZERO){
              p = -p;
            } else {
              q = -q;
            }
            r = e;
            e = d;
          }

          if ((fabs(p) < fabs(HALF * q * r)) && (p > q * (a - x))
              && (p < q * (b - x))) {  // Parabolic interpolation step
            d = p / q;
            u = x + d;
            if (((u - a) < t2) || ((b - u) < t2)) { // f must not be evaluated
                                                    // too close to a or b
              d = tol1;
              if (x >= xm) {
                d = -d;
              }
            }
          } else { // Golden-section step        
            if (x < xm) {
              e = b - x;
            } else {
              e = a - x;
            }
            d = CGOLD * e;
          }

          if (fabs(d) >= tol1) { // f must not be evaluated too close to x
              u = x + d;
          } else if(d > ZERO) { 
              u = x + tol1;
          } else {
              u = x - tol1;
          }
          fu = function_s(u, h);

          if (fu <= fx) {
            if (u >= x){
              a = x;
            } else {
              b = x;
            }
            v = w;
            fv = fw;
            w = x;
            fw = fx;
            x = u;
            fx = fu;
          } else {
            if (u < x) {
              a = u;
            } else {
              b = u;
            }
            if (fu <= fw || w == x) {
              v = w;
              w = u;
              fv = fw;
              fw = fu;
            } else if (fu <= fv || v == x || v == w) {
              v = u;
              fv = fu;
            }
          }

          xm = HALF * (a + b);
          tol1 = eps * fabs(x) + tol3;
          t2 = TWO * tol1;
        } // End main loop

        if(t == tmax){
          throw "ERROR. Max number of iterations allowed in Brent exceeded";
        }

        return x;
      }    
       
    
                
       /**
        * Calculates alpha constant along a search 
        * direction in symmetric networks using the Newton-Raphson method.
        * Alpha is a small step size constant used to calculate the next 
        * partition in the iteration. The method uses a line minimization 
        * implementing the Newton-Raphson method.
        * 
        * @param h     The current search direction
        * @return      The alpha parameter minimizing the function along the 
        *              search direction
        */
      double alpha_NR(double** h){
        double fpa, fpb, a, b, ba, bc, fbc, fba, alpha, alpha1, tol, 
               delta, sqEps;
        int nIter, maxIter;
        bool goOn = true;
        
        // Initializing variables
        delta = 0.001;          // Initial increment for line search
        tol = 1.0E-3;           // Convergence limit (3 significant digits)
        nIter = 0;              // Bracketing iterations counter
        maxIter = 100;          // Maximum number of iterations
        a= ZERO;                // Initial alpha value
        fpa = dfa(a, h);        // Derivative of the function
        sqEps = 2.0E-8;         // Square root of machine double precision
        
        
        // Bracketing root
        while (goOn && maxIter) {
          nIter++;
          
          if (fpa < 0) {
            b = a + delta;
          } else {
            b = b - delta;
          }
          
          fpb = dfa(b, h);
          
          if (fpa * fpb < 0 ){
            goOn = false;
          } else {
            a = b;
            fpa = fpb;
            delta = TWO * delta; // Exponential movement
          }          
        }
        
        // First approach to the root of the derivative by linear interpolation
        alpha1 = (fpa * b - fpb * a) / (fpa - fpb);
        //cout << endl << "a, b, alphaIni: " << a << "  " << b << "  " << alpha1 << endl;
        //cout <<endl<< "a, b, alphaIni, fpa, fpb: " << a << "  "<<b<<"  "<<alpha1<<"  "<<fpa<<"  "<<fpb<<endl;
        
        
        // Locating minimum with Newton-Raphson method
        alpha = INFINITY;
        nIter = 0;
        while (fabs(alpha - alpha1)/fmin(alpha, alpha1) > tol && 
               fabs(alpha - alpha1) > sqEps && nIter < maxIter){
          nIter++;
          alpha = alpha1;
          alpha1 = alpha - dfa(alpha, h) / d2fa2(alpha, h);
          //cout << "alpha1: " << alpha1 << "  " << dfa(alpha, h) << "  " << d2fa2(alpha, h) << endl;
        }
        if (nIter == maxIter) alpha1 = 1.0E-4; // Small value;
      
        return alpha1;
          
      }
      

      
       /**
        * Calculates alpha constant in parallel along a search 
        * direction in symmetric networks using the Newton-Raphson method.
        * Alpha is a small step size constant used to calculate the next 
        * partition in the iteration. The method uses a line minimization 
        * implementing the Newton-Raphson method.
        * 
        * @param h     The current search direction
        * @return      The alpha parameter minimizing the function along the 
        *              search direction
        */
      double alpha_NR_parallel(double** h){
        double fpa, fpb, a, b, ba, bc, fbc, fba, alpha, alpha1, tol, 
               delta, sqEps;
        int nIter, maxIter;
        bool goOn = true;
        
        // Initializing variables
        delta = 0.001;            // Initial increment for line search
        tol = 1.0E-3;             // Convergence limit (3 significant digits)
        nIter = 0;                // Bracketing iterations counter
        maxIter = 100;            // Maximum number of iterations
        a= ZERO;                  // Initial alpha value
        fpa = dfa_parallel(a, h); // Derivative of the function
        sqEps = 2.0E-8;           // Square root of machine double precision
        
        
        // Bracketing root
        while (goOn && maxIter) {
          nIter++;
          
          if (fpa < 0) {
            b = a + delta;
          } else {
            b = b - delta;
          }
          
          fpb = dfa_parallel(b, h);
          
          if (fpa * fpb < 0 ){
            goOn = false;
          } else {
            a = b;
            fpa = fpb;
            delta = TWO * delta; // Exponential movement
          }          
        }
        
        
        // First approach to the root of the derivative by linear interpolation
        alpha1 = (fpa * b - fpb * a) / (fpa - fpb);
        //cout << endl << "a, b, alphaIni: " << a << "  " << b << "  " << alpha1 << endl;
        
        
        // Locating minimum with Newton-Raphson method
        alpha = INFINITY;
        nIter = 0;
        while (fabs(alpha - alpha1)/fmin(alpha, alpha1) > tol && 
               fabs(alpha - alpha1) > sqEps && nIter < maxIter){
          nIter++;
          alpha = alpha1;
          alpha1 = alpha - dfa_parallel(alpha, h) / d2fa2_parallel(alpha, h);
          //cout << "alpha: " << alpha << "  " << dfa(alpha, h) << "  " << d2fa2(alpha, h) << endl;
        }
        if (nIter == maxIter) alpha1 = 1.0E-4; // Small value
      
        return alpha1;
          
      }
      
    
      
       /**
        * Calculates alpha constant (exact line search) along a search 
        * direction in networks using Brent's method. Alpha is a small 
        * step size constant used to calculate the next membership matrix 
        * in the iteration. This method uses a line minimization implementing 
        * Brent's method which brackets the minimum by the golden rule and 
        * locates the minimum  using inverse parabolic interpolation or, if 
        * it does not work, golden ratio search.
        * 
        * @param h     The current search direction
        * @return      The alpha parameter minimizing the function along the 
        *              search direction
        */
      double alpha_Brent_parallel(double** h){
        double ax, bx, cx, x, fa, fb, fc, fx, aux, tol, golden = 1.6180339887,
               CGOLD, d, e, eps, xm, p, q, r, tol1, t2, u, v, w, fu, fv , fw, 
               tol3, a, b;

        tol = 3e-8; //Square root of machine double precision
        int t, tmax = 100;
        CGOLD = HALF * (3.0 - sqrt(5.0));
        d = ZERO;

        /* Bracketing the minimum using golden ratio */
        // Initializing
        ax = ZERO;
        bx = ONE;
        fa = function_s_parallel(ax, h);
        fb = function_s_parallel(bx, h);

        if (fb > fa) {  // We always want fb < fa (we always 
          aux = ax;     // search in the decreasing direction)
          ax = bx;
          bx = aux;
          aux = fa;
          fa = fb;
          fb = aux;
        }

        cx = bx + golden * (bx-ax); // Increasing by the golden ratio of the
                                    // a-b segment
        fc = function_s_parallel(cx, h) ;

        while (fc < fb) { // Searching until the function increases
          ax = bx;
          fa = fb;
          bx = cx;
          fb = fc;
          cx = bx + golden * (bx-ax);
          fc = function_s_parallel(cx, h);
        }

        /* Minimum bracketed */

        // Computing minimum applying Brent's approach

        a = (ax < cx ? ax : cx);
        b = (ax > cx ? ax : cx);

        eps = 1.2e-16;       //Machine precision
        tol1 = eps + ONE;
        eps = sqrt(eps);

        x = w = v = bx;
        e = ZERO;
        fw = fv = fx = function_s_parallel(x, h);
        tol3 = tol / 3.0;

        xm = HALF * (a + b);
        tol1 = eps * fabs(x) + tol3;
        t2 = TWO * tol1;

        // Main loop
        for(t = 0; fabs(x - xm) > (t2 - HALF * (b - a)) && t < tmax; t++){
          p = q = r = ZERO;

          if(fabs(e) > tol1){ // Fitting a parabola
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = TWO * (q - r);
            if(q > ZERO){
              p = -p;
            } else {
              q = -q;
            }
            r = e;
            e = d;
          }

          if ((fabs(p) < fabs(HALF * q * r)) && (p > q * (a - x))
              && (p < q * (b - x))) {  // Parabolic interpolation step
            d = p / q;
            u = x + d;
            if (((u - a) < t2) || ((b - u) < t2)) { // f must not be evaluated
                                                    // too close to a or b
              d = tol1;
              if (x >= xm) {
                d = -d;
              }
            }
          } else { // Golden-section step        
            if (x < xm) {
              e = b - x;
            } else {
              e = a - x;
            }
            d = CGOLD * e;
          }

          if (fabs(d) >= tol1) { // f must not be evaluated too close to x
              u = x + d;
          } else if(d > ZERO) { 
              u = x + tol1;
          } else {
              u = x - tol1;
          }
          fu = function_s_parallel(u, h);

          if (fu <= fx) {
            if (u >= x){
              a = x;
            } else {
              b = x;
            }
            v = w;
            fv = fw;
            w = x;
            fw = fx;
            x = u;
            fx = fu;
          } else {
            if (u < x) {
              a = u;
            } else {
              b = u;
            }
            if (fu <= fw || w == x) {
              v = w;
              w = u;
              fv = fw;
              fw = fu;
            } else if (fu <= fv || v == x || v == w) {
              v = u;
              fv = fu;
            }
          }

          xm = HALF * (a + b);
          tol1 = eps * fabs(x) + tol3;
          t2 = TWO * tol1;
        } // End main loop

        if(t == tmax){
          throw "ERROR. Max number of iterations allowed in brent method exceeded";
        }

        return x;
      }    

      
    
      /**
       * Computing gradient vector for the current value of the U membership 
       * matrix.
       */
      inline void gradient(){
        int l, k;
        
        for(l = 0; l < N; l++){
          for(k = 0; k < c; k++){
            grad[l][k] = dulk(l, k);   // derivative 
          } 
        }
      }
       

      
      /**
       * Computing numerically the gradient vector for the current value of the 
       * U membership matrix.
       */
      inline void gradientNum(){
        int l, k;
        double aux, fmas, fmenos, delta, twoDeltaInv;
        delta = 0.0001;    // Increment for numerical derivatives
        twoDeltaInv = ONE / (TWO * delta);
        
        for(l = 0; l < N; l++){
          for(k = 0; k < c; k++){
            aux = U[l][k];
            U[l][k] = aux + delta;
            fmas = function();
            U[l][k] = aux - delta;
            fmenos = function();
            U[l][k] = aux;
            grad[l][k] = (fmas - fmenos) * twoDeltaInv;   // derivative 
          } 
        }
      }
            
 
      
      /**
       * Parallel computing of the gradient vector for the current value 
       * of the U membership matrix.
       */
      inline void gradient_parallel(){
        int l, k;
        
        #pragma omp parallel for num_threads (nThreads) private (k)
        for(l = 0; l < N; l++){
          for(k = 0; k < c; k++){
            grad[l][k] = dulk(l, k);   // derivative 
          } 
        }
      }
  
     
      
      /**
       * Method that computes the function used in the calculation 
       * of alpha along the search direction. This version is for an undirected,
       * symmetric network.
       * @param alpha     The current alpha value
       * @param h         The current search direction
       * @return          The value of the function for the given input
       */
      inline double function_s (double alpha, double** h){
        double f, aux, aux1;
        int j;
        
        nF++; // Increasing the counter of function computation calls
              
  
        f = aux1 = ZERO;

        for (int i = 0; i < N; i++){
          aux = siih(i, alpha, h);
          f += aux * aux;
          
          for (j = i + 1; j < N; j++){
            aux = sijh(i, j, alpha, h);
            aux1 += aux * aux;
          }
        }

        f += TWO * (aux1 + A->get_weightsSum());
        
        aux1 = ZERO;
        for (j = 0; j < M; j++){
          aux1 += sijh(A->r(j), A->c(j), alpha, h) * A->v(j);
        }
        f -= FOUR * aux1;
        return f;
        
      }
     

      
      /**
       * Function that calculates in parallel the function used in the 
       * calculation of alpha along the search direction. This version is for 
       * an undirected, symmetric, network.
       * @param alpha     The current alpha value
       * @param h         The current search direction
       * @return          The value of the function for the given input
       */
      inline double function_s_parallel (double alpha, double** h){
        double f, aux, aux1;
        int j;
        
        nF++; // Increasing the counter of function computation calls
              
  
        f = aux1 = ZERO;
        
        #pragma omp parallel num_threads (nThreads)
        {
          #pragma omp for private (j, aux) reduction(+: f, aux1)
          for (int i = 0; i < N; i++){
            aux = siih(i, alpha, h);
            f += aux * aux;
            for (j = i + 1; j < N; j++){
              aux = sijh(i, j, alpha, h);
              aux1 += aux * aux;
            }
          }
          
          #pragma omp single
          {
            f += TWO * (aux1 + A->get_weightsSum());
            aux1 = ZERO;
          }
          #pragma omp for reduction(+: aux1)
          for (j = 0; j < M; j++){
            aux1 += sijh(A->r(j), A->c(j), alpha, h) * A->v(j);
          }
        }
        f -= FOUR * aux1;
        return f;
        
      }
   
        
      
      /**
       * Function that calculates the topological error functional D.
       * @return  The value of the function for the current U membership matrix
       */
      inline double function(){
        double f, aux, aux1;
        int j;
        
        nF++; // Increasing the counter of function computation calls
        
        f = aux1 = ZERO;
  
        for (int i = 0; i < N; i++){
          aux = sii(i);
          f += aux * aux;
          
          for (j = i + 1; j < N; j++){
            aux = sij(i, j);
            aux1 += aux * aux;
          }
        }
        f += TWO * (aux1 + A->get_weightsSum());
        
        aux1 = ZERO;
        for (j = 0; j < M; j++){
          aux1 += sij(A->r(j), A->c(j)) * A->v(j);
        }
        f -= FOUR * aux1;
        return f;
      }
 
        
      
      /**
       * Function that calculates the derivative of the topological error 
       * functional D with respect to alpha for the current search 
       * direction.
       * @param alpha  The current alpha value
       * @param h      The current search direction
       * @return       The value of the derivative for the current search 
       *               direction
       */
      inline double dfa(double alpha, double** h){
        double f, aux;
        int j;
                
        f = aux = ZERO;
  
        for (int i = 0; i < N; i++){
          aux += siih(i, alpha, h) * dsiia(i, alpha, h);
          
          for (j = i + 1; j < N; j++){
            f += sijh(i, j, alpha, h) * dsija(i, j, alpha, h);
          }
        }
        
        for (j = 0; j < M; j++){
          f -= dsija(A->r(j), A->c(j), alpha, h) * A->v(j);
        }
                
        return TWO * aux + FOUR * f;
      }

         
      /**
       * Function that calculates in parallel the derivative of the topological  
       * error functional D with respect to alpha for the current search 
       * direction.
       * @param alpha  The current alpha value
       * @param h      The current search direction
       * @return       The value of the derivative for the current search 
       *               direction
       */
      inline double dfa_parallel(double alpha, double** h){
        double f, aux;
        int j;
                
        f = aux = ZERO;
  
        #pragma omp parallel num_threads (nThreads)
        {
          #pragma omp for reduction(+: f, aux) private (j)
          for (int i = 0; i < N; i++){
            aux += siih(i, alpha, h) * dsiia(i, alpha, h);

            for (j = i + 1; j < N; j++){
              f += sijh(i, j, alpha, h) * dsija(i, j, alpha, h);
            }
          }
        
        
          #pragma omp for reduction(-: f)
          for (j = 0; j < M; j++){
            f -= dsija(A->r(j), A->c(j), alpha, h) * A->v(j);
          }
        }
                
        return TWO * aux + FOUR * f;
      }
      
  
      
      /**
       * Function that calculates the second derivative of the topological error 
       * functional D with respect to alpha for the current search direction.
       * @param alpha  The current alpha value
       * @param h      The current search direction
       * @return       The value of the derivative for the current search 
       *               direction
       */
      inline double d2fa2(double alpha, double** h){
        double f, aux, aux1;
        int j;
                
        f = aux1 = ZERO;
  
        for (int i = 0; i < N; i++){
          aux = dsiia(i, alpha, h);
          aux1 += aux * aux;
          aux1 += siih(i, alpha, h) * d2siia2(i, h);
          
          for (j = i + 1; j < N; j++){
            aux = dsija(i, j, alpha, h);
            f += aux * aux;
            f += sijh(i, j, alpha, h) * d2sija2(i, j, h);
          }
        }
        
        for (j = 0; j < M; j++){
          f -= d2sija2(A->r(j), A->c(j), h) * A->v(j);
        }
        
        return TWO * aux1 + FOUR * f;
      }

      
      /**
       * Function that calculates in parallel the second derivative of the 
       * topological error functional D with respect to alpha for the current 
       * search direction.
       * @param alpha  The current alpha value
       * @param h      The current search direction
       * @return       The value of the derivative for the current search 
       *               direction
       */
      inline double d2fa2_parallel(double alpha, double** h){
        double f, aux, aux1;
        int j;
                
        f = aux1 = ZERO;

        #pragma omp parallel num_threads (nThreads)
        {
          #pragma omp for reduction(+: f, aux1) private (aux, j)
          for (int i = 0; i < N; i++){
            aux = dsiia(i, alpha, h);
            aux1 += aux * aux;
            aux1 += siih(i, alpha, h) * d2siia2(i, h);
            for (j = i + 1; j < N; j++){
              aux = dsija(i, j, alpha, h);
              f += aux * aux;
              f += sijh(i, j, alpha, h) * d2sija2(i, j, h);
            }
          }
        
          #pragma omp for reduction(-: f)
          for (j = 0; j < M; j++){
            f -= d2sija2(A->r(j), A->c(j), h) * A->v(j);
          }
        }
        
        return TWO * aux1 + FOUR * f;
      }

      
      /**
       * Function that calculates the topological error functional D for a 
       * symmetric network.
       * @return  The value of the function for the current U membership matrix
       */
      inline double function_parallel(){
        double f, aux, aux1;
        int j;
        
        nF++; // Increasing the counter of function computation calls
        
        f = aux1 = ZERO;
        
        #pragma omp parallel num_threads (nThreads)
        {
          #pragma omp for private (j, aux) reduction(+: f, aux1)
          for (int i = 0; i < N; i++){
            aux = sii(i);
            f += aux * aux;

            for (j = i + 1; j < N; j++){
              aux = sij(i, j);
              aux1 += aux * aux;
            }
          }
          #pragma omp single
          {
            f += TWO * (aux1 + A->get_weightsSum());
            aux1 = ZERO;
          }
        
          #pragma omp for reduction(+: aux1)
          for (j = 0; j < M; j++){
            aux1 += sij(A->r(j), A->c(j)) * A->v(j);
          }
        }
        f -= FOUR * aux1;
        return f;
      }

  
      /**
       * Function that calculates the self-relationship for node i
       * @param i The current node index  
       * @return  The self-relationship value for node i
       */
      inline double sii(int i){
        double aux, pii, sumim;
        pii = sumim = ZERO;
        
        for (int m = 0; m < c; m++){
          aux = U[i][m];
          pii += aux * aux;
          sumim += aux;
        }
        sumim = ONE - sumim;
        return pii + sumim * sumim;
      }
      
      
      /**
       * Function that calculates the relationship between nodes i and j
       * @param i The first node index  
       * @param j The second node index  
       * @return  The relationship value between nodes i and j
       */
      inline double sij(int i, int j){
        double pij, sumim, sumjm;
        pij = sumim = sumjm = ZERO;
        
        for (int m = 0; m < c; m++){
          pij += U[i][m] * U[j][m];
          sumim += U[i][m];
          sumjm += U[j][m];
        }
        return pij + (ONE - sumim) * (ONE - sumjm);
      }
         
            
         
     
      /**
       * Function that calculates the self-relationship for node i along 
       * the search direction
       * @param i      The current node index  
       * @param alpha  The current alpha value
       * @param h      The current search direction
       * @return       The self-relationship value for node i
       */
      inline double siih(int i, double alpha, double** h){
        double aux, pii, sumim;
        pii = sumim = ZERO;
        
        for (int m = 0; m < c; m++){
          aux = U[i][m] + alpha * h[i][m];
          pii += aux * aux;
          sumim += aux;
        }
        sumim = ONE - sumim;
        return pii + sumim * sumim;
      }
             

      
      /**
       * Function that calculates the relationship between nodes i and j along 
       * the search direction
       * @param i The first node index  
       * @param j The second node index  
       * @return  The relationship value between nodes i and j
       */
      inline double sijh(int i, int j, double alpha, double** h){
        double pij, sumim, sumjm, aux1, aux2;
        pij = sumim = sumjm = ZERO;
        
        for (int m = 0; m < c; m++){
          aux1 = U[i][m] + alpha * h[i][m];
          aux2 = U[j][m] + alpha * h[j][m];
          pij += aux1 * aux2;
          sumim += aux1;
          sumjm += aux2;
        }
        return pij + (ONE - sumim) * (ONE - sumjm);
      }
      
           
           
      /**
       * Function that calculates the derivative of self-relationship for 
       * node i along the search direction with respect to alpha
       * @param i      The current node index  
       * @param alpha  The current alpha value
       * @param h      The current search direction
       * @return       The derivative of the self-relationship for node i
       */
      inline double dsiia(int i, double alpha, double** h){
        double aux, aux1, p, sumh, sumim;
        p = sumh = sumim = ZERO;
        
        for (int m = 0; m < c; m++){
          aux1 = h[i][m];
          aux = U[i][m] + alpha * aux1;
          p += aux * aux1;
          sumh += aux1;
          sumim += aux;
        }
        sumim = ONE - sumim;
        return TWO * (p - sumh * sumim);
      }
      
      
      /**
       * Function that calculates the derivative of relationship between 
       * nodes i and j along the search direction with respect to alpha
       * @param i      The first node index  
       * @param j      The second node index  
       * @param alpha  The current alpha value
       * @param h      The current search direction
       * @return       The derivative of the relationship between nodes i and j
       */
      inline double dsija(int i, int j, double alpha, double** h){
        double p1, p2, sumhi, sumhj, sumim, sumjm, auxi, auxj, aux1i, aux1j;
        
        p1 = p2 = sumhi = sumhj = sumim = sumjm = ZERO;
        
        for (int m = 0; m < c; m++){
          auxi = h[i][m];
          auxj = h[j][m];
          aux1i = U[i][m] + alpha * auxi;
          aux1j = U[j][m] + alpha * auxj;
          p1 += aux1i * auxj;
          p2 += aux1j * auxi;
          sumhi += auxi;
          sumhj += auxj;
          sumim += aux1i;
          sumjm += aux1j;
        }
        return p1 + p2 - (ONE - sumim) * sumhj - (ONE - sumjm) * sumhi;
      }
      
          
      
      /**
       * Function that calculates the second derivative of self-relationship for 
       * node i along the search direction with respect to alpha
       * @param i    The current node index  
       * @param h    The current search direction
       * @return     The second derivative of self-relationship value for node i
       */
      inline double d2siia2(int i, double** h){
        double aux, aux1, sum, sum1;
        sum = sum1 = ZERO;
        
        for (int m = 0; m < c; m++){
          aux1 = h[i][m];
          aux = aux1 * aux1;
          sum += aux;
          sum1 += aux1;
        }
        return TWO * (sum + sum1 * sum1);
      }
      
     
      /**
       * Function that calculates the second derivative of the relationship 
       * between nodes i and j along the search direction with respect to alpha
       * @param i      The first node index  
       * @param j      The second node index  
       * @param h      The current search direction
       * @return       The second derivative of the relationship between nodes
       *               i and j
       */
      inline double d2sija2(int i, int j, double** h){
        double auxi, auxj, p, sumi, sumj;
        p = sumi = sumj = ZERO;
        
        for (int m = 0; m < c; m++){
          auxi = h[i][m];
          auxj = h[j][m];
          p += auxi * auxj;
          sumi += auxi;
          sumj += auxj;
        }
        return TWO * (p + sumi * sumj);
      }
      
   
      
      /**
       * Function that calculates the derivative of the error functional, D,  
       * with respect to the membership coefficient of node l in community k.
       * @param l The node for the derivative  
       * @param k The community for the derivative  
       * @return  The derivative of the function with respect to
       *          the membership coefficient of node l in community k
       */
      inline double dulk(int l, int k){
        return TWO * (t1(l, k) ) + FOUR * (t2(l, k) - t3(l, k));
      }
      
  
      
      /**
       * Function that calculates the first term of the derivative of the error 
       * functional, D, with respect to the membership coefficient of node 
       * l in community k.
       * @param l The node for the derivative  
       * @param k The community for the derivative  
       * @return  The first term of the derivative of the function with respect 
       *          to the membership coefficient of node l in community k
       */
      inline double t1(int l, int k){
        double sumlm;
        
        sumlm = ZERO;
        for (int m = 0; m < c; m++){
          sumlm += U[l][m];
        }
        return TWO * sii(l) * (U[l][k]- ONE + sumlm);
      }
      
      
      
      /**
       * Function that calculates the second term of the derivative of the error
       * functional, D, with respect to the membership coefficient of node l in 
       * community k.
       * @param l The node for the derivative  
       * @param k The community for the derivative  
       * @return  The second term of the derivative of the function with respect 
       *          to the membership coefficient of node l in community k
       */
      inline double t2(int l, int k){
        double t, sum;
        int j, m;
        
        t = ZERO;
        
        for (j = l + 1; j < N; j++){
          sum = ZERO;
          for (m = 0; m < c; m++){
            sum += U[j][m];
          }
          t += sij(l, j)*(U[j][k] - ONE + sum);
        }
        
        for (j = 0; j < l; j++){
          sum = ZERO;
          for (m = 0; m < c; m++){
            sum += U[j][m];
          }
          t += sij(j, l)*(U[j][k] - ONE + sum);
        }
        return t;
      }
      
            
     
      /**
       * Function that calculates the third term of the derivative of the  
       * error functional, D, with respect to the membership coefficient of node
       * l in community k.
       * @param l The node for the derivative  
       * @param k The community for the derivative  
       * @return  The third term of the derivative of the function with respect 
       *          to the membership coefficient of node l in community k
       */
      inline double t3(int l, int k){
        double t, sum;
        int j, m;
        
        t = ZERO;
        
        for (int i = 0; i < M; i++){  // Running on edges
          if (A->r(i) == l) {
            j = A->c(i);
            sum = ZERO;
            for (m = 0; m < c; m++){
              sum += U[j][m];
            }
            t += A->v(i) * (U[j][k] - ONE + sum);
          }
          
          if (A->c(i) == l){
            j = A->r(i);
            sum = ZERO;
            for (m = 0; m < c; m++){
              sum += U[j][m];
            }
            t += A->v(i) * (U[j][k] - ONE + sum);
          }
          
        }
        return t;
      }
      
  };

}    

#endif /* FUZZYCOMMUNITIES_H */

