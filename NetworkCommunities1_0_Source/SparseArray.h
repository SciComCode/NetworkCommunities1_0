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
* Release date: August, 2018
*
* Purpose: Definition (header) and implementation file for class SparseArray
*
*****************************************************************************/

/**
 * @file
 * Definition (header) and implementation file for class SparseArray.
 * This class is part of the package networkCommunities.
*/

#ifndef SPARSEARRAY_H
#define SPARSEARRAY_H


using namespace std;

namespace NetworkCommunities {
  /**
 * \brief Defines a sparse array
 *
 * This class represents a sparse array where only non-zero elements are 
 * stored.
 * @author SciCom research group-E.S.I. Universidad de Castilla-La Mancha.
 *         Paseo de la Universidad 4, 13004-Ciudad Real. SPAIN
 * 
 * @date August, 2018
 */
  

  class SparseArray {
    private:
      /**
       * Structure element. This stores the value of the row, colum and 
       * value given to an item. In the current case each element represent an 
       * edge of the network. The row field stores the origin node of the edge,
       * the colum field the destination node and the value field stores the
       * weight of the edge.
       */
      struct element{
        int row;
        int column;
        double value;
      };
      
      element* matrix; // Array of elements
      int M;           // Number of edges added to the network
      double sum;      // Sum of squared values (weights) added to the array
      int index;       // Current index of array matrix

    public:
      
      /**
       * Constructor
       * @param n The number of rows (edges) in the array
       */
      SparseArray(int n){
        M = n; 
        sum = 0.0;
        index = -1;
        matrix = new element[M]; 
      }
      
     
      /**
       * Method for adding the i, j element to the sparse array. 
       * @param i The row of the element
       * @param j The column of the element
       * @param value The value of the (i, j) element
       */
      void put (int i, int j, double value){
        index++;
        matrix[index].row = i;
        matrix[index].column = j;
        matrix[index].value = value;
        sum += value * value;
      }
      
      
      /**
       * Method for retrieving the row of the ith element from the sparse array
       * @param i The edge index of the element
       * @return The row of the ith element (edge)
       */
      int r (int i){
        return matrix[i].row;
      }
      
      
      /**
       * Method for retrieving the column of the ith element from the 
       * sparse array
       * @param i The edge index of the element
       * @return The column of the ith element (edge)
       */
      int c (int i){
        return matrix[i].column;
      }
      

      /**
       * Method for retrieving the value of the ith element from the 
       * sparse array
       * @param i The edge index of the element
       * @return The value (weight) of the ith element (edge)
       */
      double v (int i){
        return matrix[i].value;
      }      
      
      
      /**
       * Method for getting the number of edges (non-zero elements) in the 
       * network
       * @return The number of edges in the network.
       */
      int get_nEdges(){
        return M;
      }
      
      
      /**
       * Method for getting the sum of squared weights of non-zero elements in
       * the network
       * @return The sum of squared weights of non-zero elements.
       */
      double get_weightsSum(){
        return sum;
      }
      
            
      /**
       * Destructor method
       */
      ~SparseArray(){
        delete[] matrix;
      }
  };
  
}

#endif /* SPARSEARRAY_H */

