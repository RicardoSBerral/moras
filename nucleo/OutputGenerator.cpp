#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>

using namespace std;

/**
 * Determines if a double value is close to zero.
 */
bool isZero(double &value) {
  return (fabs(value) < 0.0001) ? true : false;
}

/**
 * Allocates memory for a matrix of LxN dimension.
 */
template<class T>
T **getMatrix(int L, int N) {
  T **matrix = new T*[L];
  
  for(int i = 0 ; i < L; ++i) {
    matrix[i] = new T[N];
  }
  
  return matrix;
}

/**
 * Delete memory allocated for a matrix
 */
template<class T>
void deleteMatrix(T **matrix, int L) {
  for(int i = 0; i < L; i++) {
    delete [] matrix[i];
  }
  
  delete [] matrix;
}

template<class T>
void fillMatrix(T **matrix, int L, T value, T diagonalValue) {
  for (int i = 0 ; i < L; ++i) {
    for (int j = 0; j< L; ++j) {
      matrix[i][j] = (i == j) ? diagonalValue : value;
     }
  }
}

template<class T>
void printMatrix(T **matrix, int L, int N) {
  for(int i = 0; i < L; ++i) {
    for(int j = 0; j < N; ++j) {
      cout << matrix[j][i] << " ";
    }
    cout << endl;
  }
  cout << endl;
}

/**
 * Returns a random number between [min, max] range
 */
float random_between(float min, float max) {
    double range = max - min;
    return min + (range*rand()/RAND_MAX);
}

/** 
 * This class generates random output to use with 'combining classifier' using
 * a probability p of success and a pairwise matrix Q.
 */
class OutputGenerator {
protected:
  // Number of classifiers
  int L;

  // Number of elements
  int N;

  // The probability of success for each classifier
  double *p;
  
  // The pairwise matrix
  double **Q;
  
  // Indicator of whether Q and p and owned by the object
  bool MemoryOwned;

  /**
   * Calculates the DISCR value
   */
  double getDiscr(double Qik, double pi, double pk) {
    double preDiscr = 1 - Qik + 2 * Qik * (pi - pk);
    return (preDiscr * preDiscr)
           - 8 * Qik * (1 - pi) * pk * (Qik - 1);
  }
  
  /**
   * Calculates P1 and P2 values from Q and p
   */
  void calculate(double **P1, double **P2) {
    for(int i = 0; i < L; ++i) {
      for(int k = 0; k < L; ++k) {
        if (i != k) {
          if (isZero(Q[i][k])) {
            P2[i][k] = p[k];            
            P1[i][k] = 1 - p[k];
          } else {
              double Discr = getDiscr(Q[i][k], p[i], p[k]);
              P2[i][k] = (-(1 - Q[i][k] + 2*Q[i][k] * (p[i] - p[k])) + sqrt(Discr)) / (4*Q[i][k]*(1 - p[i]));
              P1[i][k] = 1 - P2[i][k] - p[k]/p[i] + P2[i][k]/p[i];            
          }
        } else {
          // These values will not be used.
          P2[i][k] = 0;
          P1[i][k] = 0;
        }
      }
    }
  }
  
  /**
   * Permutes an array of L elements
   */
  void permutate(int *order, int L) {
    for(int i = 0; i < L - 1; ++i) {
      int j = lroundf(random_between(i, L - 1));
      
      int swap = order[i];
      order[i] = order[j];
      order[j] = swap;
      //cout << order[i] << " ";
    }
    //cout << order[L - 1] << " ";
    //cout << endl;
  }
  
public:
  // Constructor
  OutputGenerator(int L, int N, double *p, double **Q) : L(L), N(N), p(p), Q(Q)  {
    MemoryOwned = false;
  }
  OutputGenerator(int L, int N, double p, double Q) : L(L), N(N)  {
    this->p = new double[L];
    this->Q = getMatrix<double>(L, L);

    for(int j = 0; j < L; ++j) {
      this->p[j] = p;
    }

    // fill Q matrix with a fixed value
    fillMatrix(this->Q, L, Q, 1.0);    

    MemoryOwned = true;
  }

  ~OutputGenerator()  {
    if (MemoryOwned) {
      deleteMatrix(Q, L);
      delete [] p;
    }
  }
  
  int **generate() {
    int **outputs = getMatrix<int>(L, N);
    return generate(outputs);
  }
  
  /**
   * Generates classifier outputs of fixed accuracy and diversity
   */
  int **generate(int **outputs) {
    double **P1 = getMatrix<double>(L, L);  // The P1 probabilities (probability of a 1 turns to 0) 
    double **P2 = getMatrix<double>(L, L);  // The P2 probabilities (probability of a 0 turns to 1)
    int *order;                // Followed order when generating the classifier outputs
    //int outputs[L];              // Outputs for one element and for each classifier
    
    order = new int[L];

    // Initialize order array
    for(int t = 0; t < L; ++t) 
      order[t] = t;
          
    // Calculate P1 and P2 values from Q and p
    calculate(P1, P2);
    //cout << endl << "P1: " << P1[0][1];
    //cout << endl << "P2: " << P2[0][1] << endl;
    
    // For each element
    for (int j = 0; j < N; ++j) {
      // Permutate generating order
      permutate((int*)order, L);
      // If a random number is greater than probability p of first element in order array, 
      // set output to 1
      outputs[order[0]][j] = (random_between(0, 1) <= p[order[0]]) ? 1 : 0;
            
      // For each remaining classifier
      for (int t = 1; t < L; ++t) {
        int prevIndex = order[t-1];
        int currIndex = order[t];
        double currentP1 = P1[prevIndex][currIndex];
        double currentP2 = P2[prevIndex][currIndex];
        
        //cout << endl << "P1: " << currentP1 << endl;
        //cout << endl << "P2: " << currentP2 << endl;
        
        double changingProbability = (outputs[prevIndex][j] == 1) ? currentP1 : currentP2;
        
        // If generated number is greather than probability, output does not change.
        outputs[currIndex][j] = (random_between(0, 1) <= changingProbability) ? !outputs[prevIndex][j] : outputs[prevIndex][j];
      }      
    }
    
    deleteMatrix(P1, L);
    deleteMatrix(P2, L);
    delete [] order;

    return outputs;
  }
  
  /**
   * Calculates the accuracy of a classifier from its
   *  outputs vector
   */
  static double calculateAccuracy(int *output, int N) {
    double accuracy = 0.0;
    
    for(int i = 0; i < N; ++i) {
      if (output[i] == 1) ++accuracy;
    }
    
    return accuracy / N;
  }
  
  /**
   * Calculates the accuracy of a classifier from its
   *  outputs vector
   */
  static double calculateAccuracy(int **output, int N, int M) {
    double accuracy = 0.0;
        int *o;

    o = new int[N];

    for(int i = 0; i < N; ++i) {
      o[i] = 0;
      for(int j = 0; j < M; ++j) {
         o[i] += output[j][i];
      }
      o[i] = o[i] > M/2 ? 1 : 0;
    }

    accuracy = calculateAccuracy(o, N);
    delete [] o;

    return accuracy;
  }
  
  
  /**
   * Realizes test 1 described in Kuncheva paper:
   * - Tests with probabilities and pairwises fixed
   * - L = 3 classifiers
   * - N = 200 samples
   */ 
  static void test1(int L = 3, int N = 200) {
    static int P_LENGTH = 4;
    static double P_ARRAY[] = {0.6, 0.7, 0.8, 0.9};
    //static double P_ARRAY[] = {0.6};
    static int Q_LENGTH = 2;//21;
    static double Q_ARRAY[] = {-1.0, 0.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 
                             0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

    int N00, N01, N10, N11;
    double Qexp;
    double *p = new double[L];
    double **Q = getMatrix<double>(L, L);
    int **outputs = getMatrix<int>(L, N);
    double **probabilities = getMatrix<double>(Q_LENGTH , L);
        
    // For each experiment...
    for(int i = 0; i < P_LENGTH; ++i) {
      // Fill p vector with a fixed`value
      for(int j = 0; j < L; ++j) {
        p[j] = P_ARRAY[i];
      }
      
      // For each Q value in experiment...
      for(int j = 0; j < Q_LENGTH; ++j) {
        // fill Q matrix with a fixed value
        fillMatrix(Q, L, Q_ARRAY[j], 1.0);
        
        OutputGenerator o(L, N, (double*)p, Q);
        o.generate(outputs);
        
        printMatrix(outputs, N, L);
        
        // Calculate accuracy for each classifier
        for(int z = 0; z < L; ++z) {
          probabilities[j][z] = calculateAccuracy(outputs[z], N);
          cout << probabilities[j][z] << '\t';
          cout << Q_ARRAY[j] << '\t' << p[z] -  probabilities[j][z] << endl;
        }
        // Calculate the obtained average Q
        Qexp = 0.0;
        for(int ic = 0; ic < L; ++ic) {
          for(int jc = ic + 1; jc < L; ++jc) {
            N00 = N01 = N10 = N11 = 0; 
            for(int io = 0; io < N; ++io) {
              if (outputs[ic][io]==1 && outputs[jc][io]==1) N11++;
              else if (outputs[ic][io]==0 && outputs[jc][io]==1) N01++;
              else if (outputs[ic][io]==1 && outputs[jc][io]==0) N10++;
              else if (outputs[ic][io]==0 && outputs[jc][io]==0) N00++;
            }
            double N11_N00 = (double)N11*N00;
            double N01_N10 = (double)N01*N10;
            Qexp += (N11_N00-N01_N10)/(N11_N00+N01_N10);
          }
        }
        Qexp = Qexp/((L*(L-1))/2);
        cout << Q_ARRAY[j] << '\t' << Qexp << '\t' << calculateAccuracy(outputs, N, L) << "<-----------------------------------" << endl;
        //break;
      }
      cout << "----------------" << endl;
      //break;
    }
    
    deleteMatrix(probabilities, Q_LENGTH);
    deleteMatrix(Q, L);
    deleteMatrix(outputs, L);
    delete [] p;
  }
  
  /**
   * Generate a random Q matrix and p vector and test the 
   * OutputGenerator.generate() method
    */
  static void randomTest() {
    int L = 3;
    int N = 15;
    double *p = new double[L];
    double **Q = getMatrix<double>(L, L);

    for(int i = 0; i < L; ++i) {
      p[i] = random_between(0.01, 0.99);
      printf("%f ", p[i]);
      for(int z = i; z < L; ++z) {
        Q[i][z] = (i == z) ? 1 : random_between(0.1, 0.9);
        Q[z][i] = Q[i][z];
      }
    }

    cout << endl << endl;

    for(int i = 0; i < L; ++i) {
      for(int z = 0; z < L; ++z) {
        printf("%f ", Q[i][z]);
      }
      cout << endl;
    }
    cout << endl;

    OutputGenerator o(L, N, (double*)p, Q);

    int **outputs = o.generate();
    
    printMatrix(outputs, N, L);

    deleteMatrix(Q, L);
    deleteMatrix(outputs, L);
    delete [] p;
  }
};

/*int main(int argn, char **argv) {
  if (argn != 4) {
    cout << argv[0] << " [L] [N] [repeticiones]" << endl;
    return 1;
  }

  srand(time(NULL));
  int limit = atoi(argv[3]);
  for(int i = 0; i < limit; ++i) {  
    //OutputGenerator::randomTest();
    OutputGenerator::test1(atoi(argv[1]), atoi(argv[2]));
  }
  return 0;
}*/
