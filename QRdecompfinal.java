public class QRdecompfinal {
    public static double dotProduct(double[] v1, double[] v2) {
        double product = 0;
        for (int i = 0; i < v1.length; i++) {
            product += v1[i] * v2[i];
        }
        return product;
    }
    public static double[][] matrixProduct(double[][] matrix1, double[][] matrix2) {
        double sum = 0;
        double[][] product = new double[matrix1.length][matrix1[1].length];
        for(int i = 0; i < matrix1.length; i++) {
            for(int j = 0; j < matrix1[0].length; j++) {
                product[i][j] = multiplyMatricesCell(matrix1,matrix2,i,j);
            }
        }
        return product;
    }

    public static double multiplyMatricesCell(double[][] matrix1, double[][] matrix2, int row, int col) {
        double cell = 0;
        for (int i = 0; i < matrix2.length; i++) {
            cell += matrix1[row][i] * matrix2[i][col];
        }
        return cell;
    }
 
    public static double[] projection(double[] v1, double[] v2) {
        double[] result = new double[v1.length];
        double scalar = dotProduct(v1, v2) / dotProduct(v2, v2);
        for (int i = 0; i < v2.length; i++) {
            result[i] = v2[i] * scalar;
        }
        return result;
    }
 
    public static void vectorSubtract(double[] v1, double[] v2) {
        for (int i = 0; i < v1.length; i++) {
            v1[i] -= v2[i];
        }
    }
 
    public static double[] copyColumn(double[][] matrix, int col) {
        double[] result = new double[matrix[0].length];
        for (int row = 0; row < result.length; row++) {
            result[row] = matrix[row][col];
        }
        return result;
    }

    public static double[] normalize(double[] v1) {
        double norm = 0;
        double[] result = v1;
        for(int i = 0; i < result.length; i++) {
            norm += result[i] * result[i];
        }
        norm = Math.sqrt(norm);
        for(int i = 0; i < result.length; i++) {
            result[i] /= norm;
        }
        return result;
    }
 
    public static double[][] getOrthoMatrix(double[][] matrix) {
        int m = matrix.length;
        int n = matrix[0].length;
        double[][] cols = new double[n][m];
        cols[0] = normalize(copyColumn(matrix, 0));
        for (int k = 1; k < m; k++) {
            // u_k = a_k - sum j=1..{k-1} proj_{u_j} a_k
            double[] a_k = copyColumn(matrix, k);
            double[] u_k = a_k.clone();
            for (int j = 0; j < k; j++) {
                // subtract proj_{u_j} a_k
                double[] proj = projection(a_k, cols[j]);
                vectorSubtract(u_k, proj);
            }
            // normalize u_k
            cols[k] = normalize(u_k);
        }
        // transpose `cols` and return
        double[][] result = new double[m][n];
        for (int row = 0; row < m; row++) {
            for (int col = 0; col < n; col++) {
                result[row][col] = cols[col][row];
            }
        }
        return result;
    }

    public static double[][] getUpperTriangularMatrix(double[][] matrix, double[][] orthoMatrix) {
        double[][] upperTriangularMatrix = new double[matrix.length][matrix[0].length];
        int zeroIndex = 0;
        int zeroTemp = 0;
        for(int i = 0; i < matrix.length; i++) {
            for(int j = 0; j < matrix[0].length; j++) {
                upperTriangularMatrix[i][j] = dotProduct(copyColumn(orthoMatrix, i),copyColumn(matrix, j));

            }
        }
        for(int i = 0; i < matrix.length; i++) {
            zeroTemp = zeroIndex;
            for(int j = 0; j < matrix[0].length; j++) {
                if(zeroTemp > 0) {
                    zeroTemp--;
                    upperTriangularMatrix[i][j] = 0.0;
                }
            }
            zeroIndex++;
        }
        return upperTriangularMatrix;
    }

    public static boolean isOrtho(double[][] matrix) {
        for(int i = 0; i < matrix.length; i++) {
            for(int j = i + 1; j < matrix.length; j++) {
                if(dotProduct(copyColumn(matrix,i),copyColumn(matrix,j)) > 0.001) {
                    System.out.println(dotProduct(copyColumn(matrix,i),copyColumn(matrix,j)));
                    return false;
                }
            }
        }
        return true;
    }
    public static double[][] schurForm(double[][] Q, double[][] R) {
        double[][] schurA = new double[Q.length][Q[1].length];
        while(1<2) {
            for(int i = 0; i < schurA.length; i++) {
                for(int j = 0; j < schurA[1].length; j++) {
                    schurA[i][j] = matrixProduct(R,Q)[i][j];
                }
            }
            for(int i = 0; i < Q.length; i++) {
                for(int j = 0; j < Q[1].length; j++) {
                    Q[i][j] = getOrthoMatrix(schurA)[i][j];
                }
            }
            for(int i = 0; i < R.length; i++) {
                for(int j = 0; j < R[1].length; j++) {
                    R[i][j] = getUpperTriangularMatrix(schurA,getOrthoMatrix(schurA))[i][j];
                }
            }
            for(int i = 0; i < schurA.length; i++) {
                for(int j = 0; j < schurA[0].length; j++) {
                    System.out.print(schurA[i][j] + " ");
                }
                System.out.println();
                try { 
                    Thread.sleep(0);
                } catch(Exception e) {
                    System.out.println(e);
                }
            }
            System.out.println(" ");
        }
        
    }

    public static void main(String[] args) {
        double[][] matrix = {
            {1,2,9},
            {12,11,2},
            {0,0,4}
        };
        double[][] orthoMatrix = getOrthoMatrix(matrix);
        System.out.println("Q");
        for(int i = 0; i < orthoMatrix.length; i++) {
            for(int j = 0; j < orthoMatrix[0].length; j++) {
                System.out.print(orthoMatrix[i][j] + " ");
            }
            System.out.println();
        }
        System.out.println();
        System.out.println("R");
        for(int i = 0; i < orthoMatrix.length; i++) {
            for(int j = 0; j < orthoMatrix[0].length; j++) {
                System.out.print(getUpperTriangularMatrix(matrix, orthoMatrix)[i][j] + " ");
            }
            System.out.println();
        }
        System.out.println();
        System.out.println("Q*R");
        for(int i = 0; i < orthoMatrix.length; i++) {
            for(int j = 0; j < orthoMatrix[0].length; j++) {
                System.out.print(matrixProduct(orthoMatrix, getUpperTriangularMatrix(matrix, orthoMatrix))[i][j] + " ");
            }
            System.out.println();
        }
        System.out.println();
        System.out.println(isOrtho(orthoMatrix));
        double[][] upperTriangleMatrix = getUpperTriangularMatrix(matrix, orthoMatrix);
        schurForm(orthoMatrix, upperTriangleMatrix);
    }
}