This code generates filter-banks for partially designed multi-dimensional compactly supported Parseval framelets. The main function is funmin.m

# 1. funmin: The supporting mathematical theory can be found in ...

# 2. Specifically, if A is an s1 x s2 x ... x sn multi-dimensional low-pass filter matrix with positive elements adding up to 1, we first vectorize the matrix so as to create a 1 x numel(A) vector 'a'.

[example: the 3x3 matrix 1/16 * (1 2 1 ; 2 4 2 ; 1 2 1) becomes 
the 1x9 vector 1/16 * (1 2 1 2 4 2 1 2 1) by column vectorization.]

# 3. Given a collection of #L pre-designed s1 x s2 x ... x sn high-pass filter matrices (their coeffients must be adding up to 0), we vectorize each matrix and create an L x numel(A) matrix B1.

[example: the 3x3 matrix (-1 -2 -1 ; 0 0 0 ; 1 2 1) becomes 
the 1x9 vector (-1 0 1 -2 0 2 -1 0 1) by column vectorization.]

# 4. The rows of the matrix B1 are rescaled (optimized) to obtain a matrix newB1 so that the augmented matrix Q = (a ; newB1) has all its singular values between (0,1] but also as large as possible.

# 5. Finally, a matrix B2 is constructed such that " a' * a + newB1' * newB1 + B2' * B2 = diag(a) " and this guarantees that the collection of s1 x s2 x ... x sn high-pass filter coefficients in newB1 and B2 create a multi-dimensional filter bank defining a compactly supported Parseval framelet. 

The input arguments of the function funmin are the matrix A and the collection B1.
