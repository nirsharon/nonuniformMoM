% script name: "test_vectorize_A_B"

% generate A, B, gamma
make_simple_A_B

% apply both directions
vvec     = A_B_to_VecAB(A, B);
[A1, B1] = VecAB_to_A_B(vvec, gamma);

% compare to origin
isequal(A,A1)
isequal(B,B1)