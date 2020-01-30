function A = setupMatrix(u, prod)
    A = SparseTensor();
    A = A.setFromTensorProd(u, prod);
    A = A.getMatrix();
    A = full(A);
end
