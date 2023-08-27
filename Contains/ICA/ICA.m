function [EigenVals, A, s] = ICA(X, P)

  [U, L, ~] = svd(X*X');

  z = L(1:P, 1:P)^(-0.5) * U(:, 1:P)' * X;

  [~, zColumns] = size(z);


  C = z*z';
  C = (C+C')/2;

  for i = 1:1:24%zColumns-zRows

    CI = z(:,1:zColumns-i)*z(:,1+i:zColumns)';
    CI = (CI+CI')/2;

    C = [C CI];

  end

%   C1 = z(:,1:zColumns-25)*z(:,26:zColumns)';
%   C1 = (C1+C1')/2;
%   C2 = z(:,1:zColumns-50)*z(:,51:zColumns)';
%   C2 = (C2+C2')/2;
%   C3 = z(:,1:zColumns-75)*z(:,76:zColumns)';
%   C3 = (C3+C3')/2;
%   C4 = z(:,1:zColumns-100)*z(:,101:zColumns)';
%   C4 = (C4+C4')/2;
%   C5 = z(:,1:zColumns-125)*z(:,126:zColumns)';
%   C5 = (C5+C5')/2;
%   C6 = z(:,1:zColumns-150)*z(:,151:zColumns)';
%   C6 = (C6+C6')/2;
% 
%   C = [C C1 C2 C3 C4 C5 C6];



  [W, ~] = rjd(C);

  A = U(:, 1:P) * L(1:P, 1:P).^(0.5) * W;

  s = W' * z;

  EigenVals = diag(L(1:P, 1:P));


end
